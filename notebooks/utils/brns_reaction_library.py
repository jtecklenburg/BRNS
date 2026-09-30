from __future__ import annotations

import ctypes
import json
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np


class BRNSReactionLibrary:
    """Thin ctypes wrapper around a model-specific BRNS reaction shared library."""

    def __init__(self, library_path: str | Path, metadata_path: str | Path | None = None):
        self.library_path = Path(library_path).resolve()
        if not self.library_path.exists():
            raise FileNotFoundError(f"Shared library not found: {self.library_path}")

        self.metadata_path = Path(metadata_path).resolve() if metadata_path else None
        self.metadata = self._load_metadata(self.metadata_path)

        self.lib = ctypes.CDLL(str(self.library_path))
        self._configure_signatures()

        self.n_species = self._get_count(self.lib.brns_get_n_species)
        self.n_reactions = self._get_count(self.lib.brns_get_n_reactions)
        self.n_parameters = self._get_count(self.lib.brns_get_n_parameters)

        self.species_names = list(self.metadata.get("species_names", []))
        self.parameter_names = list(self.metadata.get("parameter_names", []))
        self.default_parameter_values = np.asarray(
            self.metadata.get("default_parameter_values", np.zeros(self.n_parameters)),
            dtype=np.float64,
        )
        if self.default_parameter_values.size == 0 and self.n_parameters > 0:
            self.default_parameter_values = np.zeros(self.n_parameters, dtype=np.float64)
        if self.default_parameter_values.size != self.n_parameters:
            raise ValueError(
                f"Parameter vector length mismatch: metadata has {self.default_parameter_values.size}, "
                f"library expects {self.n_parameters}"
            )

    @classmethod
    def from_build_dir(cls, build_dir: str | Path) -> "BRNSReactionLibrary":
        requested_dir = Path(build_dir).resolve()

        def try_load_from_dir(candidate_dir: Path) -> "BRNSReactionLibrary" | None:
            metadata_path = candidate_dir / "brns_library_metadata.json"
            if metadata_path.exists():
                metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
                library_path = Path(
                    metadata.get("library_path", candidate_dir / "lib" / metadata.get("library_basename", ""))
                )
                if not library_path.is_absolute():
                    library_path = (candidate_dir / library_path).resolve()
                if library_path.exists():
                    return cls(library_path, metadata_path=metadata_path)

            candidates = sorted((candidate_dir / "lib").glob("libbrns_reactions_*.so"))
            if candidates:
                return cls(candidates[0], metadata_path=metadata_path if metadata_path.exists() else None)
            return None

        direct_hit = try_load_from_dir(requested_dir)
        if direct_hit is not None:
            return direct_hit

        build_root = requested_dir.parent
        stem = requested_dir.name
        variants = {stem, stem.rstrip("_"), f"{stem.rstrip('_')}_", stem.replace("-", "_"), stem.replace("_", "-")}
        variants = {v for v in variants if v}

        candidate_dirs: list[Path] = []
        if build_root.exists():
            for variant in variants:
                variant_dir = build_root / variant
                if variant_dir != requested_dir and variant_dir.exists():
                    candidate_dirs.append(variant_dir)

            for sibling in sorted(build_root.iterdir()):
                if not sibling.is_dir() or sibling == requested_dir or sibling in candidate_dirs:
                    continue
                if sibling.name.startswith(stem) or stem.startswith(sibling.name):
                    candidate_dirs.append(sibling)

        for candidate_dir in candidate_dirs:
            loaded = try_load_from_dir(candidate_dir)
            if loaded is not None:
                return loaded

        searched_dirs = [requested_dir] + candidate_dirs
        searched_str = "\n".join(f"- {path}" for path in searched_dirs)
        raise FileNotFoundError(
            "No BRNS reaction library found. Searched directories:\n"
            f"{searched_str}\n"
            "Expected either brns_library_metadata.json or lib/libbrns_reactions_*.so"
        )

    def _load_metadata(self, metadata_path: Path | None) -> dict:
        if metadata_path is None or not metadata_path.exists():
            return {}
        return json.loads(metadata_path.read_text(encoding="utf-8"))

    def _configure_signatures(self) -> None:
        int_ptr = ctypes.POINTER(ctypes.c_int)
        dbl_ptr = ctypes.POINTER(ctypes.c_double)

        self.lib.brns_get_n_species.argtypes = [int_ptr]
        self.lib.brns_get_n_species.restype = None
        self.lib.brns_get_n_reactions.argtypes = [int_ptr]
        self.lib.brns_get_n_reactions.restype = None
        self.lib.brns_get_n_parameters.argtypes = [int_ptr]
        self.lib.brns_get_n_parameters.restype = None

        self.lib.brns_react_cell.argtypes = [
            dbl_ptr,
            dbl_ptr,
            dbl_ptr,
            ctypes.c_int,
            ctypes.c_double,
            int_ptr,
            int_ptr,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            dbl_ptr,
        ]
        self.lib.brns_react_cell.restype = None

    def _get_count(self, func) -> int:
        value = ctypes.c_int()
        func(ctypes.byref(value))
        return int(value.value)

    def _as_species_vector(self, values: Sequence[float], *, name: str) -> np.ndarray:
        arr = np.ascontiguousarray(values, dtype=np.float64)
        if arr.shape != (self.n_species,):
            raise ValueError(f"{name} must have shape ({self.n_species},), got {arr.shape}")
        return arr

    def _as_parameter_vector(self, parameter_vector: Sequence[float] | None) -> np.ndarray:
        if self.n_parameters == 0:
            return np.ascontiguousarray(np.zeros(1), dtype=np.float64)
        if parameter_vector is None:
            parameter_vector = self.default_parameter_values
        arr = np.ascontiguousarray(parameter_vector, dtype=np.float64)
        if arr.shape != (self.n_parameters,):
            raise ValueError(
                f"parameter_vector must have shape ({self.n_parameters},), got {arr.shape}"
            )
        return arr

    def _as_fixed_boundary_vector(self, fixed_boundary: Sequence[int] | None) -> np.ndarray:
        fixed_arr = np.ascontiguousarray(
            np.zeros(self.n_species, dtype=np.int32) if fixed_boundary is None else fixed_boundary,
            dtype=np.int32,
        )
        if fixed_arr.shape != (self.n_species,):
            raise ValueError(f"fixed_boundary must have shape ({self.n_species},), got {fixed_arr.shape}")
        return fixed_arr

    def react_cell(
        self,
        after_transport: Sequence[float],
        *,
        before_transport: Sequence[float] | None = None,
        dt: float,
        fixed_boundary: Sequence[int] | None = None,
        position: tuple[float, float, float] = (0.0, 0.0, 0.0),
        porosity: float = 1.0,
        water_saturation: float = 1.0,
        parameter_vector: Sequence[float] | None = None,
        max_substeps: int = 32,
        retry_return_codes: tuple[int, ...] = (2, 3),
    ) -> np.ndarray:
        after_arr = self._as_species_vector(after_transport, name="after_transport")
        before_arr = self._as_species_vector(
            before_transport if before_transport is not None else after_transport,
            name="before_transport",
        )
        output_arr = np.zeros(self.n_species, dtype=np.float64)
        fixed_arr = self._as_fixed_boundary_vector(fixed_boundary)
        parameter_arr = self._as_parameter_vector(parameter_vector)

        output_arr, return_code_value = self._react_cell_once(
            after_arr=after_arr,
            before_arr=before_arr,
            dt=float(dt),
            fixed_arr=fixed_arr,
            position=position,
            porosity=float(porosity),
            water_saturation=float(water_saturation),
            parameter_arr=parameter_arr,
        )

        if return_code_value == 0:
            return output_arr

        if return_code_value in retry_return_codes and max_substeps > 1:
            current = after_arr.copy()
            n_substeps = 2
            while n_substeps <= max_substeps:
                dt_sub = float(dt) / float(n_substeps)
                success = True
                for _ in range(n_substeps):
                    current, code = self._react_cell_once(
                        after_arr=current,
                        before_arr=current,
                        dt=dt_sub,
                        fixed_arr=fixed_arr,
                        position=position,
                        porosity=float(porosity),
                        water_saturation=float(water_saturation),
                        parameter_arr=parameter_arr,
                    )
                    if code != 0:
                        success = False
                        break
                if success:
                    return current
                n_substeps *= 2

        raise RuntimeError(f"brns_react_cell failed with return code {return_code_value}")

    def _react_cell_with_buffers(
        self,
        *,
        after_arr: np.ndarray,
        before_arr: np.ndarray,
        dt: float,
        fixed_arr: np.ndarray,
        position: tuple[float, float, float],
        porosity: float,
        water_saturation: float,
        parameter_arr: np.ndarray,
        max_substeps: int,
        retry_return_codes: tuple[int, ...],
        out_arr: np.ndarray,
        work_current: np.ndarray,
        work_trial: np.ndarray,
    ) -> int:
        return_code_value = self._react_cell_once_into(
            after_arr=after_arr,
            before_arr=before_arr,
            dt=dt,
            fixed_arr=fixed_arr,
            position=position,
            porosity=porosity,
            water_saturation=water_saturation,
            parameter_arr=parameter_arr,
            output_arr=out_arr,
        )

        if return_code_value == 0:
            return 0

        if return_code_value in retry_return_codes and max_substeps > 1:
            work_current[:] = after_arr
            n_substeps = 2
            while n_substeps <= max_substeps:
                dt_sub = float(dt) / float(n_substeps)
                success = True
                current = work_current
                trial = work_trial
                for _ in range(n_substeps):
                    code = self._react_cell_once_into(
                        after_arr=current,
                        before_arr=current,
                        dt=dt_sub,
                        fixed_arr=fixed_arr,
                        position=position,
                        porosity=porosity,
                        water_saturation=water_saturation,
                        parameter_arr=parameter_arr,
                        output_arr=trial,
                    )
                    if code != 0:
                        success = False
                        break
                    current, trial = trial, current
                if success:
                    out_arr[:] = current
                    return 0
                n_substeps *= 2

        return int(return_code_value)

    def _react_cell_once(
        self,
        *,
        after_arr: np.ndarray,
        before_arr: np.ndarray,
        dt: float,
        fixed_arr: np.ndarray,
        position: tuple[float, float, float],
        porosity: float,
        water_saturation: float,
        parameter_arr: np.ndarray,
    ) -> tuple[np.ndarray, int]:
        output_arr = np.zeros(self.n_species, dtype=np.float64)
        return output_arr, self._react_cell_once_into(
            after_arr=after_arr,
            before_arr=before_arr,
            dt=dt,
            fixed_arr=fixed_arr,
            position=position,
            porosity=porosity,
            water_saturation=water_saturation,
            parameter_arr=parameter_arr,
            output_arr=output_arr,
        )

    def _react_cell_once_into(
        self,
        *,
        after_arr: np.ndarray,
        before_arr: np.ndarray,
        dt: float,
        fixed_arr: np.ndarray,
        position: tuple[float, float, float],
        porosity: float,
        water_saturation: float,
        parameter_arr: np.ndarray,
        output_arr: np.ndarray,
    ) -> int:
        return_code = ctypes.c_int(0)
        self.lib.brns_react_cell(
            after_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            before_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            output_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.c_int(self.n_species),
            ctypes.c_double(dt),
            fixed_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            ctypes.byref(return_code),
            ctypes.c_double(float(position[0])),
            ctypes.c_double(float(position[1])),
            ctypes.c_double(float(position[2])),
            ctypes.c_double(porosity),
            ctypes.c_double(water_saturation),
            parameter_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        )
        return int(return_code.value)

    def react_batch_inplace(
        self,
        concentrations: np.ndarray,
        *,
        previous_state: np.ndarray | None = None,
        dt: float,
        active_mask: Iterable[bool] | None = None,
        positions: np.ndarray | None = None,
        porosity: float = 1.0,
        water_saturation: float = 1.0,
        parameter_vector: Sequence[float] | None = None,
        max_substeps_per_cell: int = 1,
        retry_return_codes: tuple[int, ...] = (2, 3),
        on_cell_failure: str = "keep_transport",
        verbose_failures: bool = True,
        out: np.ndarray | None = None,
        fixed_boundary: np.ndarray | None = None,
        work_current: np.ndarray | None = None,
        work_trial: np.ndarray | None = None,
    ) -> np.ndarray:
        conc = np.asarray(concentrations, dtype=np.float64)
        if conc.ndim != 2 or conc.shape[1] != self.n_species:
            raise ValueError(
                f"concentrations must have shape (n_cells, {self.n_species}), got {conc.shape}"
            )
        if not conc.flags.c_contiguous:
            raise ValueError("concentrations must be C-contiguous for react_batch_inplace")

        if previous_state is None:
            prev = conc
        else:
            prev = np.asarray(previous_state, dtype=np.float64)
            if prev.shape != conc.shape:
                raise ValueError(f"previous_state must match concentrations shape {conc.shape}, got {prev.shape}")
            if not prev.flags.c_contiguous:
                raise ValueError("previous_state must be C-contiguous for react_batch_inplace")

        if active_mask is None:
            mask = np.ones(conc.shape[0], dtype=bool)
        else:
            mask = np.asarray(list(active_mask), dtype=bool)
            if mask.shape != (conc.shape[0],):
                raise ValueError(f"active_mask must have shape ({conc.shape[0]},), got {mask.shape}")

        if positions is None:
            pos = np.zeros((conc.shape[0], 3), dtype=np.float64)
        else:
            pos = np.ascontiguousarray(positions, dtype=np.float64)
            if pos.shape != (conc.shape[0], 3):
                raise ValueError(f"positions must have shape ({conc.shape[0]}, 3), got {pos.shape}")

        if out is None:
            out_arr = conc
        else:
            out_arr = np.asarray(out, dtype=np.float64)
            if out_arr.shape != conc.shape:
                raise ValueError(f"out must have shape {conc.shape}, got {out_arr.shape}")
            if not out_arr.flags.c_contiguous:
                raise ValueError("out must be C-contiguous for react_batch_inplace")
            if out_arr is not conc:
                out_arr[:] = conc

        parameter_arr = self._as_parameter_vector(parameter_vector)
        fixed_arr = self._as_fixed_boundary_vector(fixed_boundary)

        work_current_arr = np.empty(self.n_species, dtype=np.float64) if work_current is None else np.asarray(work_current, dtype=np.float64)
        work_trial_arr = np.empty(self.n_species, dtype=np.float64) if work_trial is None else np.asarray(work_trial, dtype=np.float64)
        if work_current_arr.shape != (self.n_species,):
            raise ValueError(f"work_current must have shape ({self.n_species},), got {work_current_arr.shape}")
        if work_trial_arr.shape != (self.n_species,):
            raise ValueError(f"work_trial must have shape ({self.n_species},), got {work_trial_arr.shape}")
        if not work_current_arr.flags.c_contiguous:
            raise ValueError("work_current must be C-contiguous for react_batch_inplace")
        if not work_trial_arr.flags.c_contiguous:
            raise ValueError("work_trial must be C-contiguous for react_batch_inplace")

        failure_count = 0
        for cell_index in range(conc.shape[0]):
            if not mask[cell_index]:
                out_arr[cell_index] = conc[cell_index]
                continue

            code = self._react_cell_with_buffers(
                after_arr=conc[cell_index],
                before_arr=prev[cell_index],
                dt=dt,
                fixed_arr=fixed_arr,
                position=tuple(pos[cell_index]),
                porosity=float(porosity),
                water_saturation=float(water_saturation),
                parameter_arr=parameter_arr,
                max_substeps=max_substeps_per_cell,
                retry_return_codes=retry_return_codes,
                out_arr=out_arr[cell_index],
                work_current=work_current_arr,
                work_trial=work_trial_arr,
            )
            if code == 0:
                continue

            failure_count += 1
            if on_cell_failure == "keep_previous":
                out_arr[cell_index] = prev[cell_index]
            elif on_cell_failure == "raise":
                raise RuntimeError(f"brns_react_cell failed with return code {code}")
            else:
                out_arr[cell_index] = conc[cell_index]

            if verbose_failures and failure_count <= 5:
                print(
                    f"Warning: chemistry failed at cell {cell_index} (dt={dt}). "
                    f"Applied fallback='{on_cell_failure}'. Details: brns_react_cell failed with return code {code}"
                )

        if verbose_failures and failure_count > 5:
            print(f"Warning: chemistry failures encountered in {failure_count} cells (only first 5 shown).")
        return out_arr

    def react_batch(
        self,
        concentrations: np.ndarray,
        *,
        previous_state: np.ndarray | None = None,
        dt: float,
        active_mask: Iterable[bool] | None = None,
        positions: np.ndarray | None = None,
        porosity: float = 1.0,
        water_saturation: float = 1.0,
        parameter_vector: Sequence[float] | None = None,
        max_substeps_per_cell: int = 1,
        retry_return_codes: tuple[int, ...] = (2, 3),
        on_cell_failure: str = "keep_transport",
        verbose_failures: bool = True,
    ) -> np.ndarray:
        conc = np.ascontiguousarray(concentrations, dtype=np.float64)
        prev = None if previous_state is None else np.ascontiguousarray(previous_state, dtype=np.float64)
        out = np.empty_like(conc)
        return self.react_batch_inplace(
            conc,
            previous_state=prev,
            dt=dt,
            active_mask=active_mask,
            positions=positions,
            porosity=porosity,
            water_saturation=water_saturation,
            parameter_vector=parameter_vector,
            max_substeps_per_cell=max_substeps_per_cell,
            retry_return_codes=retry_return_codes,
            on_cell_failure=on_cell_failure,
            verbose_failures=verbose_failures,
            out=out,
        )
