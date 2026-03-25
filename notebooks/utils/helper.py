"""Hilfsfunktionen zum Lesen und Vorverarbeiten von BRNS-Ausgabedateien (.dat)."""

from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np

ArrayLike = np.ndarray
DatasetMap = Dict[str, Dict[str, ArrayLike]]


def load_dat(path: Path) -> Optional[ArrayLike]:
    """Laedt eine .dat-Datei als 2D-Array.

    Bei 1D-Daten wird auf Form (n, 1) normalisiert.
    """
    try:
        arr = np.loadtxt(path)
        return arr.reshape(-1, 1) if arr.ndim == 1 else arr
    except Exception as exc:
        print('  Ladefehler', path.name, exc)
        return None


def collect_common_dat_datasets(ref_dir: Path, py_dir: Path) -> Tuple[DatasetMap, List[str], List[str], List[str]]:
    """Sammelt gemeinsame .dat-Dateien aus Referenz- und Python-Ordner.

    Returns:
        datasets: Mapping Dateiname -> {'ref': ndarray, 'py': ndarray}
        common: sortierte gemeinsame Dateinamen
        only_ref: nur in Referenz vorhandene Dateinamen
        only_py: nur in Python vorhandene Dateinamen
    """
    ref_names = {f.name for f in ref_dir.glob('*.dat')}
    py_names = {f.name for f in py_dir.glob('*.dat')}

    common = sorted(ref_names & py_names)
    only_ref = sorted(ref_names - py_names)
    only_py = sorted(py_names - ref_names)

    datasets: DatasetMap = {}
    for fname in common:
        ref = load_dat(ref_dir / fname)
        py = load_dat(py_dir / fname)
        if ref is not None and py is not None:
            datasets[fname] = {'ref': ref, 'py': py}

    return datasets, common, only_ref, only_py


def collect_dat_datasets(result_dir: Path) -> Dict[str, ArrayLike]:
    """Sammelt alle .dat-Dateien aus einem Ergebnisordner.

    Args:
        result_dir: Pfad zum Ergebnisordner mit .dat-Dateien.

    Returns:
        datasets: Mapping Dateiname -> ndarray
    """
    datasets: Dict[str, ArrayLike] = {}
    
    if not result_dir.exists():
        print(f'Warnung: Ordner nicht gefunden: {result_dir}')
        return datasets
    
    for fpath in sorted(result_dir.glob('*.dat')):
        arr = load_dat(fpath)
        if arr is not None:
            datasets[fpath.name] = arr
    
    return datasets


def build_result_dirs(example: str, build_root: Optional[Path] = None) -> Tuple[Path, Path]:
    """Ermittelt Referenz- und Python-Ergebnisordner fuer ein Beispiel.

    Args:
        example: Beispielname, z.B. 'equilibrium', 'single_species', 'multiple_species'.
        build_root: Optionales Build-Root. Standard ist BRNSPackage/build relativ zum Repo.

    Returns:
        (ref_dir, py_dir)
    """
    if build_root is None:
        build_root = Path(__file__).resolve().parents[2] / 'build'

    base = build_root / example / 'results'
    ref_dir = base / 'reference'
    py_dir = base / 'python'
    return ref_dir, py_dir


def prepare_example_datasets(
    example: str,
    build_root: Optional[Path] = None,
    verbose: bool = True,
) -> DatasetMap:
    """Laedt und bereitet gemeinsame .dat-Datasets fuer ein Beispiel vor.

    Diese Funktion kapselt Pfadaufloesung, Dateiabgleich und Laden der Daten,
    damit das Notebook sich auf die Visualisierung konzentrieren kann.
    """
    ref_dir, py_dir = build_result_dirs(example, build_root=build_root)
    datasets, common, only_ref, only_py = collect_common_dat_datasets(ref_dir, py_dir)

    if verbose:
        print('Beispiel :', example)
        print('Referenz :', ref_dir)
        print('Python   :', py_dir)
        print('Ref OK   :', ref_dir.exists())
        print('Py  OK   :', py_dir.exists())
        print('Gemeinsame .dat-Dateien:', len(common))
        if only_ref:
            print('Nur Referenz:', only_ref)
        if only_py:
            print('Nur Python  :', only_py)
        print('Geladen:', len(datasets), 'Dateien')

    return datasets


def split_snapshots(arr: ArrayLike) -> List[ArrayLike]:
    """Teilt ein Array in Zeit-Snapshots anhand eines Resets der Tiefenachse."""
    if arr is None or len(arr) == 0:
        return []

    if arr.shape[1] >= 2:
        x = arr[:, 1]  # Tiefe
    else:
        x = np.arange(len(arr), dtype=float)

    # Snapshot-Grenzen: wenn x wieder "nach unten springt" (z.B. von 50 auf 0)
    reset_idx = np.where(np.diff(x) < 0)[0] + 1
    starts = np.r_[0, reset_idx]
    ends = np.r_[reset_idx, len(arr)]

    return [arr[s:e] for s, e in zip(starts, ends) if e > s]
