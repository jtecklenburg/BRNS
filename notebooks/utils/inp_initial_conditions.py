from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

ArrayLike = np.ndarray
InpDatasetMap = Dict[str, ArrayLike]


def load_inp_file(path: Path) -> ArrayLike:
    """Load a single `.inp` file as a 2D array with two columns.

    Expected file format per row:
    - column 0: concentration/value
    - column 1: depth
    """
    arr = np.loadtxt(path, dtype=float)
    arr = arr.reshape(1, -1) if arr.ndim == 1 else arr

    if arr.shape[1] < 2:
        raise ValueError(f"Expected at least two columns in {path.name}, got shape {arr.shape}")

    return arr[:, :2]


def collect_inp_datasets(model_dir: Path) -> InpDatasetMap:
    """Load all `.inp` files in a model directory.

    Returns:
        Mapping filename -> ndarray with shape (n_rows, 2)
    """
    datasets: InpDatasetMap = {}

    if not model_dir.exists():
        raise FileNotFoundError(f"Model directory not found: {model_dir}")

    for file_path in sorted(model_dir.glob("*.inp")):
        datasets[file_path.name] = load_inp_file(file_path)

    return datasets


def summarize_inp_datasets(datasets: InpDatasetMap) -> list[dict[str, float | int | str]]:
    """Create a compact summary table for loaded `.inp` datasets."""
    summary: list[dict[str, float | int | str]] = []
    for filename, arr in sorted(datasets.items()):
        values = arr[:, 0]
        depth = arr[:, 1]
        summary.append(
            {
                "file": filename,
                "rows": int(arr.shape[0]),
                "value_min": float(np.min(values)),
                "value_max": float(np.max(values)),
                "depth_min": float(np.min(depth)),
                "depth_max": float(np.max(depth)),
            }
        )
    return summary


def plot_inp_profiles(
    datasets: InpDatasetMap,
    *,
    ncols: int = 3,
    figsize_per_subplot: Tuple[float, float] = (4.8, 3.2),
    xscale: str = "linear",
    invert_depth_axis: bool = True,
) -> tuple[plt.Figure, np.ndarray]:
    """Plot all `.inp` profiles in a compact subplot grid.

    Args:
        datasets: Mapping filename -> 2-column array [value, depth].
        ncols: Number of subplot columns.
        figsize_per_subplot: Width/height per subplot.
        xscale: Matplotlib x-axis scale (e.g., 'linear', 'log').
        invert_depth_axis: If True, depth increases downward.
    """
    if not datasets:
        raise ValueError("No datasets provided for plotting.")

    n_items = len(datasets)
    ncols = max(1, int(ncols))
    nrows = int(np.ceil(n_items / ncols))

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(figsize_per_subplot[0] * ncols, figsize_per_subplot[1] * nrows),
        squeeze=False,
    )
    axes_flat = axes.flatten()

    for ax, (filename, arr) in zip(axes_flat, sorted(datasets.items())):
        values = arr[:, 0]
        depth = arr[:, 1]
        species = Path(filename).stem

        ax.plot(values, depth, color="steelblue", linewidth=2.0)
        ax.set_title(species)
        ax.set_xlabel("Initial value")
        ax.set_ylabel("Depth")
        ax.set_xscale(xscale)
        ax.grid(True, alpha=0.3)
        if invert_depth_axis:
            ax.invert_yaxis()

    for ax in axes_flat[n_items:]:
        ax.set_visible(False)

    fig.suptitle("Initial conditions from .inp files", fontsize=14)
    fig.tight_layout(rect=[0, 0.02, 1, 0.98])
    return fig, axes


def plot_single_inp_profile(
    datasets: InpDatasetMap,
    filename: str,
    *,
    ax: Optional[plt.Axes] = None,
    xscale: str = "linear",
    invert_depth_axis: bool = True,
) -> plt.Axes:
    """Plot one selected profile from loaded `.inp` datasets."""
    if filename not in datasets:
        available = ", ".join(sorted(datasets.keys()))
        raise KeyError(f"{filename} not found. Available files: {available}")

    arr = datasets[filename]
    values = arr[:, 0]
    depth = arr[:, 1]
    species = Path(filename).stem

    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4))

    ax.plot(values, depth, color="steelblue", linewidth=2.2)
    ax.set_title(f"Initial profile: {species}")
    ax.set_xlabel("Initial value")
    ax.set_ylabel("Depth")
    ax.set_xscale(xscale)
    ax.grid(True, alpha=0.3)
    if invert_depth_axis:
        ax.invert_yaxis()

    return ax
