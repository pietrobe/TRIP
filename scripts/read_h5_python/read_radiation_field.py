"""
HDF5 Radiation Field Reader
============================
Reads and plots data from output_3D_field_mpi.h5 produced by the MPI radiation
transfer code.  The file layout is:

  /emergent_angular_grid/   – angular quadrature metadata
  /frequencies_grid/        – frequency axis
  /radiation_field/         – I, QI_pc, UI_pc, VI_pc  (float32)
                              shape: (Nx, Ny, Nz, N_inc, N_az, N_freq)
                                      [8,  8, 17,   8,   16,   96 ]
"""

from __future__ import annotations

import argparse
import dataclasses
from pathlib import Path
from typing import Tuple

import h5py
import numpy as np

import matplotlib

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from wavelength_utils import freq_to_angstrom, vacuum_to_air




# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclasses.dataclass
class AngularGrid:
    """Container for the emergent angular quadrature."""
    N_directions: int
    N_inclination_angles: int
    N_azimuthal_angles: int
    inclination_angles: np.ndarray    # shape (N_inclination_angles,)
    azimuthal_angles: np.ndarray      # shape (N_azimuthal_angles,)
    inclinations_indices: np.ndarray  # shape (N_inclination_angles,)
    azimuthal_indices: np.ndarray     # shape (N_azimuthal_angles,)

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"AngularGrid(N_directions={self.N_directions}, "
            f"N_inc={self.N_inclination_angles}, N_az={self.N_azimuthal_angles})"
        )


@dataclasses.dataclass
class RadiationFieldSlice:
    """
    Four Stokes-like 3D arrays at a single spatial point (i, j, k).

    Each array has shape  (N_inclination, N_azimuthal, N_frequencies).
    """
    i: int
    j: int
    k: int
    I:     np.ndarray   # Stokes I
    QI_pc: np.ndarray   # Q/I  (percent)
    UI_pc: np.ndarray   # U/I  (percent)
    VI_pc: np.ndarray   # V/I  (percent)

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"RadiationFieldSlice(ijk=({self.i},{self.j},{self.k}), "
            f"shape={self.I.shape})"
        )


@dataclasses.dataclass
class RadiationFieldMap2D:
    """
    Four Stokes-like 2D maps on a selected spatial plane for fixed direction and frequency.

    ``axis`` selects the plane: ``"xy"`` (fixed z), ``"xz"`` (fixed y), or ``"yz"`` (fixed x).

    Each array has shape (N_dim1, N_dim2).
    """
    k: int
    freq_idx: int
    inc_idx: int
    az_idx: int
    axis: str = "xy"
    I: np.ndarray = None
    QI_pc: np.ndarray = None
    UI_pc: np.ndarray = None
    VI_pc: np.ndarray = None
    height: float | None = None
    mu: float | None = None
    chi: float | None = None
    y_vals: np.ndarray | None = None
    delta: float | None = None

    def __repr__(self) -> str:  # pragma: no cover
        return (
            "RadiationFieldMap2D("
            f"axis={self.axis}, k={self.k}, freq_idx={self.freq_idx}, "
            f"dir=({self.inc_idx},{self.az_idx}), "
            f"shape={self.I.shape})"
        )


# ---------------------------------------------------------------------------
# Reader functions
# ---------------------------------------------------------------------------

def read_angular_grid(filepath: str | Path) -> AngularGrid:
    """
    Read the ``/emergent_angular_grid`` group and return an :class:`AngularGrid`.

    Parameters
    ----------
    filepath:
        Path to the HDF5 file.

    Returns
    -------
    AngularGrid
    """
    with h5py.File(filepath, "r") as f:
        grp = f["emergent_angular_grid"]
        return AngularGrid(
            N_directions=int(grp["N_directions"][0]),
            N_inclination_angles=int(grp["N_inclination_angles"][0]),
            N_azimuthal_angles=int(grp["N_azimuthal_angles"][0]),
            inclination_angles=grp["inclination_angles"][:],
            azimuthal_angles=grp["azimuthal_angles"][:],
            inclinations_indices=grp["inclinations_indices"][:],
            azimuthal_indices=grp["azimuthal_indices"][:],
        )


def read_frequencies(filepath: str | Path) -> np.ndarray:
    """
    Read ``/frequencies_grid/frequencies`` and return a 1-D NumPy array.

    Parameters
    ----------
    filepath:
        Path to the HDF5 file.

    Returns
    -------
    np.ndarray  shape ``(N_freq,)``
    """
    with h5py.File(filepath, "r") as f:
        return f["frequencies_grid/frequencies"][:]


def read_radiation_field(
    filepath: str | Path,
    i: int,
    j: int,
    k: int,
) -> RadiationFieldSlice:
    """
    Read the radiation field at spatial grid point ``(i, j, k)`` and return
    four 3-D NumPy arrays, one per Stokes component.

    The datasets in the file have shape
    ``(Nx, Ny, Nz, N_inclination, N_azimuthal, N_frequencies)``.
    Only the slice ``[i, j, k, :, :, :]`` is loaded into memory, giving arrays
    of shape ``(N_inclination, N_azimuthal, N_frequencies)``.

    Parameters
    ----------
    filepath:
        Path to the HDF5 file.
    i, j, k:
        Spatial indices (must be within the dataset extents).

    Returns
    -------
    RadiationFieldSlice
    """
    with h5py.File(filepath, "r") as f:
        rf = f["radiation_field"]

        # Validate indices against the actual dataset shape
        shape = rf["I"].shape   # (Nx, Ny, Nz, N_inc, N_az, N_freq)
        Nx, Ny, Nz = shape[:3]
        if not (0 <= i < Nx and 0 <= j < Ny and 0 <= k < Nz):
            raise IndexError(
                f"Indices (i={i}, j={j}, k={k}) out of bounds for "
                f"spatial grid ({Nx}, {Ny}, {Nz})."
            )

        # HDF5 hyperslab: only load the slice we need
        try:
            I_arr = rf["I"][i, j, k, :, :, :]
        except OSError:
            N_inc, N_az, N_freq = shape[3:]
            I_arr = np.full((N_inc, N_az, N_freq), np.nan, dtype=np.float32)
        QI_arr = rf["QI_pc"][i, j, k, :, :, :]
        UI_arr = rf["UI_pc"][i, j, k, :, :, :]
        VI_arr = rf["VI_pc"][i, j, k, :, :, :]

    return RadiationFieldSlice(
        i=i, j=j, k=k,
        I=I_arr, QI_pc=QI_arr, UI_pc=UI_arr, VI_pc=VI_arr,
    )


def map_filed(
    filepath: str | Path,
    fixed_idx: int,
    freq_idx: int,
    inc_idx: int,
    az_idx: int,
    axis: str = "xy",
) -> RadiationFieldMap2D:
    """
    Build 2D maps over a spatial plane at fixed third index, frequency, and direction.

    Parameters
    ----------
    filepath:
        Path to the HDF5 file.
    fixed_idx:
        Index of the fixed spatial dimension (z for xy, y for xz, x for yz).
    freq_idx:
        Frequency index.
    inc_idx, az_idx:
        Direction indices in inclination and azimuthal grids.
    axis:
        Plane to extract: ``"xy"`` (default), ``"xz"``, or ``"yz"``.

    Returns
    -------
    RadiationFieldMap2D
        2D maps for I, Q/I, U/I, V/I.
    """
    with h5py.File(filepath, "r") as f:
        rf = f["radiation_field"]

        shape = rf["I"].shape  # (Nx, Ny, Nz, N_inc, N_az, N_freq)
        Nx, Ny, Nz, N_inc, N_az, N_freq = shape

        if not (0 <= freq_idx < N_freq):
            raise IndexError(
                f"freq_idx={freq_idx} out of bounds [0, {N_freq - 1}].")
        if not (0 <= inc_idx < N_inc):
            raise IndexError(
                f"inc_idx={inc_idx} out of bounds [0, {N_inc - 1}].")
        if not (0 <= az_idx < N_az):
            raise IndexError(f"az_idx={az_idx} out of bounds [0, {N_az - 1}].")

        y_vals = None

        if axis == "xy":
            if not (0 <= fixed_idx < Nz):
                raise IndexError(
                    f"k={fixed_idx} out of bounds [0, {Nz - 1}] for axis=xy.")
            I_map = rf["I"][:, :, fixed_idx, inc_idx, az_idx, freq_idx]
            QI_map = rf["QI_pc"][:, :, fixed_idx, inc_idx, az_idx, freq_idx]
            UI_map = rf["UI_pc"][:, :, fixed_idx, inc_idx, az_idx, freq_idx]
            VI_map = rf["VI_pc"][:, :, fixed_idx, inc_idx, az_idx, freq_idx]
            expected = (Nx, Ny)
        elif axis == "xz":
            if not (0 <= fixed_idx < Ny):
                raise IndexError(
                    f"j={fixed_idx} out of bounds [0, {Ny - 1}] for axis=xz.")
            heights = f["geometry_3D/heights"][:]
            y_vals = heights[::-1].copy()
            I_map = rf["I"][:, fixed_idx, :, inc_idx, az_idx, freq_idx][:, ::-1]
            QI_map = rf["QI_pc"][:, fixed_idx, :, inc_idx, az_idx, freq_idx][:, ::-1]
            UI_map = rf["UI_pc"][:, fixed_idx, :, inc_idx, az_idx, freq_idx][:, ::-1]
            VI_map = rf["VI_pc"][:, fixed_idx, :, inc_idx, az_idx, freq_idx][:, ::-1]
            expected = (Nx, Nz)
        elif axis == "yz":
            if not (0 <= fixed_idx < Nx):
                raise IndexError(
                    f"i={fixed_idx} out of bounds [0, {Nx - 1}] for axis=yz.")
            heights = f["geometry_3D/heights"][:]
            y_vals = heights[::-1].copy()
            I_map = rf["I"][fixed_idx, :, :, inc_idx, az_idx, freq_idx][:, ::-1]
            QI_map = rf["QI_pc"][fixed_idx, :, :, inc_idx, az_idx, freq_idx][:, ::-1]
            UI_map = rf["UI_pc"][fixed_idx, :, :, inc_idx, az_idx, freq_idx][:, ::-1]
            VI_map = rf["VI_pc"][fixed_idx, :, :, inc_idx, az_idx, freq_idx][:, ::-1]
            expected = (Ny, Nz)
        else:
            raise ValueError(f"Unknown axis '{axis}'; expected 'xy', 'xz', or 'yz'.")

        try:
            delta = float(f["geometry_3D/delta"][0])
        except (KeyError, IndexError):
            delta = None

    assert I_map.shape == expected, (
        f"Expected shape {expected}, got {I_map.shape}")

    return RadiationFieldMap2D(
        k=fixed_idx,
        freq_idx=freq_idx,
        inc_idx=inc_idx,
        az_idx=az_idx,
        axis=axis,
        I=I_map,
        QI_pc=QI_map,
        UI_pc=UI_map,
        VI_pc=VI_map,
        y_vals=y_vals,
        delta=delta,
    )


map_field = map_filed


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_radiation_field(
    rf_slice: RadiationFieldSlice,
    frequencies: np.ndarray,
    angular_grid: AngularGrid | None = None,
    max_directions: int = 16,
    figsize: Tuple[float, float] = (14, 10),
    as_percent: bool = False,
) -> plt.Figure:
    """
    Plot all four Stokes components in a 2 × 2 grid.

    Each sub-plot shows the field as a function of frequency.  Each curve
    corresponds to one (inclination, azimuth) direction pair.  Only the first
    ``max_directions`` direction pairs are drawn to keep the figure readable.

    Parameters
    ----------
    rf_slice:
        Output of :func:`read_radiation_field`.
    frequencies:
        1-D array of frequencies (output of :func:`read_frequencies`).
    angular_grid:
        Optional.  If supplied, inclination angles are used as legend labels.
    max_directions:
        Maximum number of direction curves drawn per panel.
    figsize:
        Figure size in inches.

    Returns
    -------
    matplotlib Figure
    """
    percent_label = r"$Q/I$ [%]"
    components = [
        ("I",     rf_slice.I,     r"Stokes $I$",        "C0"),
        ("Q/I",   rf_slice.QI_pc, percent_label,          "C1"),
        ("U/I",   rf_slice.UI_pc, r"$U/I$ [%]", "C2"),
        ("V/I",   rf_slice.VI_pc, r"$V/I$ [%]", "C3"),
    ]

    N_inc, N_az, N_freq = rf_slice.I.shape
    assert N_freq == len(frequencies), (
        f"Frequency axis mismatch: field has {N_freq} points but "
        f"frequencies array has {len(frequencies)}."
    )

    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.38, wspace=0.32)

    # Build a flat list of (inc_idx, az_idx) pairs, truncated to max_directions
    dir_pairs = [
        (ii, ia)
        for ii in range(N_inc)
        for ia in range(N_az)
    ][:max_directions]

    # Colour map: one shade per inclination angle
    inc_cmap = plt.cm.viridis(np.linspace(0.1, 0.9, N_inc))

    for panel_idx, (label, data, ylabel, _) in enumerate(components):
        ax = fig.add_subplot(gs[panel_idx // 2, panel_idx % 2])

        for (ii, ia) in dir_pairs:
            curve = data[ii, ia, :]   # shape (N_freq,)
            lw = 1.2
            alpha = 0.75

            # Build label only for first azimuthal index (avoid legend clutter)
            if ia == 0 and angular_grid is not None:
                inc_deg = np.degrees(angular_grid.inclination_angles[ii])
                lbl = rf"$\theta={inc_deg:.1f}°$"
            elif ia == 0:
                lbl = rf"inc {ii}"
            else:
                lbl = None

            ax.plot(
                frequencies,
                curve,
                color=inc_cmap[ii],
                lw=lw,
                alpha=alpha,
                label=lbl,
            )

        ax.set_title(
            rf"{label}  |  point $({rf_slice.i},{rf_slice.j},{rf_slice.k})$",
            fontsize=11,
        )
        ax.set_xlabel("Frequency", fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.tick_params(labelsize=8)
        ax.grid(True, lw=0.4, alpha=0.4)

        if panel_idx == 0:
            handles, labels = ax.get_legend_handles_labels()
            if handles:
                ax.legend(
                    handles, labels,
                    fontsize=7,
                    loc="upper right",
                    ncol=2,
                    framealpha=0.6,
                )

    fig.suptitle(
        f"Radiation field at spatial point (i={rf_slice.i}, "
        f"j={rf_slice.j}, k={rf_slice.k})",
        fontsize=13,
        y=0.98,
    )
    return fig


def plot_radiation_field_direction(
    rf_slice: RadiationFieldSlice,
    frequencies: np.ndarray,
    angular_grid: AngularGrid,
    inclination_rad: float,
    azimuthal_rad: float,
    figsize: Tuple[float, float] = (14, 10),
    as_percent: bool = False,
) -> plt.Figure:
    """
    Plot all four Stokes components in a 2 × 2 grid for a **single direction**.

    The direction is specified by ``inclination_rad`` and ``azimuthal_rad``
    (both in radians).  The nearest available grid angles are selected
    automatically; the actual angles used are shown in the panel titles and
    all labels remain in radians.

    Parameters
    ----------
    rf_slice:
        Output of :func:`read_radiation_field`.
    frequencies:
        1-D array of frequencies (output of :func:`read_frequencies`).
    angular_grid:
        :class:`AngularGrid` used to resolve angle indices.
    inclination_rad:
        Desired inclination angle in radians.
    azimuthal_rad:
        Desired azimuthal angle in radians.
    figsize:
        Figure size in inches.

    Returns
    -------
    matplotlib Figure

    Raises
    ------
    ValueError
        If ``angular_grid`` contains no angles.
    """
    if angular_grid.N_inclination_angles == 0 or angular_grid.N_azimuthal_angles == 0:
        raise ValueError("angular_grid contains no angle entries.")

    # --- find nearest inclination index ---
    ii = int(np.argmin(np.abs(angular_grid.inclination_angles - inclination_rad)))
    ia = int(np.argmin(np.abs(angular_grid.azimuthal_angles - azimuthal_rad)))

    theta_actual = angular_grid.inclination_angles[ii]   # radians
    phi_actual = angular_grid.azimuthal_angles[ia]     # radians

    # --- extract 1-D frequency curves ---
    I_curve = rf_slice.I[ii, ia, :]
    QI_curve = rf_slice.QI_pc[ii, ia, :]
    UI_curve = rf_slice.UI_pc[ii, ia, :]
    VI_curve = rf_slice.VI_pc[ii, ia, :]

    percent_label = r"$Q/I$ [%]"
    components = [
        (r"Stokes $I$", I_curve,  r"$I$",           "C0"),
        (r"$Q/I$",      QI_curve, percent_label,  "C1"),
        (r"$U/I$",      UI_curve, r"$U/I$ [%]",    "C2"),
        (r"$V/I$",      VI_curve, r"$V/I$ [%]",    "C3"),
    ]

    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.32)

    wl_Angstroms = freq_to_angstrom(frequencies)
    wl_Angstroms = vacuum_to_air(wl_Angstroms)

    for panel_idx, (title, curve, ylabel, color) in enumerate(components):
        ax = fig.add_subplot(gs[panel_idx // 2, panel_idx % 2])

        ax.plot(
            wl_Angstroms,
            curve,
            color=color,
            lw=1.8,
            label=(
                rf"$\theta={theta_actual:.4f}$ rad,  "
                rf"$\chi={phi_actual:.4f}$ rad"
            ),
        )

        ax.set_title(title, fontsize=11)
        ax.set_xlabel("Wavelength [Å]", fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.tick_params(labelsize=8)
        ax.grid(True, lw=0.4, alpha=0.4)
        # ax.legend(fontsize=8, loc="best", framealpha=0.6)

    # Warn if the requested angles differed noticeably from the grid
    inc_diff = abs(theta_actual - inclination_rad)
    az_diff = abs(phi_actual - azimuthal_rad)

    fig.suptitle(
        rf"Radiation field  —  spatial point $({rf_slice.i},{rf_slice.j},{rf_slice.k})$"
        "\n"
        rf"Angular grid: $\theta={theta_actual:.4f}$ rad, $\mu={np.cos(theta_actual):.4f}$, $\chi={phi_actual:.4f}$ rad"
        + (
            rf"   ⚠ $\Delta\theta={inc_diff:.3f}$, $\Delta\chi={az_diff:.3f}$ rad"
            if (inc_diff > 0.05 or az_diff > 0.05) else ""
        ),
        fontsize=10,
        y=0.99,
    )
    return fig


_AXIS_LABELS = {
    "xy": ("x", "y"),
    "xz": ("x", "height"),
    "yz": ("y", "height"),
}

_FIXED_LABELS = {
    "xy": "k",
    "xz": "j",
    "yz": "i",
}


def plot_radiation_field_maps(
    rf_map: RadiationFieldMap2D,
    figsize: Tuple[float, float] = (12, 10),
    as_percent: bool = False,
    log_I: bool = False,
) -> plt.Figure:
    """Plot I, Q/I, U/I, V/I spatial maps in a 2 × 2 layout."""
    from matplotlib.colors import LogNorm

    percent_label = r"$Q/I$ [%]"
    components = [
        (r"Stokes $I$", rf_map.I),
        (percent_label, rf_map.QI_pc),
        (r"$U/I$ [%]", rf_map.UI_pc),
        (r"$V/I$ [%]", rf_map.VI_pc),
    ]

    fig, axes = plt.subplots(2, 2, figsize=figsize, constrained_layout=True)

    xlabel, ylabel = _AXIS_LABELS.get(rf_map.axis, ("x", "y"))

    dx = rf_map.delta if rf_map.delta is not None else 1.0

    for panel_idx, (ax, (title, arr)) in enumerate(zip(axes.ravel(), components)):
        norm = LogNorm() if (log_I and panel_idx == 0) else None
        ny, nx = arr.T.shape
        if rf_map.y_vals is not None:
            extent = (0, (nx - 1) * dx, rf_map.y_vals[0], rf_map.y_vals[-1])
        else:
            extent = (0, (nx - 1) * dx, 0, (ny - 1) * dx)
        im = ax.imshow(arr.T, origin="lower", aspect="auto",
                       cmap="viridis", extent=extent, norm=norm)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel(xlabel, fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.tick_params(labelsize=8)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    fixed_name = _FIXED_LABELS.get(rf_map.axis, "k")

    if rf_map.axis == "xy" and rf_map.height is not None:
        fixed_str = f"z={rf_map.height:.1f}"
    else:
        fixed_str = f"{fixed_name}={rf_map.k}"

    if rf_map.mu is not None and rf_map.chi is not None:
        title = (
            f"Spatial maps at {fixed_str}, "
            f"freq_idx={rf_map.freq_idx}, "
            f"dir=(\u03bc={rf_map.mu:.4f}, \u03c7={rf_map.chi:.4f} rad)"
        )
    else:
        title = (
            f"Spatial maps at {fixed_str}, "
            f"freq_idx={rf_map.freq_idx}, "
            f"dir=(inc={rf_map.inc_idx}, az={rf_map.az_idx})"
        )
    fig.suptitle(title, fontsize=12,)
    return fig


# ---------------------------------------------------------------------------
# NaN detection
# ---------------------------------------------------------------------------

_COORD_LABELS = {
    "xy": ("i", "j"),
    "xz": ("i", "k"),
    "yz": ("j", "k"),
}


def _report_nan(rf_map: RadiationFieldMap2D) -> None:
    """Print warnings for NaN/Inf/zero pixels in each Stokes component."""
    labels = _COORD_LABELS.get(rf_map.axis, ("d0", "d1"))
    components = [
        ("I",     rf_map.I),
        ("Q/I",   rf_map.QI_pc),
        ("U/I",   rf_map.UI_pc),
        ("V/I",   rf_map.VI_pc),
    ]
    for name, arr in components:
        n_nan = int(np.isnan(arr).sum())
        n_inf = int(np.isinf(arr).sum())
        report_zero = name == "I"
        n_zero = int((arr == 0).sum()) if report_zero else 0
        if n_nan == 0 and n_inf == 0 and n_zero == 0:
            continue
        parts = []
        if n_nan > 0:
            parts.append(f"{n_nan} NaN")
        if n_inf > 0:
            parts.append(f"{n_inf} Inf")
        if n_zero > 0:
            parts.append(f"{n_zero} zero")
        bad = ~np.isfinite(arr) | ((arr == 0) if report_zero else False)
        n_bad = int(bad.sum())
        idxs = np.argwhere(bad)
        n_show = min(n_bad, 10)
        examples = ", ".join(
            f"({lbls[0]}={idx[0]}, {lbls[1]}={idx[1]})"
            for idx, lbls in zip(idxs[:n_show], [labels] * n_show)
        )
        if n_bad > n_show:
            examples += f", ... ({n_bad - n_show} more)"
        print(f"  ⚠ Abnormal in {name}: {', '.join(parts)} at {examples}")


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

def _build_parser() -> "argparse.ArgumentParser":
    import argparse

    parser = argparse.ArgumentParser(
        prog="hdf5_radiation_field",
        description=(
            "Plot the radiation field at a spatial point or as 2D spatial maps.\n\n"
            "Example:\n"
            "  python hdf5_radiation_field.py \\\n"
            "      --file output/output_3D_field_mpi.h5 \\\n"
            "      --coords 3 3 8 \\\n"
            "      --dir 2 4\n\n"
            "Map mode:\n"
            "  python hdf5_radiation_field.py --map 8 --dir 2 4 --freq-idx 10"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-f",
        "--file",
        metavar="PATH",
        default="output/output_3D_field_mpi.h5",
        help="Path to the HDF5 file  (default: output/output_3D_field_mpi.h5)",
    )
    parser.add_argument(
        "-c",
        "--coords",
        nargs=3,
        type=int,
        metavar=("I", "J", "K"),
        default=[3, 3, 8],
        help="Spatial grid indices i j k  (default: 3 3 8)",
    )
    parser.add_argument(
        "-d",
        "--dir",
        nargs=2,
        type=int,
        metavar=("INC_IDX", "AZ_IDX"),
        default=None,
        help=(
            "Direction to plot as integer indices into the angular grid: "
            "inclination index and azimuthal index. "
            "If omitted, indices (0, 0) are used."
        ),
    )
    parser.add_argument(
        "-a",
        "--adir",
        nargs=2,
        type=float,
        metavar=("MU", "CHI"),
        default=None,
        help=(
            "Direction to plot using angular values (mu, chi), where "
            "mu = cos(theta) and chi is in radians. "
            "The nearest angular-grid indices are selected automatically."
        ),
    )
    parser.add_argument(
        "--map",
        metavar="IDX",
        type=int,
        default=None,
        help=(
            "If set, produce 2D maps on the spatial plane selected by --axis. "
            "IDX is the fixed index along the third axis "
            "(z-index for xy, y-index for xz, x-index for yz). "
            "Direction is taken from --dir/--adir and frequency index from --freq-idx."
        ),
    )

    parser.add_argument(
        "--axis",
        choices=["xy", "xz", "yz"],
        default="xy",
        help="Spatial plane for --map mode (default: xy).",
    )

    parser.add_argument(
        "--freq-idx",
        metavar="N",
        type=int,
        default=0,
        help="Frequency index used in --map mode (default: 0).",
    )

    parser.add_argument(
        "--output",
        metavar="FILE",
        default=None,
        help=(
            "Save the figure to this path instead of (only) displaying it. "
            "The format is inferred from the extension, e.g. plot.png or plot.pdf."
        ),
    )
    parser.add_argument(
        "--percent",
        action="store_true",
        default=False,
        help="Show Q/I, U/I, V/I labels with percent symbol.",
    )
    parser.add_argument(
        "--log",
        action="store_true",
        default=False,
        help="Plot Stokes I in logarithmic scale.",
    )
    return parser


def _print_h5_structure(filepath: str | Path) -> None:
    """Walk the HDF5 tree and print all groups/datasets in a compact tree format."""
    with h5py.File(filepath, "r") as f:
        def _show(name: str, indent: str = "") -> None:
            item = f[name]
            if isinstance(item, h5py.Group):
                label = name.rstrip("/").split("/")[-1]
                if label:
                    print(f"{indent}{label}/")
                for child in item:
                    _show(f"{name}/{child}" if name else child, indent + ("  " if label else ""))
            elif isinstance(item, h5py.Dataset):
                label = name.split("/")[-1]
                total = np.prod(item.shape) if item.shape else 1
                val_str = ""
                if total <= 20:
                    data = item[()]
                    if item.dtype.kind == 'f':
                        vals = np.atleast_1d(data)
                        if vals.size <= 6:
                            val_str = "  [" + ", ".join(f"{v:10.4e}" for v in vals) + "]"
                        else:
                            val_str = f"  [{vals[0]:10.4e}  ...  {vals[-1]:10.4e}]  N={vals.size}"
                    elif item.dtype.kind in ('i', 'u'):
                        vals = np.atleast_1d(data)
                        if vals.size <= 12:
                            val_str = "  " + str(vals.tolist())
                        else:
                            val_str = f"  [{vals[0]}  ...  {vals[-1]}]  N={vals.size}"
                print(f"{indent}{label:25s}  {str(item.shape):25s}  {str(item.dtype):10s}{val_str}")

        print(f"\n{'─' * 75}")
        _show("/")
        print(f"{'─' * 75}\n")


def main() -> None:
    import argparse

    parser = _build_parser()
    args = parser.parse_args()

    path = Path(args.file)
    if not path.exists():
        parser.error(f"File not found: {path.resolve()}")

    _print_h5_structure(path)

    i, j, k = args.coords

    print(f"Reading  : {path.resolve()}")
    print(f"Coords   : i={i}, j={j}, k={k}")

    ang = read_angular_grid(path)
    freqs = read_frequencies(path)

    print(f"  Angular grid  : {ang}")
    print(
        f"  Frequencies   : {freqs.shape}, min={freqs.min():.4g}, max={freqs.max():.4g}")

    if args.dir is not None and args.adir is not None:
        parser.error("Use either --dir or --adir, not both.")

    # Resolve direction indices — default to (0, 0)
    if args.adir is not None:
        mu_req = float(args.adir[0])
        chi_req = float(args.adir[1])

        if not (-1.0 <= mu_req <= 1.0):
            parser.error("--adir MU must be in [-1, 1].")

        mu_grid = np.cos(ang.inclination_angles)
        ii = int(np.argmin(np.abs(mu_grid - mu_req)))
        ia = int(np.argmin(np.abs(ang.azimuthal_angles - chi_req)))

        theta_req = float(np.arccos(mu_req))
        theta_plot = float(ang.inclination_angles[ii])
        chi_plot = float(ang.azimuthal_angles[ia])
    else:
        ii, ia = args.dir if args.dir is not None else (0, 0)
        theta_req = float(ang.inclination_angles[ii])
        theta_plot = float(ang.inclination_angles[ii])
        chi_plot = float(ang.azimuthal_angles[ia])

    N_inc = ang.N_inclination_angles
    N_az = ang.N_azimuthal_angles
    if not (0 <= ii < N_inc):
        parser.error(
            f"--dir inclination index {ii} out of range [0, {N_inc - 1}]")
    if not (0 <= ia < N_az):
        parser.error(
            f"--dir azimuthal index {ia} out of range [0, {N_az - 1}]")

    if args.adir is not None:
        print(
            "Direction: from --adir "
            f"(mu={mu_req:.4f}, chi={chi_req:.4f} rad) -> "
            f"inc_idx={ii}, az_idx={ia}"
        )
    elif args.dir is None:
        print(f"  --dir/--adir not set, using indices (0, 0)")

    print(
        f"Direction: inc_idx={ii}, az_idx={ia}  ->  "
        f"theta={theta_plot:.4f} rad, chi={chi_plot:.4f} rad, mu={np.cos(theta_plot):.4f}"
    )

    if args.map is not None:
        axis = args.axis
        fixed_idx = int(args.map)
        freq_idx = int(args.freq_idx)

        if not (0 <= freq_idx < len(freqs)):
            parser.error(
                f"--freq-idx {freq_idx} out of range [0, {len(freqs) - 1}]")

        rf_map = map_filed(path, fixed_idx, freq_idx, ii, ia, axis=axis)

        with h5py.File(path, "r") as f:
            heights = f["geometry_3D/heights"][:]
        if axis == "xy":
            rf_map.height = float(heights[fixed_idx])
        rf_map.mu = float(np.cos(theta_plot))
        rf_map.chi = float(chi_plot)

        if args.percent:
            rf_map = RadiationFieldMap2D(
                k=rf_map.k,
                freq_idx=rf_map.freq_idx,
                inc_idx=rf_map.inc_idx,
                az_idx=rf_map.az_idx,
                axis=rf_map.axis,
                I=rf_map.I,
                QI_pc=rf_map.QI_pc * 100.0,
                UI_pc=rf_map.UI_pc * 100.0,
                VI_pc=rf_map.VI_pc * 100.0,
                height=rf_map.height,
                mu=rf_map.mu,
                chi=rf_map.chi,
                y_vals=rf_map.y_vals,
                delta=rf_map.delta,
            )

        _report_nan(rf_map)

        print(f"  Radiation map : {rf_map}")
        print(f"  freq[{freq_idx}] = {float(freqs[freq_idx]):.6e}")

        fig = plot_radiation_field_maps(rf_map, as_percent=args.percent, log_I=args.log)
    else:
        rf = read_radiation_field(path, i, j, k)

        if args.percent:
            rf = RadiationFieldSlice(
                i=rf.i, j=rf.j, k=rf.k,
                I=rf.I,
                QI_pc=rf.QI_pc * 100.0,
                UI_pc=rf.UI_pc * 100.0,
                VI_pc=rf.VI_pc * 100.0,
            )

        print(f"  Radiation slice: {rf}")
        fig = plot_radiation_field_direction(rf, freqs, ang, theta_req, chi_req, as_percent=args.percent) if args.adir is not None else plot_radiation_field_direction(rf, freqs, ang, theta_plot, chi_plot, as_percent=args.percent)

    if args.output:
        out = Path(args.output)
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"Figure saved -> {out.resolve()}")

    plt.show()


if __name__ == "__main__":
    main()
