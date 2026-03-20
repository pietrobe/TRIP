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

import dataclasses
from pathlib import Path
from typing import Tuple

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


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
            N_directions          = int(grp["N_directions"][0]),
            N_inclination_angles  = int(grp["N_inclination_angles"][0]),
            N_azimuthal_angles    = int(grp["N_azimuthal_angles"][0]),
            inclination_angles    = grp["inclination_angles"][:],
            azimuthal_angles      = grp["azimuthal_angles"][:],
            inclinations_indices  = grp["inclinations_indices"][:],
            azimuthal_indices     = grp["azimuthal_indices"][:],
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
        I_arr     = rf["I"]    [i, j, k, :, :, :]
        QI_arr    = rf["QI_pc"][i, j, k, :, :, :]
        UI_arr    = rf["UI_pc"][i, j, k, :, :, :]
        VI_arr    = rf["VI_pc"][i, j, k, :, :, :]

    return RadiationFieldSlice(
        i=i, j=j, k=k,
        I=I_arr, QI_pc=QI_arr, UI_pc=UI_arr, VI_pc=VI_arr,
    )


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_radiation_field(
    rf_slice: RadiationFieldSlice,
    frequencies: np.ndarray,
    angular_grid: AngularGrid | None = None,
    max_directions: int = 16,
    figsize: Tuple[float, float] = (14, 10),
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
    components = [
        ("I",     rf_slice.I,     r"Stokes $I$",        "C0"),
        ("Q/I",   rf_slice.QI_pc, r"$Q/I$ [%]",         "C1"),
        ("U/I",   rf_slice.UI_pc, r"$U/I$ [%]",         "C2"),
        ("V/I",   rf_slice.VI_pc, r"$V/I$ [%]",         "C3"),
    ]

    N_inc, N_az, N_freq = rf_slice.I.shape
    assert N_freq == len(frequencies), (
        f"Frequency axis mismatch: field has {N_freq} points but "
        f"frequencies array has {len(frequencies)}."
    )

    fig = plt.figure(figsize=figsize)
    gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.38, wspace=0.32)

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
            lw    = 1.2
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
    ia = int(np.argmin(np.abs(angular_grid.azimuthal_angles   - azimuthal_rad)))

    theta_actual = angular_grid.inclination_angles[ii]   # radians
    phi_actual   = angular_grid.azimuthal_angles[ia]     # radians

    # --- extract 1-D frequency curves ---
    I_curve     = rf_slice.I    [ii, ia, :]
    QI_curve    = rf_slice.QI_pc[ii, ia, :]
    UI_curve    = rf_slice.UI_pc[ii, ia, :]
    VI_curve    = rf_slice.VI_pc[ii, ia, :]

    components = [
        (r"Stokes $I$", I_curve,  r"$I$",        "C0"),
        (r"$Q/I$",      QI_curve, r"$Q/I$ [%]",  "C1"),
        (r"$U/I$",      UI_curve, r"$U/I$ [%]",  "C2"),
        (r"$V/I$",      VI_curve, r"$V/I$ [%]",  "C3"),
    ]

    fig = plt.figure(figsize=figsize)
    gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.32)
    
    wl_Angstroms = frequencies * 1e10 / 2.99792458e8  # convert from Hz to Angstroms

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
    az_diff  = abs(phi_actual   - azimuthal_rad)

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


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

def _build_parser() -> "argparse.ArgumentParser":
    import argparse

    parser = argparse.ArgumentParser(
        prog="hdf5_radiation_field",
        description=(
            "Plot the radiation field at a spatial point for a single direction.\n\n"
            "Example:\n"
            "  python hdf5_radiation_field.py \\\n"
            "      --file output/output_3D_field_mpi.h5 \\\n"
            "      --coords 3 3 8 \\\n"
            "      --dir 2 4"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--file",
        metavar="PATH",
        default="output/output_3D_field_mpi.h5",
        help="Path to the HDF5 file  (default: output/output_3D_field_mpi.h5)",
    )
    parser.add_argument(
        "--coords",
        nargs=3,
        type=int,
        metavar=("I", "J", "K"),
        default=[3, 3, 8],
        help="Spatial grid indices i j k  (default: 3 3 8)",
    )
    parser.add_argument(
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
        "--output",
        metavar="FILE",
        default=None,
        help=(
            "Save the figure to this path instead of (only) displaying it. "
            "The format is inferred from the extension, e.g. plot.png or plot.pdf."
        ),
    )
    return parser


def main() -> None:
    import argparse

    parser = _build_parser()
    args   = parser.parse_args()

    path = Path(args.file)
    if not path.exists():
        parser.error(f"File not found: {path.resolve()}")

    i, j, k = args.coords

    print(f"Reading  : {path.resolve()}")
    print(f"Coords   : i={i}, j={j}, k={k}")

    ang   = read_angular_grid(path)
    freqs = read_frequencies(path)
    rf    = read_radiation_field(path, i, j, k)

    print(f"  Angular grid  : {ang}")
    print(f"  Frequencies   : {freqs.shape}, min={freqs.min():.4g}, max={freqs.max():.4g}")
    print(f"  Radiation slice: {rf}")

    # Resolve direction indices — default to (0, 0)
    ii, ia = args.dir if args.dir is not None else (0, 0)

    N_inc = ang.N_inclination_angles
    N_az  = ang.N_azimuthal_angles
    if not (0 <= ii < N_inc):
        parser.error(f"--dir inclination index {ii} out of range [0, {N_inc - 1}]")
    if not (0 <= ia < N_az):
        parser.error(f"--dir azimuthal index {ia} out of range [0, {N_az - 1}]")

    theta = float(ang.inclination_angles[ii])
    chi   = float(ang.azimuthal_angles[ia])

    if args.dir is None:
        print(f"  --dir not set, using indices (0, 0)")
    print(f"Direction: inc_idx={ii}, az_idx={ia}  ->  "
          f"theta={theta:.4f} rad, chi={chi:.4f} rad, mu={np.cos(theta):.4f}")

    fig = plot_radiation_field_direction(rf, freqs, ang, theta, chi)

    if args.output:
        out = Path(args.output)
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"Figure saved -> {out.resolve()}")

    plt.show()
if __name__ == "__main__":
    main()
