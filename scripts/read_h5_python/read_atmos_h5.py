"""
atmosphere_reader.py
--------------------
Reader for atmosphere_mpi.h5 HDF5 files.

Provides:
  - read_metadata(path)  → (Geometry, FrequencyGrid, Atom)
  - read_point(path, i, j, k) → AtmosPoint
  - read_field(path, field_name) → np.ndarray  (shape 128×128×128)
"""

from __future__ import annotations

import h5py
import numpy as np
from dataclasses import dataclass
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


# Custom colormaps
CMAP_BLACKRED = LinearSegmentedColormap.from_list("blackred", ["black", "red"])


# Physical constants (CGS/Gauss units)
K_B = 1.380649e-16              # Boltzmann constant [erg/K]
ATOMIC_MASS_UNIT = 1.660539e-24  # atomic mass unit [g]
PLANCK_H = 6.62607015e-27       # Planck constant [erg·s]
SPEED_OF_LIGHT = 2.99792458e10  # speed of light [cm/s]


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class Geometry:
    N_x: int
    N_y: int
    N_z: int
    delta: float          # horizontal grid spacing [cm or m, as stored]
    height: np.ndarray    # 1-D array of height values, shape (N_z,)


@dataclass
class FrequencyGrid:
    N_frequencies: int
    frequency: np.ndarray   # 1-D array of frequencies, shape (N_frequencies,)


@dataclass
class Atom:
    atomic_number: int
    atomic_mass: float
    E_lower: float
    E_upper: float
    g_lower: float
    g_upper: float
    Aul: float
    jl2: int              # 2*J for lower level
    ju2: int              # 2*J for upper level


@dataclass
class AtmosPoint:
    """All atmospheric quantities at a single grid point (i, j, k)."""
    index_i: int
    index_j: int
    index_k: int
    temperature: float
    vmicro: float
    damping: float
    pop_lower_level: float
    pop_upper_level: float
    Cul: float  # inelastic collisional rate
    Qel: float  # elastic collisional rate
    magnetic_field_x: float
    magnetic_field_y: float
    magnetic_field_z: float
    bulk_velocity_x: float
    bulk_velocity_y: float
    bulk_velocity_z: float
    c_scat_opacity_sigma_c: float
    c_therm_emissivity_epsilon_c: float
    c_tot_opacity_K_c: float


# ---------------------------------------------------------------------------
# Field name mapping
#
# Keys   = Python / AtmosPoint attribute names (what the code uses)
# Values = actual HDF5 compound field names    (what the file stores)
#
# If the HDF5 file was written with different names, only this dict needs
# to be updated -- the rest of the code adapts automatically.
# Run  inspect_fields(path)  to print the names actually present in a file.
# ---------------------------------------------------------------------------

FIELD_MAP: dict[str, str] = {
    "index_i"                     : "index_i",
    "index_j"                     : "index_j",
    "index_k"                     : "index_k",
    "temperature"                 : "temperature",
    "vmicro"                      : "vmicro",
    "damping"                     : "damping",
    "pop_lower_level"             : "pop_lower_level",
    "pop_upper_level"             : "pop_upper_level",
    "Cul"                         : "Cul",
    "Qel"                         : "Qel",
    "magnetic_field_x"            : "magnetic_field_x",
    "magnetic_field_y"            : "magnetic_field_y",
    "magnetic_field_z"            : "magnetic_field_z",
    "bulk_velocity_x"             : "bulk_velocity_x",
    "bulk_velocity_y"             : "bulk_velocity_y",
    "bulk_velocity_z"             : "bulk_velocity_z",
    "c_scat_opacity_sigma_c"      : "c_scat_opacity_sigma_c",
    "c_therm_emissivity_epsilon_c": "c_therm_emissivity_epsilon_c",
    "c_tot_opacity_K_c"           : "c_tot_opacity_K_c",
}

# HDF5-side names, used to validate read_field / read_field_slice.
ATMOS_FIELDS: tuple[str, ...] = tuple(FIELD_MAP.values())


def inspect_fields(path: str) -> tuple[str, ...]:
    """
    Return the compound field names exactly as stored in *path*.

    Use this to diagnose mismatches between ``FIELD_MAP`` and the file,
    then update the values in ``FIELD_MAP`` accordingly.

    Example
    -------
    >>> print(inspect_fields("atmosphere_mpi.h5"))
    ('index_i', 'index_j', ..., 'Cul', 'Qel', ...)
    """
    with h5py.File(path, "r") as f:
        dtype = f["atmos"].dtype
    names = dtype.names
    print("Fields found in file:")
    for n in names:
        print(f"  {n!r}")
    return names


def inspect_structure(path: str) -> None:
    """
    Print the top-level structure of the HDF5 file.
    """
    with h5py.File(path, "r") as f:
        def print_tree(name, obj):
            print(f"  {name}")
            if hasattr(obj, "shape"):
                print(f"    shape: {obj.shape}, dtype: {obj.dtype}")
        print("File structure:")
        f.visititems(print_tree)


# ---------------------------------------------------------------------------
# Helper: scalar extraction from length-1 HDF5 datasets
# ---------------------------------------------------------------------------

def _scalar(ds) -> float | int:
    """Return a plain Python scalar from an HDF5 dataset."""
    val = ds[()]
    if isinstance(val, np.ndarray):
        val = val.item()
    if np.issubdtype(type(val), np.integer):
        return int(val)
    return float(val)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def _get_group(f, *names):
    for name in names:
        if name in f:
            return f[name]
    raise KeyError(f"None of {names} found in file")


def _find_closest_index(arr: np.ndarray, value: float) -> int:
    """
    Find the index of the closest value in a sorted array.

    Parameters
    ----------
    arr : np.ndarray
        1-D sorted array of values.
    value : float
        Value to search for.

    Returns
    -------
    int
        Index of the closest value.
    """
    return int(np.argmin(np.abs(arr - value)))


def read_metadata(path: str) -> tuple[Geometry, FrequencyGrid, Atom]:
    """
    Open *path* and read the geometry, frequency grid, and atom groups.

    Parameters
    ----------
    path : str
        Path to the HDF5 file (e.g. "atmosphere_mpi.h5").

    Returns
    -------
    geometry : Geometry
    freq_grid : FrequencyGrid
    atom : Atom
    """
    with h5py.File(path, "r") as f:
        shape = f["atmos"].shape
        N_x, N_y, N_z = shape

        geo = _get_group(f, "geometry", "Geometry", "grid")
        if "height" in geo:
            height = geo["height"][:]
        elif "z" in geo:
            height = geo["z"][:]
        else:
            height = np.arange(N_z, dtype=float)

        delta = 1.0
        if "delta" in geo:
            delta = _scalar(geo["delta"])
        elif "cell_size" in geo:
            delta = _scalar(geo["cell_size"])

        geometry = Geometry(
            N_x    = N_x,
            N_y    = N_y,
            N_z    = N_z,
            delta  = delta,
            height = height,
        )

        try:
            fg = _get_group(f, "frequency_grid", "FrequencyGrid", "frequencies")
            if "frequency" in fg:
                frequency = fg["frequency"][:]
            elif "freq" in fg:
                frequency = fg["freq"][:]
            else:
                frequency = np.array([])

            N_freq = len(frequency)
            if "N_frequencies" in fg:
                N_freq = _scalar(fg["N_frequencies"])

            freq_grid = FrequencyGrid(
                N_frequencies = N_freq,
                frequency     = frequency,
            )
        except KeyError:
            freq_grid = FrequencyGrid(N_frequencies = 0, frequency = np.array([]))

        at = _get_group(f, "atom", "Atom", "atom_data")
        if at.shape == () or at.shape == (1,):
            if at.shape == (1,):
                at = at[0]
            atom = Atom(
                atomic_number = int(at["atomic_number"]),
                atomic_mass   = float(at["atomic_mass"]),
                E_lower       = float(at["E_lower"]),
                E_upper       = float(at["E_upper"]),
                g_lower       = float(at["g_lower"]),
                g_upper       = float(at["g_upper"]),
                jl2           = int(at["jl2"]),
                ju2           = int(at["ju2"]),
                Aul           = float(at["Aul"]),
            )
        else:
            atom = Atom(
                atomic_number = _scalar(at["atomic_number"]) if "atomic_number" in at else 0,
                atomic_mass   = _scalar(at["atomic_mass"]) if "atomic_mass" in at else 0.0,
                E_lower       = _scalar(at["E_lower"]) if "E_lower" in at else 0.0,
                E_upper       = _scalar(at["E_upper"]) if "E_upper" in at else 0.0,
                g_lower       = _scalar(at["g_lower"]) if "g_lower" in at else 1.0,
                g_upper       = _scalar(at["g_upper"]) if "g_upper" in at else 1.0,
                Aul           = _scalar(at["Aul"]) if "Aul" in at else 0.0,
                jl2           = _scalar(at["jl2"]) if "jl2" in at else 0,
                ju2           = _scalar(at["ju2"]) if "ju2" in at else 0,
            )

    return geometry, freq_grid, atom


def read_point(path: str, i: int, j: int, k: int) -> AtmosPoint:
    """
    Read all atmospheric quantities at grid point (i, j, k).

    The HDF5 dataset ``atmosphere/atmos_data`` has shape (128, 128, 128)
    and uses a compound dtype, so indexing ``ds[i, j, k]`` returns a
    numpy void with named fields.

    Parameters
    ----------
    path : str
        Path to the HDF5 file.
    i, j, k : int
        Zero-based grid indices along x, y, z axes respectively.

    Returns
    -------
    AtmosPoint
    """
    with h5py.File(path, "r") as f:
        ds   = f["atmos"]
        Nx, Ny, Nz = ds.shape
        if not (0 <= i < Nx and 0 <= j < Ny and 0 <= k < Nz):
            raise IndexError(
                f"Index ({i}, {j}, {k}) out of bounds for dataset "
                f"with shape ({Nx}, {Ny}, {Nz})."
            )
        rec = ds[i, j, k]   # numpy void (compound record)

    # Use FIELD_MAP so HDF5 name changes only need updating in one place.
    fm = FIELD_MAP
    return AtmosPoint(
        index_i                     = int(rec[fm["index_i"]]),
        index_j                     = int(rec[fm["index_j"]]),
        index_k                     = int(rec[fm["index_k"]]),
        temperature                 = float(rec[fm["temperature"]]),
        vmicro                      = float(rec[fm["vmicro"]]),
        damping                     = float(rec[fm["damping"]]),
        pop_lower_level             = float(rec[fm["pop_lower_level"]]),
        pop_upper_level             = float(rec[fm["pop_upper_level"]]),
        Cul                         = float(rec[fm["Cul"]]),
        Qel                         = float(rec[fm["Qel"]]),
        magnetic_field_x            = float(rec[fm["magnetic_field_x"]]),
        magnetic_field_y            = float(rec[fm["magnetic_field_y"]]),
        magnetic_field_z            = float(rec[fm["magnetic_field_z"]]),
        bulk_velocity_x             = float(rec[fm["bulk_velocity_x"]]),
        bulk_velocity_y             = float(rec[fm["bulk_velocity_y"]]),
        bulk_velocity_z             = float(rec[fm["bulk_velocity_z"]]),
        c_scat_opacity_sigma_c      = float(rec[fm["c_scat_opacity_sigma_c"]]),
        c_therm_emissivity_epsilon_c= float(rec[fm["c_therm_emissivity_epsilon_c"]]),
        c_tot_opacity_K_c           = float(rec[fm["c_tot_opacity_K_c"]]),
    )


def read_field_slice(
    path: str,
    field_name: str,
    *,
    i: int | None = None,
    j: int | None = None,
    k: int | None = None,
) -> np.ndarray:
    """
    Extract a 2-D slice of a named compound field by fixing one coordinate.

    Exactly one of *i*, *j*, or *k* must be provided.  The returned array
    spans the two free axes in their natural (x, y, z) order:

    +-----------+------------------+------------------+
    | Fixed axis| Returned shape   | Axis meaning     |
    +===========+==================+==================+
    | i (x)     | (N_y, N_z)       | rows=y, cols=z   |
    +-----------+------------------+------------------+
    | j (y)     | (N_x, N_z)       | rows=x, cols=z   |
    +-----------+------------------+------------------+
    | k (z)     | (N_x, N_y)       | rows=x, cols=y   |
    +-----------+------------------+------------------+

    Parameters
    ----------
    path : str
        Path to the HDF5 file.
    field_name : str
        Name of the field inside the compound dtype (see ``ATMOS_FIELDS``).
    i : int, optional
        Fix the x-axis index.
    j : int, optional
        Fix the y-axis index.
    k : int, optional
        Fix the z-axis index.

    Returns
    -------
    np.ndarray, shape (N_a, N_b)
        2-D slice of the requested field.

    Raises
    ------
    ValueError
        If *field_name* is invalid, or if not exactly one axis is fixed.
    IndexError
        If the fixed index is out of bounds.
    """
    if field_name not in ATMOS_FIELDS:
        raise ValueError(
            f"'{field_name}' is not a valid field name.\n"
            f"Valid fields are: {ATMOS_FIELDS}"
        )

    fixed = {name: val for name, val in (("i", i), ("j", j), ("k", k))
             if val is not None}
    if len(fixed) != 1:
        raise ValueError(
            "Exactly one of 'i', 'j', or 'k' must be provided, "
            f"but got: {fixed or 'none'}."
        )

    axis_name, idx = next(iter(fixed.items()))

    with h5py.File(path, "r") as f:
        ds = f["atmos"]
        Nx, Ny, Nz = ds.shape

        bounds = {"i": Nx, "j": Ny, "k": Nz}
        limit  = bounds[axis_name]
        if not (0 <= idx < limit):
            raise IndexError(
                f"Index {axis_name}={idx} is out of bounds "
                f"for axis of size {limit}."
            )

        # Use h5py field extraction + HDF5-level slicing so only the
        # requested slab is read from disk.
        field_ds = ds.fields(field_name)
        if axis_name == "i":
            arr = field_ds[idx, :, :]   # shape (Ny, Nz)
        elif axis_name == "j":
            arr = field_ds[:, idx, :]   # shape (Nx, Nz)
        else:                           # k
            arr = field_ds[:, :, idx]   # shape (Nx, Ny)

    return arr


def read_field(path: str, field_name: str) -> np.ndarray:
    """
    Extract a named field from the compound dataset as a 3-D NumPy array.

    Parameters
    ----------
    path : str
        Path to the HDF5 file.
    field_name : str
        Name of the field inside the compound dtype.
        Must be one of: ``ATMOS_FIELDS``.

    Returns
    -------
    np.ndarray, shape (N_x, N_y, N_z)
        Array containing the requested field for every grid point.

    Raises
    ------
    ValueError
        If *field_name* is not a recognised compound field.
    """
    if field_name not in ATMOS_FIELDS:
        raise ValueError(
            f"'{field_name}' is not a valid field name.\n"
            f"Valid fields are: {ATMOS_FIELDS}"
        )

    with h5py.File(path, "r") as f:
        ds = f["atmos"]
        # h5py supports direct field extraction via ds.fields(name)[...]
        # which avoids loading the entire compound dataset into memory.
        arr = ds.fields(field_name)[...]

    return arr   # shape (128, 128, 128), dtype matches the HDF5 field type


def rescale_atmos(input_path: str, output_path: str, step_k: int, force_even_k: bool = False, rescale_h_zero: bool = False) -> None:
    """
    Rescale an atmosphere HDF5 file by picking one in every N planes along the height (k) axis.

    This reads the input file, selects every *step_k*-th height plane (indices 0, step_k, 2*step_k, ...),
    constructs a new HDF5 file with the same structure but reduced in the k-dimension,
    updates the index_k field in the compound dataset, and writes the new height vector.

    Parameters
    ----------
    input_path : str
        Path to the input HDF5 file.
    output_path : str
        Path to the output HDF5 file.
    step_k : int
        Step size for selecting height planes (e.g., step_k=2 keeps every 2nd plane).
    force_even_k : bool
        If True, force the resulting number of heights to be even (default: False).
    rescale_h_zero : bool
        If True, rescale heights so the minimum height is zero (default: False).
        When False, the original height values are preserved.
    """
    if step_k < 1:
        raise ValueError("step_k must be >= 1")

    with h5py.File(input_path, "r") as f_in:
        geometry, freq_grid, atom = read_metadata(input_path)
        N_x = geometry.N_x
        N_y = geometry.N_y
        N_z_orig = geometry.N_z
        height_orig = geometry.height

        N_z_new = (N_z_orig + step_k - 1) // step_k
        
        if force_even_k:
            if N_z_new % 2 != 0:
                N_z_new -= 1
            k_indices = np.arange(0, N_z_orig, step_k)[1:N_z_new+1]
            new_height = height_orig[k_indices]
        else:
            k_indices = np.arange(0, N_z_orig, step_k)[:N_z_new]
            new_height = height_orig[k_indices]
        
        if rescale_h_zero:
            new_height = new_height - new_height.min()
            
        print(f"Original N_z: {N_z_orig}, New N_z: {N_z_new}, step_k: {step_k}, force_even_k: {force_even_k}, rescale_h_zero: {rescale_h_zero}")

        with h5py.File(output_path, "w") as f_out:
            f_out.create_group("geometry")
            f_out["geometry/N_height"] = np.array([N_z_new], dtype=np.int64)
            f_out["geometry/N_ij"] = np.array([N_x], dtype=np.int64)
            f_out["geometry/delta"] = np.array([geometry.delta], dtype=np.float64)
            f_out["geometry/height"] = new_height.astype(np.float32)

            f_out.create_group("frequency_grid")
            f_out["frequency_grid/N_frequencies"] = np.array([freq_grid.N_frequencies], dtype=np.int32)
            if freq_grid.N_frequencies > 0:
                f_out["frequency_grid/frequency"] = freq_grid.frequency

            atom_dtype = np.dtype([
                ("atomic_number", "<i4"),
                ("atomic_mass", "<f8"),
                ("E_lower", "<f8"),
                ("E_upper", "<f8"),
                ("g_lower", "<f8"),
                ("g_upper", "<f8"),
                ("jl2", "<i4"),
                ("ju2", "<i4"),
                ("Aul", "<f8"),
            ])
            f_out.create_dataset("atom_data", shape=(1,), dtype=atom_dtype)
            at = f_out["atom_data"]
            at[0] = (
                atom.atomic_number,
                atom.atomic_mass,
                atom.E_lower,
                atom.E_upper,
                atom.g_lower,
                atom.g_upper,
                atom.jl2,
                atom.ju2,
                atom.Aul,
            )

            f_out.create_group("continuum")

            atmos_dtype = f_in["atmos"].dtype
            atmos_ds_out = f_out.create_dataset(
                "atmos",
                shape=(N_x, N_y, N_z_new),
                dtype=atmos_dtype,
                compression="gzip",
                compression_opts=4,
            )

            for i in range(N_x):
                for j in range(N_y):
                    slab_orig = f_in["atmos"][i, j, k_indices]
                    slab_new = np.empty_like(slab_orig)
                    for fname in atmos_dtype.names:
                        if fname in ("index_i", "index_j"):
                            slab_new[fname] = slab_orig[fname]
                        elif fname == "index_k":
                            slab_new[fname] = np.arange(N_z_new, dtype=slab_orig[fname].dtype)
                        else:
                            slab_new[fname] = slab_orig[fname]
                    atmos_ds_out[i, j, :] = slab_new

            if "continuum/continuum_data" in f_in:
                cont_dtype = f_in["continuum/continuum_data"].dtype
                cont_ds_out = f_out.create_dataset(
                    "continuum/continuum_data",
                    shape=(N_x, N_y, N_z_new),
                    dtype=cont_dtype,
                    compression="gzip",
                    compression_opts=4,
                )
                for i in range(N_x):
                    for j in range(N_y):
                        cont_ds_out[i, j, :] = f_in["continuum/continuum_data"][i, j, k_indices]

            for name in ("c_scat_opacity_sigma_c", "c_therm_emissivity_epsilon_c", "c_tot_opacity_K_c"):
                if f"continuum/{name}" in f_in:
                    ds_in = f_in[f"continuum/{name}"]
                    shape_in = ds_in.shape
                    new_shape = (N_x, N_y, N_z_new, shape_in[3])
                    ds_out = f_out.create_dataset(
                        f"continuum/{name}",
                        shape=new_shape,
                        dtype=ds_in.dtype,
                        compression="gzip",
                        compression_opts=4,
                    )
                    for i in range(N_x):
                        for j in range(N_y):
                            ds_out[i, j, :, :] = ds_in[i, j, k_indices, :]


def analyze_height_grid(path: str) -> None:
    """
    Print height grid statistics: min/max height and min/max spacing between adjacent heights.

    Parameters
    ----------
    path : str
        Path to the HDF5 file.
    """
    RESET = "\033[0m"
    BOLD = "\033[1m"
    BLUE = "\033[94m"
    GREEN = "\033[92m"
    CYAN = "\033[96m"
    YELLOW = "\033[93m"

    geometry, _, _ = read_metadata(path)
    h = geometry.height

    h_min = h.min()
    h_max = h.max()
    dh = np.diff(h)
    dh_min = dh.min()
    dh_max = dh.max()

    print(f"\n{'─' * 50}")
    print(f"{BOLD}Height grid analysis:{RESET}")
    print(f"  {BLUE}Number of heights:{RESET} {len(h)}")
    print(f"  {BLUE}Height range:{RESET}")
    print(f"    min: {GREEN}{h_min:.4e} cm{RESET} ({CYAN}{h_min/1e5:.4f} km{RESET})")
    print(f"    max: {GREEN}{h_max:.4e} cm{RESET} ({CYAN}{h_max/1e5:.4f} km{RESET})")
    print(f"  {BLUE}Spacing (adjacent heights):{RESET}")
    print(f"    min Δh: {YELLOW}{dh_min:.4e} cm{RESET}")
    print(f"    max Δh: {GREEN}{dh_max:.4e} cm{RESET}")


def check_damping_stats(path: str, plot: bool = False, axis: str = "xy", idx: int = 0) -> np.ndarray:
    """
    Compute global statistics of the damping factor comparison.

    Parameters
    ----------
    path : str
        Path to the HDF5 file.
    plot : bool
        If True, plot a slice of the relative error.
    axis : str
        Axis for the slice plot ("xy", "xz", or "yz").
    idx : int
        Index along the fixed axis for the slice plot.

    Returns
    -------
    np.ndarray
        Relative error array of shape (N_x, N_y, N_z).
    """
    RESET = "\033[0m"
    BOLD = "\033[1m"
    BLUE = "\033[94m"
    GREEN = "\033[92m"
    CYAN = "\033[96m"
    YELLOW = "\033[93m"
    RED = "\033[91m"

    geometry, _, atom = read_metadata(path)
    N_x, N_y, N_z = geometry.N_x, geometry.N_y, geometry.N_z

    Aul = atom.Aul
    atomic_mass = atom.atomic_mass
    E_upper = atom.E_upper
    E_lower = atom.E_lower

    Cul = read_field(path, "Cul")
    Qel = read_field(path, "Qel")
    vmicro = read_field(path, "vmicro")
    temperature = read_field(path, "temperature")
    damping = read_field(path, "damping")

    dE_cm = E_upper - E_lower
    nu0 = dE_cm * SPEED_OF_LIGHT
    m_real = atomic_mass * ATOMIC_MASS_UNIT
    dnd = (nu0 / SPEED_OF_LIGHT) * np.sqrt(2 * K_B * temperature / m_real + vmicro**2)
    a_theoretical = (Aul + Cul + Qel) / (4 * np.pi * dnd)

    valid_mask = ~np.isnan(a_theoretical) & ~np.isnan(damping) & (a_theoretical != 0)

    a_theo_valid = a_theoretical[valid_mask]
    damping_valid = damping[valid_mask]

    rel_err = np.zeros_like(damping)
    rel_err[valid_mask] = (a_theo_valid - damping_valid) / a_theo_valid * 100

    if valid_mask.sum() == 0:
        print(f"\n{RED}Error: No valid points for comparison.{RESET}")
        return rel_err

    valid_err = rel_err[valid_mask]

    min_err = valid_err.min()
    max_err = valid_err.max()
    mean_err = valid_err.mean()
    std_err = valid_err.std()

    min_idx = np.unravel_index(np.nanargmin(rel_err), rel_err.shape)
    max_idx = np.unravel_index(np.nanargmax(rel_err), rel_err.shape)

    a_theo_at_min = a_theoretical[min_idx]
    a_theo_at_max = a_theoretical[max_idx]
    damping_at_min = damping[min_idx]
    damping_at_max = damping[max_idx]

    a_theo_mean = a_theo_valid.mean()
    damping_mean = damping_valid.mean()

    print(f"\n{'─' * 50}")
    print(f"{BOLD}Damping factor comparison:{RESET}")
    print(f"  {BLUE}Formula used:{RESET}")
    print(f"    a_theoretical = (Aul + Cul + Qel) / (4π × dnd)")
    print(f"    ν₀ = (E_upper - E_lower) × c")
    print(f"    dnd = (ν₀ / c) × √(2×k_B×T/m + vmicro²)")
    print(f"    m = atomic_mass × atomic_mass_unit")
    print(f"  {BLUE}Relative error formula:{RESET}")
    print(f"    err = (a_theoretical - damping) / a_theoretical × 100 [%]")
    print(f"  {BLUE}Points analyzed:{RESET} {N_x} × {N_y} × {N_z} = {N_x * N_y * N_z}")
    print(f"  {BLUE}Valid points:{RESET} {valid_mask.sum()}")
    print(f"  {BLUE}Relative error:{RESET}")
    print(f"    min: {GREEN}{min_err:.4f}{RESET} at ({min_idx[0]}, {min_idx[1]}, {min_idx[2]})")
    print(f"      theoretical: {a_theo_at_min:.4e}, file: {damping_at_min:.4e}")
    print(f"    max: {RED}{max_err:.4f}{RESET} at ({max_idx[0]}, {max_idx[1]}, {max_idx[2]})")
    print(f"      theoretical: {a_theo_at_max:.4e}, file: {damping_at_max:.4e}")
    print(f"    mean: {CYAN}{mean_err:.4f}{RESET}")
    print(f"    std: {YELLOW}{std_err:.4f}{RESET}")
    print(f"  {BLUE}Average damping:{RESET}")
    print(f"    theoretical: {GREEN}{a_theo_mean:.4e}{RESET}")
    print(f"    file: {GREEN}{damping_mean:.4e}{RESET}")

    hist_counts, bin_edges = np.histogram(valid_err, bins=30)
    bin_width = bin_edges[1] - bin_edges[0]
    max_count = hist_counts.max()

    print(f"\n{BOLD}Histogram (30 classes):{RESET}")
    print(f"{'Range':<26} {'Count':<12} {'Bar':<30}")
    print(f"{'─' * 68}")
    for i in range(len(hist_counts)):
        lo = bin_edges[i]
        hi = bin_edges[i+1]
        bar_len = int(hist_counts[i] / max_count * 25)
        bar = "█" * bar_len
        range_str = f"[{lo:>8.2f}, {hi:>8.2f})"
        print(f"{range_str:<26} {hist_counts[i]:>12,} {bar}")
    print(f"{'─' * 68}")
    print(f"Total points: {len(valid_err):,}")

    if plot:
        plot_damping_error_slice(path, rel_err, axis, idx, geometry)

    return rel_err


def plot_damping_error_slice(path: str, rel_err: np.ndarray, axis: str, idx: int, geometry: Geometry) -> None:
    """
    Plot a slice of the damping relative error.

    Parameters
    ----------
    path : str
        Path to the HDF5 file.
    rel_err : np.ndarray
        Relative error array.
    axis : str
        Axis for the slice ("xy", "xz", or "yz").
    idx : int
        Index along the fixed axis.
    geometry : Geometry
        Geometry object for axis labels.
    """
    if axis == "xy":
        slice_data = rel_err[:, :, idx]
        xlabel, ylabel = "x [km]", "y [km]"
        x_vals = np.arange(geometry.N_x) * geometry.delta / 1e5
        y_vals = np.arange(geometry.N_y) * geometry.delta / 1e5
    elif axis == "xz":
        slice_data = rel_err[:, idx, :]
        xlabel, ylabel = "x [km]", "z [km]"
        x_vals = np.arange(geometry.N_x) * geometry.delta / 1e5
        y_vals = geometry.height / 1e5
    else:
        slice_data = rel_err[idx, :, :]
        xlabel, ylabel = "y [km]", "z [km]"
        x_vals = np.arange(geometry.N_y) * geometry.delta / 1e5
        y_vals = geometry.height / 1e5

    fig, ax = plt.subplots()
    im = ax.imshow(slice_data.T, origin="lower", aspect="auto", cmap="seismic",
                   vmin=np.nanpercentile(rel_err, 2), vmax=np.nanpercentile(rel_err, 98))
    fig.colorbar(im, ax=ax, label="Relative error (%)")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"Damping relative error slice ({axis} at {idx})")
    plt.tight_layout()
    plt.show()


def plot_damping_error_vs_k(path: str, rel_err: np.ndarray, geometry: Geometry) -> None:
    """
    Plot the relative error of the damping factor vs height index k.

    Shows mean ± std band over x,y at each k, with min/max lines.

    Parameters
    ----------
    path : str
        Path to the HDF5 file.
    rel_err : np.ndarray
        Relative error array of shape (N_x, N_y, N_z).
    geometry : Geometry
        Geometry object for height axis.
    """
    mean = rel_err.mean(axis=(0, 1))
    std = rel_err.std(axis=(0, 1))
    min_err = rel_err.min(axis=(0, 1))
    max_err = rel_err.max(axis=(0, 1))

    k_indices = np.arange(geometry.N_z)
    height_km = geometry.height / 1e5

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.fill_between(height_km, mean - std, mean + std, alpha=0.3, color="steelblue", label="mean ± σ")
    ax.plot(height_km, mean, "k-", linewidth=1.5, label="mean")
    ax.plot(height_km, max_err, "r--", linewidth=0.8, label="max")
    ax.plot(height_km, min_err, "b--", linewidth=0.8, label="min")
    ax.axhline(0, color="gray", linestyle=":", linewidth=1)
    ax.set_xlabel("z [km]")
    ax.set_ylabel("Relative error (%)")
    ax.legend(loc="best")
    ax.set_title("Damping relative error vs height")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


def rescale_heights(path: str, index: int | None = None, dry_run: bool = False) -> bool:
    """
    Shift heights so they start from zero (or from height at given index).

    Parameters
    ----------
    path : str
        Path to the HDF5 file.
    index : int, optional
        If provided, subtract heights by the value at this index.
        If None, subtract by the minimum height.
    dry_run : bool
        If True, only print what would be done without actually updating the file.

    Returns
    -------
    bool
        True if the file was updated, False otherwise.
    """
    RESET = "\033[0m"
    BOLD = "\033[1m"
    BLUE = "\033[94m"
    GREEN = "\033[92m"
    CYAN = "\033[96m"
    YELLOW = "\033[93m"

    geometry, _, _ = read_metadata(path)
    h = geometry.height

    h_min = h.min()
    h_max = h.max()

    if index is not None:
        offset = float(h[index])
        source = f"height[{index}]"
    else:
        offset = h_min
        source = "minimum"

    h_new = h - offset
    new_min = h_new.min()
    new_max = h_new.max()

    print(f"\n{'─' * 50}")
    print(f"{BOLD}Height rescaling:{RESET}")
    print(f"  {BLUE}Current heights:{RESET}")
    print(f"    min: {GREEN}{h_min:.4e} cm{RESET} ({CYAN}{h_min/1e5:.4f} km{RESET})")
    print(f"    max: {GREEN}{h_max:.4e} cm{RESET} ({CYAN}{h_max/1e5:.4f} km{RESET})")
    print(f"  {BLUE}Offset (subtract):{RESET} {YELLOW}{offset:.4e} cm{RESET} (from {source})")
    print(f"  {BLUE}New heights:{RESET}")
    print(f"    min: {GREEN}{new_min:.4e} cm{RESET} ({CYAN}{new_min/1e5:.4f} km{RESET})")
    print(f"    max: {GREEN}{new_max:.4e} cm{RESET} ({CYAN}{new_max/1e5:.4f} km{RESET})")

    if dry_run:
        print(f"\n{BOLD}Dry run - no changes made.{RESET}")
        return False

    confirm = input(f"\n{BOLD}Apply changes?{RESET} (type '{GREEN}yes{RESET}' to confirm, default: no): ")
    if confirm.lower() != "yes":
        print("Cancelled.")
        return False

    with h5py.File(path, "r+") as f:
        f["geometry/height"][:] = h_new.astype(np.float32)
    print(f"{GREEN}Updated.{RESET}")
    return True


def extract_subdomain(path: str, i: int, j: int, N: int, output_path: str) -> None:
    """
    Extract a horizontal subdomain and write to a new HDF5 file.

    Parameters
    ----------
    path : str
        Path to the input HDF5 file.
    i : int
        Starting x-index.
    j : int
        Starting y-index.
    N : int
        Size of the subdomain (N × N in x-y).
    output_path : str
        Path to the output HDF5 file.
    """
    import os

    RESET = "\033[0m"
    BOLD = "\033[1m"
    BLUE = "\033[94m"
    GREEN = "\033[92m"
    CYAN = "\033[96m"
    YELLOW = "\033[93m"
    RED = "\033[91m"

    if os.path.exists(output_path):
        print(f"{RED}Error: Output file already exists: {output_path}{RESET}")
        return

    geometry, freq_grid, atom = read_metadata(path)
    N_x, N_y, N_z = geometry.N_x, geometry.N_y, geometry.N_z

    if i < 0 or j < 0 or i + N > N_x or j + N > N_y:
        print(f"{RED}Error: Subdomain ({i}, {j}) to ({i+N-1}, {j+N-1}) out of bounds ({N_x}, {N_y}){RESET}")
        return

    print(f"\n{'─' * 50}")
    print(f"{BOLD}Extract subdomain:{RESET}")
    print(f"  {BLUE}Input:{RESET} {path}")
    print(f"  {BLUE}Output:{RESET} {output_path}")
    print(f"  {BLUE}Region:{RESET} i=[{i}, {i+N}), j=[{j}, {j+N})")
    print(f"  {BLUE}Heights:{RESET} {N_z} (unchanged)")
    print(f"  {BLUE}Output size:{RESET} {GREEN}{N} × {N} × {N_z}{RESET}")

    confirm = input(f"\n{BOLD}Extract subdomain?{RESET} (type '{GREEN}yes{RESET}' to confirm, default: no): ")
    if confirm.lower() != "yes":
        print("Cancelled.")
        return

    with h5py.File(path, "r") as f_in:
        with h5py.File(output_path, "w") as f_out:
            f_out.create_group("geometry")
            f_out["geometry/N_height"] = np.array([N_z], dtype=np.int64)
            f_out["geometry/N_ij"] = np.array([N], dtype=np.int64)
            f_out["geometry/delta"] = np.array([geometry.delta], dtype=np.float64)
            f_out["geometry/height"] = geometry.height.astype(np.float32)

            f_out.create_group("frequency_grid")
            f_out["frequency_grid/N_frequencies"] = np.array([freq_grid.N_frequencies], dtype=np.int32)
            if freq_grid.N_frequencies > 0:
                f_out["frequency_grid/frequency"] = freq_grid.frequency

            atom_dtype = np.dtype([
                ("atomic_number", "<i4"),
                ("atomic_mass", "<f8"),
                ("E_lower", "<f8"),
                ("E_upper", "<f8"),
                ("g_lower", "<f8"),
                ("g_upper", "<f8"),
                ("jl2", "<i4"),
                ("ju2", "<i4"),
                ("Aul", "<f8"),
            ])
            f_out.create_dataset("atom_data", shape=(1,), dtype=atom_dtype)
            at = f_out["atom_data"]
            at[0] = (
                atom.atomic_number,
                atom.atomic_mass,
                atom.E_lower,
                atom.E_upper,
                atom.g_lower,
                atom.g_upper,
                atom.jl2,
                atom.ju2,
                atom.Aul,
            )

            f_out.create_group("continuum")

            atmos_dtype = f_in["atmos"].dtype
            atmos_ds_out = f_out.create_dataset(
                "atmos",
                shape=(N, N, N_z),
                dtype=atmos_dtype,
                compression="gzip",
                compression_opts=4,
            )

            slab_orig = f_in["atmos"][i:i+N, j:j+N, :][:]
            out_data = np.zeros((N, N, N_z), dtype=atmos_dtype)
            for fname in atmos_dtype.names:
                if fname not in ("index_i", "index_j", "index_k"):
                    out_data[fname][:] = slab_orig[fname]
            for ii in range(N):
                for jj in range(N):
                    for kk in range(N_z):
                        out_data[ii, jj, kk]["index_i"] = ii
                        out_data[ii, jj, kk]["index_j"] = jj
                        out_data[ii, jj, kk]["index_k"] = kk
            atmos_ds_out[:] = out_data

            if "continuum/continuum_data" in f_in:
                cont_dtype = f_in["continuum/continuum_data"].dtype
                cont_ds_out = f_out.create_dataset(
                    "continuum/continuum_data",
                    shape=(N, N, N_z),
                    dtype=cont_dtype,
                    compression="gzip",
                    compression_opts=4,
                )
                cont_ds_out[:] = f_in["continuum/continuum_data"][i:i+N, j:j+N, :]

            for name in ("c_scat_opacity_sigma_c", "c_therm_emissivity_epsilon_c", "c_tot_opacity_K_c"):
                if f"continuum/{name}" in f_in:
                    ds_in = f_in[f"continuum/{name}"]
                    shape_in = ds_in.shape
                    new_shape = (N, N, N_z, shape_in[3])
                    ds_out = f_out.create_dataset(
                        f"continuum/{name}",
                        shape=new_shape,
                        dtype=ds_in.dtype,
                        compression="gzip",
                        compression_opts=4,
                    )
                    ds_out[:] = ds_in[i:i+N, j:j+N, :, :]

    print(f"{GREEN}Extracted subdomain to: {output_path}{RESET}")


def scale_delta(path: str, factor: float, dry_run: bool = False) -> bool:
    """
    Multiply the geometry delta by a factor and save to the same file.

    Parameters
    ----------
    path : str
        Path to the HDF5 file.
    factor : float
        Factor to multiply the delta by.
    dry_run : bool
        If True, only print what would be done without actually updating the file.

    Returns
    -------
    bool
        True if the file was updated, False otherwise.
    """
    RESET = "\033[0m"
    BOLD = "\033[1m"
    GREEN = "\033[92m"
    BLUE = "\033[94m"
    YELLOW = "\033[93m"
    CYAN = "\033[96m"

    with h5py.File(path, "r") as f:
        delta = _scalar(f["geometry/delta"])
        new_delta = delta * factor
        geometry, _, _ = read_metadata(path)
        N_ij = geometry.N_x * geometry.N_y
        domain_size_cm = delta * (geometry.N_x - 1)
        domain_size_km = domain_size_cm / 1e5

    print(f"\n{'─' * 50}")
    print(f"{BOLD}Geometry update:{RESET}")
    print(f"  {BLUE}delta (cell size):{RESET}")
    print(f"    current: {delta:.4f} cm")
    print(f"    new:     {GREEN}{new_delta:.4f} cm{RESET} (×{factor})")
    print(f"  {BLUE}domain size:{RESET}")
    print(f"    N_x × N_y = {geometry.N_x} × {geometry.N_y} = {N_ij} points")
    print(f"    horizontal: {CYAN}{domain_size_cm:.4e} cm{RESET} ({YELLOW}{domain_size_km:.4f} km{RESET})")

    if dry_run:
        print(f"\n{BOLD}Dry run - no changes made.{RESET}")
        return False

    confirm = input(f"\n{BOLD}Apply changes?{RESET} (type '{GREEN}yes{RESET}' to confirm, default: no): ")
    if confirm.lower() != "yes":
        print("Cancelled.")
        return False

    with h5py.File(path, "r+") as f:
        f["geometry/delta"][()] = np.array([new_delta], dtype=np.float64)
    print(f"{GREEN}Updated.{RESET}")
    return True


def plot_z_profile(path: str, field: str = "temperature", geometry: Geometry | None = None, show_minmax: bool = False) -> None:
    """
    Plot the vertical profile of a field showing mean ± std, min and max over x,y at each z.

    Parameters
    ----------
    path : str
        Path to the HDF5 file.
    field : str
        Field name to plot.
    geometry : Geometry, optional
        Geometry object for height axis scaling.
    show_minmax : bool
        If True, show min/max lines (default: False).
    """
    print(f"Plotting z-profile of '{field}' …")
    data = read_field(path, field)
    mean = data.mean(axis=(0, 1))
    std = data.std(axis=(0, 1))
    min_val = data.min(axis=(0, 1))
    max_val = data.max(axis=(0, 1))
    z = geometry.height if geometry is not None else np.arange(len(mean))
    z_plot = z / 1e5

    fig, ax = plt.subplots()
    ax.fill_between(z_plot, mean - std, mean + std, alpha=0.3, label="mean ± σ")
    ax.plot(z_plot, mean, "k-", linewidth=1.5, label="mean")
    if show_minmax:
        ax.plot(z_plot, max_val, "r--", linewidth=0.8, label="max")
        ax.plot(z_plot, min_val, "b--", linewidth=0.8, label="min")
    ax.set_xlabel("z [km]")
    ax.set_ylabel(field)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 3))
    ax.legend()
    ax.set_title(f"Vertical profile: {field}")
    plt.tight_layout()
    plt.show()


def compute_z_profile(path: str, field: str, geometry: Geometry | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute average and std of a field at each height z.

    Parameters
    ----------
    path : str
        Path to the HDF5 file.
    field : str
        Field name to compute.
    geometry : Geometry, optional
        Geometry object for height axis.

    Returns
    -------
    tuple[np.ndarray, np.ndarray, np.ndarray]
        (z, mean, std) - height values in km, mean and std at each height.
    """
    data = read_field(path, field)
    mean = data.mean(axis=(0, 1))
    std = data.std(axis=(0, 1))
    z = geometry.height if geometry is not None else np.arange(len(mean))
    z_km = z / 1e5
    return z_km, mean, std


def plot_z_profile_with_std_band(path: str, field: str, geometry: Geometry | None = None, show_minmax: bool = False) -> None:
    """
    Plot the vertical profile of a field with a +/- std band.

    Parameters
    ----------
    path : str
        Path to the HDF5 file.
    field : str
        Field name to plot.
    geometry : Geometry, optional
        Geometry object for height axis scaling. If None, will be loaded from file.
    show_minmax : bool
        If True, show min/max lines (default: False).
    """
    if geometry is None:
        geometry, _, _ = read_metadata(path)

    print(f"Plotting z-profile of '{field}' with ±std band …")
    z_km, mean, std = compute_z_profile(path, field, geometry)

    if show_minmax:
        data = read_field(path, field)
        min_val = data.min(axis=(0, 1))
        max_val = data.max(axis=(0, 1))

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.fill_between(z_km, mean - std, mean + std, alpha=0.3, color="steelblue", label="mean ± σ")
    ax.plot(z_km, mean, "k-", linewidth=1.5, label="mean")
    if show_minmax:
        ax.plot(z_km, max_val, "r--", linewidth=0.8, label="max")
        ax.plot(z_km, min_val, "b--", linewidth=0.8, label="min")
    ax.set_xlabel("z [km]")
    ax.set_ylabel(field)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 3))
    ax.legend(loc="best")
    ax.set_title(f"Vertical profile: {field}")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


def plot_all_z_profiles(path: str, fields: list[str] | None = None, geometry: Geometry | None = None, show_minmax: bool = False, log: bool = False) -> None:
    """Plot multiple field profiles in a tiled layout."""
    if geometry is None:
        geometry, _, _ = read_metadata(path)

    if fields is None:
        fields = [
            "temperature",
            "vmicro",
            "damping",
            "pop_lower_level",
            "pop_upper_level",
            "Cul",
            "Qel",
            "magnetic_field_z",
            "bulk_velocity_z",
            "c_scat_opacity_sigma_c",
            "c_therm_emissivity_epsilon_c",
            "c_tot_opacity_K_c",
        ]

    LOG_SCALE_FIELDS = {
        "c_scat_opacity_sigma_c",
        "c_therm_emissivity_epsilon_c",
        "c_tot_opacity_K_c",
        "damping",
        "pop_lower_level",
        "pop_upper_level",
        "Cul",
        "Qel",
    }

    n_fields = len(fields)
    n_cols = min(3, n_fields)
    n_rows = (n_fields + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
    if n_fields == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    for idx, field in enumerate(fields):
        ax = axes[idx]
        z_km, mean, std = compute_z_profile(path, field, geometry)

        use_log_y = log and field in LOG_SCALE_FIELDS

        if use_log_y:
            ax.set_yscale("log")
            ax.set_ylim(bottom=min(mean), top=max(mean))
            lower = np.maximum(mean - std, 1e-100)
            upper = mean + std
            ax.fill_between(z_km, lower, upper, alpha=0.3, color="steelblue")
            ax.plot(z_km, mean, "k-", linewidth=1.5)
            ylabel = field
            if show_minmax:
                data = read_field(path, field)
                min_val = data.min(axis=(0, 1))
                max_val = data.max(axis=(0, 1))
                ax.plot(z_km, max_val, "r--", linewidth=0.8, label="max")
                ax.plot(z_km, min_val, "b--", linewidth=0.8, label="min")
                ax.legend(loc="best", fontsize=8)
        else:
            ax.fill_between(z_km, mean - std, mean + std, alpha=0.3, color="steelblue")
            ax.plot(z_km, mean, "k-", linewidth=1.5)
            ylabel = field
            if show_minmax:
                data = read_field(path, field)
                min_val = data.min(axis=(0, 1))
                max_val = data.max(axis=(0, 1))
                ax.plot(z_km, max_val, "r--", linewidth=0.8, label="max")
                ax.plot(z_km, min_val, "b--", linewidth=0.8, label="min")
                ax.legend(loc="best", fontsize=8)

        ax.set_xlabel("z [km]")
        ax.set_ylabel(ylabel)
        if not use_log_y:
            ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 3))
        ax.set_title(field)
        ax.grid(True, alpha=0.3)

    for idx in range(n_fields, len(axes)):
        axes[idx].set_visible(False)

    plt.tight_layout()
    plt.show()


def average_z_velocity(path: str) -> float:
    """
    Calculate the average of the z-component of the bulk velocity.

    Parameters
    ----------
    path : str
        Path to the HDF5 file.

    Returns
    -------
    float
        Average of bulk_velocity_z over all grid points.
    """
    vz = read_field(path, "bulk_velocity_z")
    return float(vz.mean())


def _format_ticks_km(ax, axis: str, delta: float, n_points: int) -> None:
    step = max(1, n_points // 5)
    ticks = np.arange(0, n_points + 1, step)
    if axis == "x":
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{t * delta / 1e5:.1f}" for t in ticks])
    else:
        ax.set_yticks(ticks)
        ax.set_yticklabels([f"{t * delta / 1e5:.1f}" for t in ticks])


def plot_slice_xy(path: str, field: str = "temperature", k: int = 0, log: bool = False, sign: bool = False, quiver: bool = False, geometry: Geometry | None = None, threshold: float | None = None, cmap: str = "blackred") -> None:
    """
    Plot a horizontal (x-y) slice of *field* from the HDF5 file at depth *k*.

    Parameters
    ----------
    path : str
        Path to the HDF5 file.
    field : str
        Compound field name to plot (default: ``'temperature'``).
        Use ``magnetic_field_abs`` to plot the magnitude of the magnetic field.
    k : int
        z-index of the x-y slice to plot (default: 0).
    log : bool
        If True, plot in log10 scale (default: False).
    sign : bool
        If True, use diverging colormap (seismic) to show sign (default: False).
    quiver : bool
        If True, overlay arrows showing magnetic field projection (default: False).
    geometry : Geometry, optional
        Geometry object for axis scaling (delta in cm → km).
    threshold : float, optional
        Maximum value for colormap scaling. Values above threshold are clamped.
    cmap : str
        Colormap name. Use "blackred" for black→red (default), "hot" for grayscale,
        "seismic" for diverging blue-white-red.
    """
    label = f"log10({field})" if log else field
    print(f"Plotting field '{field}' at k={k} (log={log}, sign={sign}, quiver={quiver}) …")
    slice_xy = _get_field_slice(path, field, axis="xy", idx=k)
    if threshold is not None:
        slice_xy = np.clip(slice_xy, None, threshold)
    if log:
        slice_xy = np.log10(slice_xy)
    if cmap == "blackred":
        cmap = CMAP_BLACKRED
    elif cmap == "seismic":
        cmap = "seismic"
    else:
        cmap = "hot"
    fig, ax = plt.subplots()
    im = ax.imshow(slice_xy.T, origin="lower", aspect="auto", cmap=cmap)
    fig.colorbar(im, ax=ax, label=label)
    ax.set_title(f"{label}  (k={k})")

    if geometry is not None:
        ax.set_xlabel("x [km]")
        ax.set_ylabel("y [km]")
        _format_ticks_km(ax, "x", geometry.delta, geometry.N_x)
        _format_ticks_km(ax, "y", geometry.delta, geometry.N_y)
    else:
        ax.set_xlabel("x index")
        ax.set_ylabel("y index")

    if quiver:
        Bx = read_field_slice(path, "magnetic_field_x", k=k)
        By = read_field_slice(path, "magnetic_field_y", k=k)
        _plot_quiver(ax, Bx, By, geometry=geometry, z_heights=None)

    plt.tight_layout()
    plt.show()


def plot_slice_xz(path: str, field: str = "temperature", j: int = 0, log: bool = False, sign: bool = False, quiver: bool = False, geometry: Geometry | None = None, threshold: float | None = None, cmap: str = "blackred") -> None:
    """
    Plot a vertical (x-z) slice of *field* from the HDF5 file at index *j*.

    Parameters
    ----------
    path : str
        Path to the HDF5 file.
    field : str
        Compound field name to plot (default: ``'temperature'``).
        Use ``magnetic_field_abs`` to plot the magnitude of the magnetic field.
    j : int
        y-index of the x-z slice to plot (default: 0).
    log : bool
        If True, plot in log10 scale (default: False).
    sign : bool
        If True, use diverging colormap (seismic) to show sign (default: False).
    quiver : bool
        If True, overlay arrows showing magnetic field projection (default: False).
    geometry : Geometry, optional
        Geometry object for axis scaling (x in km, z from height in cm → km).
    threshold : float, optional
        Maximum value for colormap scaling. Values above threshold are clamped.
    cmap : str
        Colormap name. Use "blackred" for black→red (default), "hot" for grayscale,
        "seismic" for diverging blue-white-red.
    """
    label = f"log10({field})" if log else field
    print(f"Plotting field '{field}' at j={j} (log={log}, sign={sign}, quiver={quiver}) …")
    slice_xz = _get_field_slice(path, field, axis="xz", idx=j)
    if threshold is not None:
        slice_xz = np.clip(slice_xz, None, threshold)
    if log:
        slice_xz = np.log10(slice_xz)
    if cmap == "blackred":
        cmap = CMAP_BLACKRED
    elif cmap == "seismic":
        cmap = "seismic"
    else:
        cmap = "hot"
    fig, ax = plt.subplots()
    im = ax.imshow(slice_xz.T, origin="lower", aspect="auto", cmap=cmap)
    fig.colorbar(im, ax=ax, label=label)
    ax.set_title(f"{label}  (j={j})")

    if geometry is not None:
        ax.set_xlabel("x [km]")
        ax.set_ylabel("z [km]")
        _format_ticks_km(ax, "x", geometry.delta, geometry.N_x)
        ax.set_yticks(np.linspace(0, len(geometry.height) - 1, 5))
        ax.set_yticklabels([f"{h / 1e5:.1f}" for h in np.linspace(geometry.height[0], geometry.height[-1], 5)])
    else:
        ax.set_xlabel("x index")
        ax.set_ylabel("z index")

    if quiver:
        Bx = read_field_slice(path, "magnetic_field_x", j=j)
        Bz = read_field_slice(path, "magnetic_field_z", j=j)
        _plot_quiver(ax, Bx, Bz, geometry=geometry, z_heights=geometry.height if geometry else None)

    plt.tight_layout()
    plt.show()


def plot_slice_yz(path: str, field: str = "temperature", i: int = 0, log: bool = False, sign: bool = False, quiver: bool = False, geometry: Geometry | None = None, threshold: float | None = None, cmap: str = "blackred") -> None:
    """
    Plot a vertical (y-z) slice of *field* from the HDF5 file at index *i*.

    Parameters
    ----------
    path : str
        Path to the HDF5 file.
    field : str
        Compound field name to plot (default: ``'temperature'``).
        Use ``magnetic_field_abs`` to plot the magnitude of the magnetic field.
    i : int
        x-index of the y-z slice to plot (default: 0).
    log : bool
        If True, plot in log10 scale (default: False).
    sign : bool
        If True, use diverging colormap (seismic) to show sign (default: False).
    quiver : bool
        If True, overlay arrows showing magnetic field projection (default: False).
    geometry : Geometry, optional
        Geometry object for axis scaling (y in km, z from height in cm → km).
    threshold : float, optional
        Maximum value for colormap scaling. Values above threshold are clamped.
    cmap : str
        Colormap name. Use "blackred" for black→red (default), "hot" for grayscale,
        "seismic" for diverging blue-white-red.
    """
    label = f"log10({field})" if log else field
    print(f"Plotting field '{field}' at i={i} (log={log}, sign={sign}, quiver={quiver}) …")
    slice_yz = _get_field_slice(path, field, axis="yz", idx=i)
    if threshold is not None:
        slice_yz = np.clip(slice_yz, None, threshold)
    if log:
        slice_yz = np.log10(slice_yz)
    if cmap == "blackred":
        cmap = CMAP_BLACKRED
    elif cmap == "seismic":
        cmap = "seismic"
    else:
        cmap = "hot"
    fig, ax = plt.subplots()
    im = ax.imshow(slice_yz.T, origin="lower", aspect="auto", cmap=cmap)
    fig.colorbar(im, ax=ax, label=label)
    ax.set_title(f"{label}  (i={i})")

    if geometry is not None:
        ax.set_xlabel("y [km]")
        ax.set_ylabel("z [km]")
        _format_ticks_km(ax, "x", geometry.delta, geometry.N_y)
        ax.set_yticks(np.linspace(0, len(geometry.height) - 1, 5))
        ax.set_yticklabels([f"{h / 1e5:.1f}" for h in np.linspace(geometry.height[0], geometry.height[-1], 5)])
    else:
        ax.set_xlabel("y index")
        ax.set_ylabel("z index")

    if quiver:
        By = read_field_slice(path, "magnetic_field_y", i=i)
        Bz = read_field_slice(path, "magnetic_field_z", i=i)
        _plot_quiver(ax, By, Bz, geometry=geometry, z_heights=geometry.height if geometry else None)

    plt.tight_layout()
    plt.show()


def _plot_quiver(ax, U: np.ndarray, V: np.ndarray, step: int = 8, geometry: Geometry | None = None, z_heights: np.ndarray | None = None) -> None:
    """
    Overlay arrows on an existing axes showing the vector field.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes to plot arrows on.
    U : np.ndarray
        X-component of the vector field.
    V : np.ndarray
        Y-component of the vector field.
    step : int
        Subsampling step for arrow placement (default: 8).
    geometry : Geometry, optional
        Geometry object for axis scaling (for xy slices).
    z_heights : np.ndarray, optional
        Height values in cm (for xz/yz slices).
    """
    Ny, Nx = U.shape
    y_grid, x_grid = np.mgrid[0:Ny:step, 0:Nx:step]
    u_sub = U[::step, ::step]
    v_sub = V[::step, ::step]
    norm = np.sqrt(u_sub**2 + v_sub**2)
    if norm.max() > 0:
        u_sub = u_sub / norm
        v_sub = v_sub / norm
    if geometry is not None:
        if z_heights is not None:
            x_grid = x_grid * geometry.delta / 1e5
            y_grid = z_heights[::step] / 1e5
        else:
            x_grid = x_grid * geometry.delta / 1e5
            y_grid = y_grid * geometry.delta / 1e5
    ax.quiver(x_grid, y_grid, u_sub, v_sub, norm, cmap="plasma",
              clim=[0, norm.max()], scale=20, pivot="mid")


def _get_field_slice(path: str, field: str, axis: str, idx: int) -> np.ndarray:
    """
    Internal helper to get a field slice, with special handling for
    ``magnetic_field_abs`` which computes the magnitude of the magnetic field.

    Parameters
    ----------
    path : str
        Path to the HDF5 file.
    field : str
        Field name or ``magnetic_field_abs`` for ||B||.
    axis : str
        One of ``xy``, ``xz``, ``yz``.
    idx : int
        Index along the fixed axis.

    Returns
    -------
    np.ndarray
        2-D slice of the requested field.
    """
    kwargs = {"i": idx, "j": idx, "k": idx}
    if field == "magnetic_field_abs":
        if axis == "xy":
            Bx = read_field_slice(path, "magnetic_field_x", k=idx)
            By = read_field_slice(path, "magnetic_field_y", k=idx)
            Bz = read_field_slice(path, "magnetic_field_z", k=idx)
        elif axis == "xz":
            Bx = read_field_slice(path, "magnetic_field_x", j=idx)
            By = read_field_slice(path, "magnetic_field_y", j=idx)
            Bz = read_field_slice(path, "magnetic_field_z", j=idx)
        else:
            Bx = read_field_slice(path, "magnetic_field_x", i=idx)
            By = read_field_slice(path, "magnetic_field_y", i=idx)
            Bz = read_field_slice(path, "magnetic_field_z", i=idx)
        return np.sqrt(Bx**2 + By**2 + Bz**2)

    if axis == "xy":
        return read_field_slice(path, field, k=idx)
    elif axis == "xz":
        return read_field_slice(path, field, j=idx)
    else:
        return read_field_slice(path, field, i=idx)


# ---------------------------------------------------------------------------
# Quick demo / smoke-test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse
    import sys

    parser = argparse.ArgumentParser(description="Read and inspect atmosphere HDF5 files.")
    parser.add_argument("-f", "--file", default="atmosphere_mpi.h5", help="Path to the HDF5 file")
    parser.add_argument("-i", type=int, default=0, help="x-index for point read")
    parser.add_argument("-j", type=int, default=0, help="y-index for point read")
    parser.add_argument("-k", type=int, default=0, help="z-index for point read")
    parser.add_argument("--plot", action="store_true", help="Plot example slices (requires matplotlib)")
    parser.add_argument("--axis", type=str, choices=["xy", "xz", "yz"], default="xy",
                        help="Plane to plot: xy (horizontal), xz (vertical), yz (vertical) (default: xy)")
    parser.add_argument("--idx", type=int, default=0,
                        help="Index along the sliced axis (default: 0)")
    parser.add_argument("--iheight", type=float, default=None,
                        help="Find closest height index for given height value (cm)")
    parser.add_argument("--field", type=str, default="temperature",
                        help="Field name to plot (default: temperature). Use 'magnetic_field_abs' for ||B||")
    parser.add_argument("--log", action="store_true",
                        help="Plot field in log scale")
    parser.add_argument("--sign", action="store_true",
                        help="Use diverging colormap (seismic) to show sign of positive/negative values")
    parser.add_argument("--quiver", action="store_true",
                        help="Plot arrows showing magnetic field projection on the slice")
    parser.add_argument("--avg-vz", action="store_true",
                        help="Print the average of the z-component of the bulk velocity")
    parser.add_argument("--no-scale", action="store_true",
                        help="Disable physical units scaling (use grid indices)")
    parser.add_argument("--z-profile", action="store_true",
                        help="Plot vertical profile (mean ± std) of the field vs height")
    parser.add_argument("--z-profile-minmax", action="store_true",
                        help="Show min/max lines in z-profile plot")
    parser.add_argument("--threshold", type=float, default=None,
                        help="Maximum value for colormap scaling")
    parser.add_argument("--cmap", type=str, default="blackred",
                        help="Colormap: blackred (default), hot, seismic")
    parser.add_argument("--rescale-k", type=int, default=0,
                        help="Rescale height axis by picking 1 in N planes (N = value)")
    parser.add_argument("--force-even-k", action="store_true",
                        help="Force the resulting number of heights to be even")
    parser.add_argument("--rescale-h-zero", action="store_true",
                        help="Rescale heights so the minimum height is zero (default: False)")
    parser.add_argument("--scale-delta", type=float, default=0.0,
                        help="Multiply delta by this factor in place (requires --out or in-place file)")
    parser.add_argument("--analyze-height", action="store_true",
                        help="Print height grid statistics (min/max height and spacing)")
    parser.add_argument("--rescale-heights", type=int, default=None, nargs="?", const=0,
                        help="Shift heights to start from zero. If index given, use height[index] as offset")
    parser.add_argument("--extract", nargs=3, type=int, metavar=("I", "J", "N"),
                        help="Extract subdomain: I J N - subdomain starting at (I,J) with size N×N")
    parser.add_argument("--out", type=str, default="",
                        help="Output file path for rescaled file")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would be done without making changes")
    parser.add_argument("--check-damping", action="store_true",
                        help="Compute global damping factor statistics")
    parser.add_argument("--plot-damping-err", action="store_true",
                        help="Plot relative error slice of damping factor")
    parser.add_argument("--plot-damping-vs-k", action="store_true",
                        help="Plot relative error of damping factor vs height index k")
    parser.add_argument("--plot-all-profiles", action="store_true",
                        help="Plot all z-profiles in tiled layout")
    parser.add_argument("--fields", nargs="+", default=None,
                        help="List of fields to plot (use with --plot-all-profiles)")
    args = parser.parse_args(args=None if sys.argv[1:] else ["--help"])

    hdf5_path = args.file

    print(f"Reading metadata from: {hdf5_path}\n")

    print("=== Compound fields in file ===")
    inspect_fields(hdf5_path)
    print()

    print("=== File structure ===")
    inspect_structure(hdf5_path)
    print()

    geometry, freq_grid, atom = read_metadata(hdf5_path)

    if args.rescale_k > 0:
        if not args.out:
            raise ValueError("--out is required when --rescale-k is specified")
        print(f"Rescaling height axis with step {args.rescale_k} (force_even_k={args.force_even_k}, rescale_h_zero={args.rescale_h_zero}) ...")
        rescale_atmos(hdf5_path, args.out, args.rescale_k, force_even_k=args.force_even_k, rescale_h_zero=args.rescale_h_zero)
        print(f"Output written to: {args.out}")
        exit(0)

    if args.scale_delta != 0.0:
        print(f"Scaling delta by {args.scale_delta} ...")
        scale_delta(hdf5_path, args.scale_delta, dry_run=args.dry_run)
        exit(0)

    if args.analyze_height:
        analyze_height_grid(hdf5_path)
        exit(0)

    if args.rescale_heights is not None:
        index = args.rescale_heights if args.rescale_heights > 0 else None
        rescale_heights(hdf5_path, index=index, dry_run=args.dry_run)
        exit(0)

    if args.extract:
        if not args.out:
            print("Error: --out is required with --extract")
            exit(1)
        i, j, N = args.extract
        extract_subdomain(hdf5_path, i, j, N, args.out)
        exit(0)

    if args.check_damping:
        if args.plot_damping_vs_k:
            rel_err = check_damping_stats(hdf5_path, plot=False)
            plot_damping_error_vs_k(hdf5_path, rel_err, geometry)
            exit(0)
        check_damping_stats(hdf5_path, plot=args.plot_damping_err, axis=args.axis, idx=args.idx)
        if not args.plot_damping_err:
            exit(0)

    if args.iheight is not None:
        heights = geometry.height
        idx = _find_closest_index(heights, args.iheight)
        print(f"Height {args.iheight:.4e} cm -> index {idx} (height: {heights[idx]:.4e} cm)")
        args.idx = idx

    if args.plot:
        print(f"Delta: {geometry.delta:.4e} cm ({geometry.delta/1e5:.4f} km)")
        geo = None if args.no_scale else geometry
        if args.axis == "xy":
            plot_slice_xy(hdf5_path, field=args.field, k=args.idx, log=args.log, sign=args.sign, quiver=args.quiver, geometry=geo, threshold=args.threshold, cmap=args.cmap)
        elif args.axis == "xz":
            plot_slice_xz(hdf5_path, field=args.field, j=args.idx, log=args.log, sign=args.sign, quiver=args.quiver, geometry=geo, threshold=args.threshold, cmap=args.cmap)
        elif args.axis == "yz":
            plot_slice_yz(hdf5_path, field=args.field, i=args.idx, log=args.log, sign=args.sign, quiver=args.quiver, geometry=geo, threshold=args.threshold, cmap=args.cmap)
        exit(0)

    if args.avg_vz:
        avg = average_z_velocity(hdf5_path)
        print(f"Average z bulk velocity: {avg}")
        exit(0)

    if args.z_profile:
        plot_z_profile(hdf5_path, field=args.field, geometry=geometry, show_minmax=args.z_profile_minmax)
        exit(0)

    if args.plot_all_profiles:
        plot_all_z_profiles(hdf5_path, fields=args.fields, geometry=geometry, show_minmax=args.z_profile_minmax, log=args.log)
        exit(0)

    print("=== Geometry ===")
    print(f"  Grid size : {geometry.N_x} × {geometry.N_y} × {geometry.N_z}")
    print(f"  delta     : {geometry.delta}")
    print(f"  height    : {geometry.height[:4]} … (first 4 values)")

    print("\n=== Frequency Grid ===")
    print(f"  N_frequencies : {freq_grid.N_frequencies}")
    print(f"  frequency     : {freq_grid.frequency[:4]} … (first 4 values)")

    print("\n=== Atom ===")
    print(f"  atomic_mass : {atom.atomic_mass}")
    print(f"  E_lower     : {atom.E_lower}")
    print(f"  E_upper     : {atom.E_upper}")
    print(f"  g_lower / g_upper : {atom.g_lower} / {atom.g_upper}")
    print(f"  Aul         : {atom.Aul}")
    print(f"  jl2 / ju2   : {atom.jl2} / {atom.ju2}")

    print(f"\n=== Point ({args.i}, {args.j}, {args.k}) ===")
    pt = read_point(hdf5_path, args.i, args.j, args.k)
    for (name, value) in pt.__dict__.items():
        print(f"  {name} : {value}")

    print("\n=== Field 'temperature' (shape, dtype, min, max) ===")
    import time
    _t0 = time.perf_counter()
    T = read_field(hdf5_path, "temperature")
    _elapsed = time.perf_counter() - _t0
    print(f"  shape   : {T.shape}")
    print(f"  dtype   : {T.dtype}")
    print(f"  range   : {T.min():.3f} – {T.max():.3f}")
    print(f"  elapsed : {_elapsed * 1e3:.3f} ms")

    print("\n=== read_field_slice examples (field: 'temperature') ===")

    test_i = min(geometry.N_x // 2, 4)
    test_j = min(geometry.N_y // 2, 4)
    test_k = min(10, geometry.N_z - 1)

    # fix i → y-z plane, shape (N_y, N_z)
    _t0 = time.perf_counter()
    slice_yz = read_field_slice(hdf5_path, "temperature", i=test_i)
    print(f"  fix i={test_i}  → shape {slice_yz.shape}  "
          f"range [{slice_yz.min():.3f}, {slice_yz.max():.3f}]  "
          f"({(time.perf_counter() - _t0) * 1e3:.3f} ms)")

    # fix j → x-z plane, shape (N_x, N_z)
    _t0 = time.perf_counter()
    slice_xz = read_field_slice(hdf5_path, "temperature", j=test_j)
    print(f"  fix j={test_j}  → shape {slice_xz.shape}  "
          f"range [{slice_xz.min():.3f}, {slice_xz.max():.3f}]  "
          f"({(time.perf_counter() - _t0) * 1e3:.3f} ms)")

    # fix k → x-y plane, shape (N_x, N_y)
    _t0 = time.perf_counter()
    slice_xy = read_field_slice(hdf5_path, "temperature", k=test_k)
    print(f"  fix k={test_k}  → shape {slice_xy.shape}  "
          f"range [{slice_xy.min():.3f}, {slice_xy.max():.3f}]  "
          f"({(time.perf_counter() - _t0) * 1e3:.3f} ms)")
