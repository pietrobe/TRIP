import h5py
import numpy as np
import argparse
import os


def resolve_cli_direction_indices(incl_value, azim_value, angular_grid):
    """
    Resolve CLI-provided direction arguments into valid dataset indices.

    Supported forms when angular grid metadata is available:
      1) explicit 0-based pair: (incl, azim)
      2) explicit 1-based pair: (incl, azim)
      3) flat direction index in --incl (0-based or 1-based), mapped via:
         incl = direction // N_azimuthal_angles, azim = direction % N_azimuthal_angles

    If no angular grid is available, values are returned unchanged.
    """
    if angular_grid is None:
        return incl_value, azim_value, "no angular-grid metadata; using raw indices"

    n_incl = int(angular_grid["N_inclination_angles"])
    n_azim = int(angular_grid["N_azimuthal_angles"])
    n_dirs = int(angular_grid["N_directions"])

    # 0-based pair (backward compatible behaviour)
    if 0 <= incl_value < n_incl and 0 <= azim_value < n_azim:
        return incl_value, azim_value, "interpreted as 0-based (incl, azim)"

    # 1-based pair
    if 1 <= incl_value <= n_incl and 1 <= azim_value <= n_azim:
        return incl_value - 1, azim_value - 1, "interpreted as 1-based (incl, azim)"

    # Flat direction index passed in --incl (0-based)
    if 0 <= incl_value < n_dirs:
        incl_i = incl_value // n_azim
        azim_i = incl_value % n_azim
        return incl_i, azim_i, "interpreted --incl as flat 0-based direction index"

    # Flat direction index passed in --incl (1-based)
    if 1 <= incl_value <= n_dirs:
        direction_i = incl_value - 1
        incl_i = direction_i // n_azim
        azim_i = direction_i % n_azim
        return incl_i, azim_i, "interpreted --incl as flat 1-based direction index"

    raise IndexError(
        "Invalid direction indices: "
        f"incl={incl_value}, azim={azim_value}. "
        f"Expected either (incl,azim) in [0..{n_incl - 1}]x[0..{n_azim - 1}] "
        f"or in [1..{n_incl}]x[1..{n_azim}], "
        f"or flat direction index in --incl within [0..{n_dirs - 1}] or [1..{n_dirs}]."
    )


def find_closest_direction_by_angle(incl_rad, azim_rad, angular_grid):
    """
    Search the angular grid for the direction whose (inclination, azimuth) angles
    are closest to the requested values (both expressed in radians).

    The angular distance metric is the great-circle separation on the unit sphere:

        d = arccos( sin(theta1)*sin(theta2)*cos(chi1-chi2) + cos(theta1)*cos(theta2) )

    where theta is the inclination (polar angle from z-axis) and chi is the azimuth.
    mu = cos(theta) is the direction cosine of the inclination.

    Parameters
    ----------
    incl_rad : float
        Requested inclination angle in radians  (theta).
        Can be derived from mu via  incl_rad = arccos(mu).
    azim_rad : float
        Requested azimuthal angle in radians (chi).
    angular_grid : dict
        Dict as returned by THDF_read_angular_grid_from_hdf5.

    Returns
    -------
    dict with keys:
        'incl_i'               -- 0-based inclination index of the closest direction
        'azim_i'               -- 0-based azimuthal index of the closest direction
        'direction_index'      -- 0-based flat direction index
        'incl_rad'             -- actual inclination angle theta (rad)
        'azim_rad'             -- actual azimuthal angle chi (rad)
        'incl_deg'             -- actual inclination angle theta (deg)
        'azim_deg'             -- actual azimuthal angle chi (deg)
        'mu'                   -- cos(theta) of the closest direction
        'angular_distance_rad' -- great-circle distance to the closest direction (rad)
    """
    incl_angles  = angular_grid["inclination_angles"]    # shape (N_incl,)
    azim_angles  = angular_grid["azimuthal_angles"]      # shape (N_azim,)
    incl_idx_arr = angular_grid["inclinations_indices"]  # shape (N_dirs,)
    azim_idx_arr = angular_grid["azimuthal_indices"]     # shape (N_dirs,)

    n_dirs = int(angular_grid["N_directions"])
    n_incl = int(angular_grid["N_inclination_angles"])
    n_azim = int(angular_grid["N_azimuthal_angles"])

    # inclinations_indices / azimuthal_indices are compact arrays (one entry per
    # unique angle, not one per direction) and may carry an internal offset
    # (e.g. [4,5,6,7] for 4 inclinations).  Normalise to 0-based and build the
    # full per-direction angle arrays via a Cartesian product
    # (inclination varies slowest, azimuth fastest — PORTA convention).
    incl_idx_local = incl_idx_arr[:n_incl] - incl_idx_arr[:n_incl].min()
    azim_idx_local = azim_idx_arr[:n_azim] - azim_idx_arr[:n_azim].min()

    # shape (N_dirs,): each inclination repeated n_azim times
    dir_incl = incl_angles[np.repeat(incl_idx_local, n_azim)]
    # shape (N_dirs,): azimuth pattern tiled n_incl times
    dir_azim = azim_angles[np.tile(azim_idx_local, n_incl)]

    # Great-circle distance:
    #   cos(d) = sin(theta1)*sin(theta2)*cos(phi1-phi2) + cos(theta1)*cos(theta2)
    cos_d = (
        np.sin(incl_rad) * np.sin(dir_incl) * np.cos(azim_rad - dir_azim)
        + np.cos(incl_rad) * np.cos(dir_incl)
    )
    # Clamp to [-1, 1] to guard against floating-point round-off
    cos_d = np.clip(cos_d, -1.0, 1.0)
    angular_distances = np.arccos(cos_d)  # (N_dirs,)

    best = int(np.argmin(angular_distances))
    # Recover the per-axis 0-based indices from the flat direction index
    best_incl_i = int(incl_idx_local[best // n_azim])
    best_azim_i = int(azim_idx_local[best  % n_azim])

    return {
        "incl_i":               best_incl_i,
        "azim_i":               best_azim_i,
        "direction_index":      best,
        "incl_rad":             float(dir_incl[best]),
        "azim_rad":             float(dir_azim[best]),
        "incl_deg":             float(np.degrees(dir_incl[best])),
        "azim_deg":             float(np.degrees(dir_azim[best])),
        "mu":                   float(np.cos(dir_incl[best])),
        "angular_distance_rad": float(angular_distances[best]),
    }


def THDF_read_field_from_hdf5(h5file, N_frequencies=None, x_i=0, y_i=0, incl_i=0, azimuth_i=0, N_azimuth=None):
    """
    Read a single surface output field (I,Q,U,V) for the given indices.

    Parameters:
      h5file: path or open h5py.File
      N_frequencies: optional int; if None inferred from dataset
      x_i, y_i, incl_i, azimuth_i: indices
      N_azimuth: optional int used to compute angular index as incl_i * N_azimuth + azimuth_i.
                 If omitted, falls back to the original C behaviour (incl_i * azimuth_i).

    Returns:
      dict {
        'index_i','index_j','index_incl','index_azimuth',
        'stokes_I','stokes_QI','stokes_UI','stokes_VI' } where stokes_* are numpy arrays dtype float64
    """
    close_when_done = False
    if isinstance(h5file, (str, bytes)):
        f = h5py.File(h5file, "r")
        close_when_done = True
    else:
        f = h5file

    try:
        if "/emergent_field" in f:
            grp = f["/emergent_field"]
            dsI = grp["emergent_I"]
            dsQ = grp["emergent_QI_pc"]
            dsU = grp["emergent_UI_pc"]
            dsV = grp["emergent_VI_pc"]

            if N_frequencies is None:
                N_frequencies = dsI.shape[4]

            stokes_I = np.array(
                dsI[x_i, y_i, incl_i, azimuth_i, :], dtype=np.float64).reshape((N_frequencies,))
            stokes_Q = np.array(
                dsQ[x_i, y_i, incl_i, azimuth_i, :], dtype=np.float64).reshape((N_frequencies,))
            stokes_U = np.array(
                dsU[x_i, y_i, incl_i, azimuth_i, :], dtype=np.float64).reshape((N_frequencies,))
            stokes_V = np.array(
                dsV[x_i, y_i, incl_i, azimuth_i, :], dtype=np.float64).reshape((N_frequencies,))
            index_incl    = incl_i
            index_azimuth = azimuth_i

        elif "/emergent_beams" in f:
            grp     = f["/emergent_beams"]
            dsI_grp = grp["beams_I"]
            dsQ_grp = grp["beams_QI_pc"]
            dsU_grp = grp["beams_UI_pc"]
            dsV_grp = grp["beams_VI_pc"]

            direction_names = sorted(list(dsI_grp.keys()))
            if "/beams_directions" in f:
                dgrp = f["/beams_directions"]
                if "mu" in dgrp and "chi" in dgrp:
                    mu = np.asarray(dgrp["mu"]).reshape(-1)
                    chi = np.asarray(dgrp["chi"]).reshape(-1)
                    ordered_names = [
                        f"mu_{m:.3f}_chi_{c:.3f}" for m, c in zip(mu, chi)]
                    if all(name in dsI_grp for name in ordered_names):
                        direction_names = ordered_names

            if len(direction_names) == 0:
                raise KeyError(
                    "No direction datasets found under '/emergent_beams/beams_I'")

            direction_i = incl_i
            if direction_i >= len(direction_names) and (direction_i - 1) < len(direction_names):
                direction_i = direction_i - 1
            if direction_i < 0 or direction_i >= len(direction_names):
                raise IndexError(
                    f"Direction index out of range: incl={incl_i}, available=[0..{len(direction_names)-1}]")

            direction_key = direction_names[direction_i]
            dsI = dsI_grp[direction_key]
            dsQ = dsQ_grp[direction_key]
            dsU = dsU_grp[direction_key]
            dsV = dsV_grp[direction_key]

            if N_frequencies is None:
                N_frequencies = dsI.shape[2]

            stokes_I = np.array(dsI[x_i, y_i, :], dtype=np.float64).reshape((N_frequencies,))
            stokes_Q = np.array(dsQ[x_i, y_i, :], dtype=np.float64).reshape((N_frequencies,))
            stokes_U = np.array(dsU[x_i, y_i, :], dtype=np.float64).reshape((N_frequencies,))
            stokes_V = np.array(dsV[x_i, y_i, :], dtype=np.float64).reshape((N_frequencies,))
            index_incl    = direction_i
            index_azimuth = 0

        else:
            raise KeyError(
                "Neither '/emergent_field' nor '/emergent_beams' group found")

        return {
            "index_i":        x_i,
            "index_j":        y_i,
            "index_incl":     index_incl,
            "index_azimuth":  index_azimuth,
            "stokes_I":       stokes_I,
            "stokes_QI":      stokes_Q,
            "stokes_UI":      stokes_U,
            "stokes_VI":      stokes_V,
        }
    finally:
        if close_when_done:
            f.close()


def THDF_read_angular_grid_from_hdf5(h5file):
    """
    Read the /emergent_angular_grid group and return a dict:
      {
        'N_directions': int,
        'N_inclination_angles': int,
        'N_azimuthal_angles': int,
        'inclination_angles': np.ndarray(dtype=float64),
        'azimuthal_angles': np.ndarray(dtype=float64),
        'inclinations_indices': np.ndarray(dtype=int),
        'azimuthal_indices': np.ndarray(dtype=int)
      }

    `h5file` may be a path (str/Path) or an open h5py.File.
    """
    close_when_done = False
    if isinstance(h5file, (str, bytes)):
        f = h5py.File(h5file, "r")
        close_when_done = True
    else:
        f = h5file

    try:
        grp = f["/emergent_angular_grid"]
    except KeyError:
        if close_when_done:
            f.close()
        raise KeyError("Group '/emergent_angular_grid' not found")

    # Read N_directions (if present) otherwise infer from arrays
    if "N_directions" in grp:
        N_val = grp["N_directions"][()]
        N = int(N_val.item() if hasattr(N_val, 'item')
                else N_val[0] if hasattr(N_val, '__getitem__') else N_val)
    else:
        if "inclination_angles" in grp:
            N = int(grp["inclination_angles"].shape[0])
        else:
            if close_when_done:
                f.close()
            raise KeyError(
                "Neither 'N_directions' nor 'inclination_angles' found in group")

    N_incl_val = grp["N_inclination_angles"][()]
    N_inclination_angles = int(N_incl_val.item() if hasattr(
        N_incl_val, 'item') else N_incl_val[0] if hasattr(N_incl_val, '__getitem__') else N_incl_val)
    N_azim_val = grp["N_azimuthal_angles"][()]
    N_azimuthal_angles = int(N_azim_val.item() if hasattr(
        N_azim_val, 'item') else N_azim_val[0] if hasattr(N_azim_val, '__getitem__') else N_azim_val)

    inclination_angles = np.array(grp["inclination_angles"], dtype=np.float64)
    azimuthal_angles   = np.array(grp["azimuthal_angles"],   dtype=np.float64)
    incl_indices       = np.array(grp["inclinations_indices"], dtype=np.int32)
    azim_indices       = np.array(grp["azimuthal_indices"],    dtype=np.int32)

    if close_when_done:
        f.close()

    return {
        "N_directions":          N,
        "inclination_angles":    inclination_angles,
        "azimuthal_angles":      azimuthal_angles,
        "inclinations_indices":  incl_indices,
        "azimuthal_indices":     azim_indices,
        "N_inclination_angles":  N_inclination_angles,
        "N_azimuthal_angles":    N_azimuthal_angles,
    }


def THDF_read_frequencies_grid_from_hdf5(h5file):
    """
    Read the /frequencies_grid group and return a dict:
      {'N_frequencies': int, 'frequencies': np.ndarray(dtype=float64)}

    `h5file` may be a path (str) or an open h5py.File.
    """
    close_when_done = False
    if isinstance(h5file, str):
        f = h5py.File(h5file, "r")
        close_when_done = True
    else:
        f = h5file

    try:
        grp = f["/frequencies_grid"]
    except KeyError:
        if close_when_done:
            f.close()
        raise KeyError("Group '/frequencies_grid' not found")

    if "N_frequencies" in grp:
        N_val = grp["N_frequencies"][()]
        N = int(np.asarray(N_val).reshape(-1)[0])
    else:
        N = int(grp["frequencies"].shape[0])

    frequencies = np.array(grp["frequencies"], dtype=np.float64)

    if close_when_done:
        f.close()

    return {"N_frequencies": N, "frequencies": frequencies}


def CSV_read_angular_grid(csv_path):
    """
    Read angular_grid_0_0.csv with columns:
        theta_i, chi_i, theta, mu, chi

    Returns a dict:
        {
          'theta_i' : np.ndarray(int),
          'chi_i'   : np.ndarray(int),
          'theta'   : np.ndarray(float64),   # inclination in radians
          'mu'      : np.ndarray(float64),   # cos(theta)
          'chi'     : np.ndarray(float64),   # azimuth in radians
        }
    """
    data = np.genfromtxt(csv_path, delimiter=",", names=True, dtype=None,
                         encoding="utf-8")
    return {
        "theta_i": data["theta_i"].astype(np.int32),
        "chi_i":   data["chi_i"].astype(np.int32),
        "theta":   data["theta"].astype(np.float64),
        "mu":      data["mu"].astype(np.float64),
        "chi":     data["chi"].astype(np.float64),
    }


def CSV_find_closest_direction(theta_rad, chi_rad, csv_ag):
    """
    Find the row in a CSV angular grid (as returned by CSV_read_angular_grid)
    that is closest in great-circle distance to (theta_rad, chi_rad).

    Returns (dir_index, row_dict) where dir_index is the 0-based row index and
    row_dict contains the matched theta, mu, chi values.
    """
    thetas = csv_ag["theta"]
    chis   = csv_ag["chi"]

    cos_d = (
        np.sin(theta_rad) * np.sin(thetas) * np.cos(chi_rad - chis)
        + np.cos(theta_rad) * np.cos(thetas)
    )
    cos_d = np.clip(cos_d, -1.0, 1.0)
    distances = np.arccos(cos_d)

    best = int(np.argmin(distances))
    return best, {
        "theta": float(thetas[best]),
        "mu":    float(csv_ag["mu"][best]),
        "chi":   float(chis[best]),
        "angular_distance_rad": float(distances[best]),
    }


def CSV_read_stokes_profiles(csv_path, dir_index):
    """
    Read profiles_<i>_<j>.csv and extract the four Stokes rows for direction
    dir_index.

    Layout: groups of 4 consecutive rows per direction → [I, Q, U, V].
    Each row contains N_freq values (the frequency axis).

    Trailing commas on each line (a known output artefact) are stripped before
    parsing so the file is treated as a proper CSV.

    Returns a dict:
        {
          'stokes_I'  : np.ndarray(float64),
          'stokes_QI' : np.ndarray(float64),
          'stokes_UI' : np.ndarray(float64),
          'stokes_VI' : np.ndarray(float64),
        }
    """
    # Read raw text, strip trailing commas, then parse via numpy
    with open(csv_path, "r", encoding="utf-8") as fh:
        cleaned = "\n".join(line.rstrip().rstrip(",") for line in fh)

    import io
    raw = np.genfromtxt(io.StringIO(cleaned), delimiter=",", dtype=np.float64)

    # Skip an all-NaN header row if present
    if raw.ndim > 1 and np.all(np.isnan(raw[0])):
        raw = raw[1:]

    row0 = 4 * dir_index
    if row0 + 3 >= len(raw):
        raise IndexError(
            f"dir_index={dir_index} requires rows {row0}..{row0+3} but "
            f"the profiles file has only {len(raw)} rows."
        )

    return {
        "stokes_I":  raw[row0    ],
        "stokes_QI": raw[row0 + 1],
        "stokes_UI": raw[row0 + 2],
        "stokes_VI": raw[row0 + 3],
    }


def read_map_data(h5file, mu, chi, freq_idx=None):
    """
    Read the emergent Stokes field over all spatial points (i, j) for the
    direction closest to the requested (mu, chi) and the requested frequency.

    Parameters
    ----------
    h5file : str or h5py.File
        Path to the HDF5 file or an already-open file object.
    mu : float
        Direction cosine of the inclination: mu = cos(theta).  Must be in [-1, 1].
    chi : float
        Azimuthal angle in radians.
    freq_idx : int or None
        Index along the frequency axis to extract.
        - int  → returns 2-D arrays of shape (N_x, N_y) for each Stokes parameter.
        - None → returns 3-D arrays of shape (N_x, N_y, N_freq).

    Returns
    -------
    dict with keys:
        'incl_i'      – 0-based inclination index used
        'azim_i'      – 0-based azimuthal index used
        'incl_rad'    – actual inclination angle theta (rad)
        'azim_rad'    – actual azimuthal angle chi (rad)
        'mu'          – cos(theta) of the matched direction
        'angular_distance_rad' – great-circle distance to the matched direction
        'stokes_I'    – np.ndarray  shape (N_x, N_y)        if freq_idx is int
                                     shape (N_x, N_y, N_freq) if freq_idx is None
        'stokes_QI'   – same shape as stokes_I
        'stokes_UI'   – same shape as stokes_I
        'stokes_VI'   – same shape as stokes_I
        'freq_idx'    – the freq_idx argument echoed back (int or None)
    """
    close_when_done = False
    if isinstance(h5file, (str, bytes)):
        f = h5py.File(h5file, "r")
        close_when_done = True
    else:
        f = h5file

    try:
        # ------------------------------------------------------------------ #
        # 1. Read angular grid and find the closest direction                 #
        # ------------------------------------------------------------------ #
        angular_grid = THDF_read_angular_grid_from_hdf5(f)
        theta_rad    = float(np.arccos(np.clip(mu, -1.0, 1.0)))
        closest      = find_closest_direction_by_angle(theta_rad, chi, angular_grid)
        incl_i       = closest["incl_i"]
        azim_i       = closest["azim_i"]

        # ------------------------------------------------------------------ #
        # 2. Read the full spatial map                                        #
        # ------------------------------------------------------------------ #
        if "/emergent_field" in f:
            grp = f["/emergent_field"]
            dsI = grp["emergent_I"]        # shape: (N_x, N_y, N_incl, N_azim, N_freq)
            dsQ = grp["emergent_QI_pc"]
            dsU = grp["emergent_UI_pc"]
            dsV = grp["emergent_VI_pc"]

            if freq_idx is None:
                # Full frequency axis → (N_x, N_y, N_freq)
                I = np.array(dsI[:, :, incl_i, azim_i, :], dtype=np.float64)
                Q = np.array(dsQ[:, :, incl_i, azim_i, :], dtype=np.float64)
                U = np.array(dsU[:, :, incl_i, azim_i, :], dtype=np.float64)
                V = np.array(dsV[:, :, incl_i, azim_i, :], dtype=np.float64)
            else:
                # Single frequency → (N_x, N_y)
                I = np.array(dsI[:, :, incl_i, azim_i, freq_idx], dtype=np.float64)
                Q = np.array(dsQ[:, :, incl_i, azim_i, freq_idx], dtype=np.float64)
                U = np.array(dsU[:, :, incl_i, azim_i, freq_idx], dtype=np.float64)
                V = np.array(dsV[:, :, incl_i, azim_i, freq_idx], dtype=np.float64)

        elif "/emergent_beams" in f:
            grp     = f["/emergent_beams"]
            dsI_grp = grp["beams_I"]       # each direction is a sub-dataset
            dsQ_grp = grp["beams_QI_pc"]
            dsU_grp = grp["beams_UI_pc"]
            dsV_grp = grp["beams_VI_pc"]

            # Reconstruct the direction key list (same logic as in THDF_read_field_from_hdf5)
            direction_names = sorted(list(dsI_grp.keys()))
            if "/beams_directions" in f:
                dgrp = f["/beams_directions"]
                if "mu" in dgrp and "chi" in dgrp:
                    mu_arr  = np.asarray(dgrp["mu"]).reshape(-1)
                    chi_arr = np.asarray(dgrp["chi"]).reshape(-1)
                    ordered = [f"mu_{m:.3f}_chi_{c:.3f}"
                               for m, c in zip(mu_arr, chi_arr)]
                    if all(name in dsI_grp for name in ordered):
                        direction_names = ordered

            # Flat direction index from (incl_i, azim_i)
            n_azim   = int(angular_grid["N_azimuthal_angles"])
            flat_idx = incl_i * n_azim + azim_i
            if flat_idx >= len(direction_names):
                flat_idx = len(direction_names) - 1

            key  = direction_names[flat_idx]
            # sub-dataset shape: (N_x, N_y, N_freq)
            dsI_d = dsI_grp[key]
            dsQ_d = dsQ_grp[key]
            dsU_d = dsU_grp[key]
            dsV_d = dsV_grp[key]

            if freq_idx is None:
                I = np.array(dsI_d[:, :, :], dtype=np.float64)
                Q = np.array(dsQ_d[:, :, :], dtype=np.float64)
                U = np.array(dsU_d[:, :, :], dtype=np.float64)
                V = np.array(dsV_d[:, :, :], dtype=np.float64)
            else:
                I = np.array(dsI_d[:, :, freq_idx], dtype=np.float64)
                Q = np.array(dsQ_d[:, :, freq_idx], dtype=np.float64)
                U = np.array(dsU_d[:, :, freq_idx], dtype=np.float64)
                V = np.array(dsV_d[:, :, freq_idx], dtype=np.float64)

        else:
            raise KeyError(
                "Neither '/emergent_field' nor '/emergent_beams' group found")

    finally:
        if close_when_done:
            f.close()

    return {
        "incl_i":               incl_i,
        "azim_i":               azim_i,
        "incl_rad":             closest["incl_rad"],
        "azim_rad":             closest["azim_rad"],
        "mu":                   closest["mu"],
        "angular_distance_rad": closest["angular_distance_rad"],
        "stokes_I":             I,
        "stokes_QI":            Q,
        "stokes_UI":            U,
        "stokes_VI":            V,
        "freq_idx":             freq_idx,
    }


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Read HDF5 demo outputs")
    parser.add_argument("file", nargs="?", default="output_field_mpi.h5",
                        help="HDF5 file to read (default: output_field_mpi.h5)")
    parser.add_argument("--plot", action="store_true", default=True,
                        help="Plot Stokes I,Q,U,V vs frequency")
    parser.add_argument("--save", help="If set, save plot to this filename (e.g. plot.png)")
    parser.add_argument("--x",    type=int,   default=1, help="x index (default 1)")
    parser.add_argument("--y",    type=int,   default=2, help="y index (default 2)")
    parser.add_argument("--incl", type=int,   default=2, help="inclination index (default 2)")
    parser.add_argument("--azim", type=int,   default=2, help="azimuth index (default 2)")
    # --- angle-based direction search (values in radians) ---
    parser.add_argument(
        "--incl-rad", type=float, default=None,
        help="Requested inclination angle theta in radians.  When supplied together "
             "with --azim-rad, the closest grid direction is used instead of "
             "the --incl / --azim index arguments.  Mutually exclusive with --mu.")
    parser.add_argument(
        "--azim-rad", type=float, default=None,
        help="Requested azimuthal angle chi in radians. Used with --incl-rad or --mu. "
             "If omitted when --incl-rad/--mu is given, the azimuthal angle is taken "
             "from the grid at the --azim index.")
    parser.add_argument(
        "--mu", type=float, default=None,
        help="Requested direction cosine mu = cos(theta) of the inclination.  "
             "Must be in [-1, 1].  When supplied together with --azim-rad, the closest "
             "grid direction is used.  Mutually exclusive with --incl-rad.")
    parser.add_argument(
        "--lincl", type=int, nargs="+", metavar="INCL_IDX", default=None,
        help="List of 0-based inclination indices to plot together for a fixed "
             "azimuthal angle (given via --azim-rad or --azim). Each inclination "
             "is drawn with a distinct color; legend shown only in the first panel. "
             "Mutually exclusive with --incl-rad, --mu.")
    parser.add_argument(
        "--compare", action="store_true", default=False,
        help="Overlay reference data from CSV files found in the same directory "
             "as the HDF5 file.  Looks for 'angular_grid_0_0.csv' (to match the "
             "closest direction) and 'profiles_<x>_<y>.csv' (Stokes I,Q,U,V rows "
             "at 4*dir_index : 4*dir_index+4).  Works in both single-direction "
             "and --lincl modes.")
    parser.add_argument(
        "--compare-dir", default=None, metavar="DIR",
        help="Directory containing the CSV comparison files "
             "(angular_grid_0_0.csv and profiles_<x>_<y>.csv). "
             "Defaults to the directory of the HDF5 file when --compare is used.")
    parser.add_argument(
        "--map", action="store_true", default=False,
        help="Plot 2-D spatial maps (I, Q/I, U/I, V/I) over all (x, y) points "
             "for the direction resolved from --mu/--incl-rad/--incl and --azim/--azim-rad. "
             "Use --freq-idx to select a single frequency slice (default: 0).")
    parser.add_argument(
        "--freq-idx", type=int, default=0, metavar="F",
        help="Frequency index for --map mode (default: 0).  "
             "Ignored when --map is not set.")
    args = parser.parse_args()
    file_name = args.file

    output_freq_grid = THDF_read_frequencies_grid_from_hdf5(file_name)
    print("Number of frequencies:", output_freq_grid["N_frequencies"])
    print("Frequencies (Hz):", output_freq_grid["frequencies"])

    try:
        angular_grid = THDF_read_angular_grid_from_hdf5(file_name)
    except KeyError:
        angular_grid = None
        N_azimuth    = None
        print("Warning: group '/emergent_angular_grid' not found; "
              "continuing with explicit --incl/--azim indices.")
    else:
        azimithal_angles   = angular_grid["azimuthal_angles"]
        inclination_angles = angular_grid["inclination_angles"]

        az_len   = len(azimithal_angles)
        incl_len = len(inclination_angles)

        print("Number of directions:", angular_grid["N_directions"])
        print("Inclination angles (rad):", inclination_angles)
        print("Azimuthal angles (rad):  ", azimithal_angles)
        print("Inclination indices:", angular_grid["inclinations_indices"][:incl_len])
        print("Azimuthal indices:  ", angular_grid["azimuthal_indices"][:az_len])
        print("N_inclination_angles:", angular_grid["N_inclination_angles"])
        print("N_azimuthal_angles:  ", angular_grid["N_azimuthal_angles"])
        print("N_directions:        ", angular_grid["N_directions"])

        N_azimuth = angular_grid["N_azimuthal_angles"]

    # -----------------------------------------------------------------------
    # Validate mutually exclusive top-level modes
    # -----------------------------------------------------------------------
    if args.lincl is not None and (args.incl_rad is not None or args.mu is not None):
        parser.error("--lincl is mutually exclusive with --incl-rad and --mu.")

    # -----------------------------------------------------------------------
    # Pre-load CSV comparison data (shared by both plot branches)
    # -----------------------------------------------------------------------
    h5_dir      = os.path.dirname(os.path.abspath(file_name))
    compare_dir = os.path.abspath(args.compare_dir) if args.compare_dir else h5_dir

    def _load_compare_data(theta_rad, chi_rad, x_i, y_i):
        """
        Load the CSV reference profiles for the direction closest to
        (theta_rad, chi_rad) at spatial point (x_i, y_i).

        Returns a dict:
            {
              'stokes_I', 'stokes_QI', 'stokes_UI', 'stokes_VI',
              'dir_index', 'theta', 'mu', 'chi', 'angular_distance_rad',
            }
        or None if the required files are not found.
        """
        ag_csv   = os.path.join(compare_dir, "angular_grid_0_0.csv")
        prof_csv = os.path.join(compare_dir, f"profiles_{x_i}_{y_i}.csv")

        if not os.path.isfile(ag_csv):
            print(f"[compare] angular grid CSV not found: {ag_csv}")
            return None
        if not os.path.isfile(prof_csv):
            print(f"[compare] profiles CSV not found: {prof_csv}")
            return None

        csv_ag                = CSV_read_angular_grid(ag_csv)
        dir_index, match_info = CSV_find_closest_direction(
            theta_rad, chi_rad, csv_ag)
        profiles              = CSV_read_stokes_profiles(prof_csv, dir_index)

        print(
            f"[compare] angular_grid: {ag_csv}\n"
            f"[compare] profiles:     {prof_csv}\n"
            f"[compare] closest CSV direction: dir_index={dir_index}, "
            f"theta={match_info['theta']:.4f} rad, "
            f"mu={match_info['mu']:.4f}, "
            f"chi={match_info['chi']:.4f} rad, "
            f"angular distance={np.degrees(match_info['angular_distance_rad']):.4f}°"
        )

        return {**profiles, **match_info, "dir_index": dir_index}

    # -----------------------------------------------------------------------
    # Helper: build a closest_info-style dict from 0-based axis indices
    # -----------------------------------------------------------------------
    def _info_from_indices(incl_i, azim_i, agrid):
        ia = agrid["inclination_angles"]
        aa = agrid["azimuthal_angles"]
        ir = float(ia[incl_i])
        ar = float(aa[azim_i])
        return {
            "incl_i":               incl_i,
            "azim_i":               azim_i,
            "direction_index":      incl_i * int(agrid["N_azimuthal_angles"]) + azim_i,
            "incl_rad":             ir,
            "azim_rad":             ar,
            "incl_deg":             float(np.degrees(ir)),
            "azim_deg":             float(np.degrees(ar)),
            "mu":                   float(np.cos(ir)),
            "angular_distance_rad": 0.0,
        }

    # -----------------------------------------------------------------------
    # Helper: build a legend label from a closest_info dict
    # -----------------------------------------------------------------------
    def _legend_label(info, is_nearest=False):
        label = (
            rf"$\theta$={info['incl_rad']:.4f},  "
            rf"$\mu$={info['mu']:.4f}"
        )
        if is_nearest and info["angular_distance_rad"] > 0.0:
            label += rf",  $\Delta$={np.degrees(info['angular_distance_rad']):.4f}°"
        return label

    # -----------------------------------------------------------------------
    # Helper: resolve the azimuthal angle in radians
    # -----------------------------------------------------------------------
    def _resolve_azim_rad(args, angular_grid, parser):
        if args.azim_rad is not None:
            return args.azim_rad, f"chi={args.azim_rad:.6f} rad (from --azim-rad)"
        if angular_grid is None:
            parser.error("Angle-based search requires angular-grid metadata.")
        n_azim = int(angular_grid["N_azimuthal_angles"])
        idx = args.azim
        if not (0 <= idx < n_azim):
            parser.error(
                f"--azim {idx} out of range for grid with N_azimuthal_angles={n_azim}.")
        val = float(angular_grid["azimuthal_angles"][idx])
        return val, f"chi={val:.6f} rad (from --azim index {idx})"

    # -----------------------------------------------------------------------
    # Validate mutually exclusive angle-search options (single-direction)
    # -----------------------------------------------------------------------
    if args.incl_rad is not None and args.mu is not None:
        parser.error("--incl-rad and --mu are mutually exclusive.")

    # -----------------------------------------------------------------------
    # Branch: --lincl  (multi-inclination mode)
    # -----------------------------------------------------------------------
    if args.lincl is not None:
        if angular_grid is None:
            parser.error("--lincl requires angular-grid metadata in the HDF5 file.")

        n_incl_grid = int(angular_grid["N_inclination_angles"])
        for idx in args.lincl:
            if not (0 <= idx < n_incl_grid):
                parser.error(
                    f"--lincl index {idx} out of range [0..{n_incl_grid - 1}].")

        azim_rad_fixed, azim_note = _resolve_azim_rad(args, angular_grid, parser)
        azim_resolved = int(np.argmin(
            np.abs(angular_grid["azimuthal_angles"] - azim_rad_fixed)))

        print(f"--lincl mode: fixed {azim_note}, resolved to azim_i={azim_resolved}")

        directions = []
        for incl_i in args.lincl:
            info  = _info_from_indices(incl_i, azim_resolved, angular_grid)
            field = THDF_read_field_from_hdf5(
                file_name,
                N_frequencies=output_freq_grid["N_frequencies"],
                x_i=args.x,
                y_i=args.y,
                incl_i=incl_i,
                azimuth_i=azim_resolved,
                N_azimuth=N_azimuth,
            )
            directions.append((info, field))
            print(f"  incl_i={incl_i}: "
                  f"theta={info['incl_rad']:.4f} rad ({info['incl_deg']:.2f}°), "
                  f"mu={info['mu']:.4f}")

        if args.plot:
            try:
                import matplotlib.pyplot as plt
            except Exception:
                print("matplotlib not available; cannot plot")
            else:
                from wavelength_utils import freq_to_angstrom, vacuum_to_air

                freqs       = output_freq_grid["frequencies"]
                wl_angstrom = freq_to_angstrom(freqs)
                wl_air      = vacuum_to_air(wl_angstrom)

                fig, axs = plt.subplots(2, 2, figsize=(10, 8), sharex=True)
                ax = axs.ravel()

                plot_keys = ["stokes_I", "stokes_QI", "stokes_UI", "stokes_VI"]
                ylabels   = ["Stokes I", "Stokes Q/I (%)",
                             "Stokes U/I (%)", "Stokes V/I (%)"]
                legend_handles = []

                for k, (info, field) in enumerate(directions):
                    color = f"C{k}"
                    label = _legend_label(info)
                    for panel, key in enumerate(plot_keys):
                        line, = ax[panel].plot(
                            wl_air, field[key],
                            color=color, marker=".", label=label)
                        if panel == 0:
                            legend_handles.append(line)

                    # --compare overlay for this inclination
                    if args.compare:
                        cmp = _load_compare_data(
                            info["incl_rad"], info["azim_rad"],
                            args.x, args.y)
                        if cmp is not None:
                            cmp_label = (
                                rf"ref $\theta$={cmp['theta']:.4f}, "
                                rf"$\mu$={cmp['mu']:.4f}"
                            )
                            cmp_keys = ["stokes_I", "stokes_QI",
                                        "stokes_UI", "stokes_VI"]
                            for panel, key in enumerate(cmp_keys):
                                ref_line, = ax[panel].plot(
                                    wl_air, cmp[key],
                                    color="darkorange", linestyle="--",
                                    marker="x", markersize=6,
                                    markevery=4, label=cmp_label)
                            legend_handles.append(ref_line)

                for panel, ylabel in enumerate(ylabels):
                    ax[panel].set_ylabel(ylabel)
                    ax[panel].grid(True)
                    if panel >= 2:
                        ax[panel].set_xlabel("Wavelength (Å)")

                ax[0].legend(
                    handles=legend_handles,
                    loc="upper right",
                    fontsize=8,
                    framealpha=0.85,
                    handlelength=1.2,
                )

                abs_path   = os.path.abspath(file_name)
                path_parts = abs_path.split(os.sep)[-3:]
                path_title = "/".join(path_parts)
                incl_str        = ",".join(str(i) for i in args.lincl)
                actual_azim_rad = float(angular_grid["azimuthal_angles"][azim_resolved])

                title_text = fig.suptitle(
                    f"{path_title}\nStokes vs Wavelength  "
                    rf"(x={args.x}, y={args.y}, incl=[{incl_str}], "
                    rf"$\chi$={actual_azim_rad:.4f})",
                    fontsize=14, fontweight="bold", color="darkblue")
                title_text.set_bbox(
                    dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))
                fig.tight_layout()
                fig.subplots_adjust(top=0.85)

                if args.save:
                    fig.savefig(args.save)
                    print(f"Saved plot to {args.save}")
                else:
                    plt.show()

    # -----------------------------------------------------------------------
    # Branch: single-direction mode  (--incl-rad / --mu / --incl)
    # -----------------------------------------------------------------------
    else:
        if args.mu is not None:
            if not (-1.0 <= args.mu <= 1.0):
                parser.error(f"--mu must be in [-1, 1], got {args.mu}")
            search_incl_rad = float(np.arccos(args.mu))
            angle_search    = True
        elif args.incl_rad is not None:
            search_incl_rad = args.incl_rad
            angle_search    = True
        else:
            angle_search = False

        closest_info = None

        if angle_search:
            if angular_grid is None:
                print("Warning: angle-based search requires angular-grid metadata; "
                      "falling back to index-based resolution.")
                incl_resolved, azim_resolved, resolution_note = \
                    resolve_cli_direction_indices(args.incl, args.azim, angular_grid)
            else:
                search_azim_rad, azim_source_note = _resolve_azim_rad(
                    args, angular_grid, parser)

                closest_info  = find_closest_direction_by_angle(
                    search_incl_rad, search_azim_rad, angular_grid)
                incl_resolved = closest_info["incl_i"]
                azim_resolved = closest_info["azim_i"]
                requested_str = (
                    f"mu={args.mu:.6f} -> theta={search_incl_rad:.6f} rad"
                    if args.mu is not None
                    else f"theta={args.incl_rad:.6f} rad"
                )
                resolution_note = (
                    f"closest direction to ({requested_str}, {azim_source_note}): "
                    f"grid point "
                    f"(theta={closest_info['incl_rad']:.6f} rad = "
                    f"{closest_info['incl_deg']:.3f}°, "
                    f"mu={closest_info['mu']:.6f}, "
                    f"chi={closest_info['azim_rad']:.6f} rad = "
                    f"{closest_info['azim_deg']:.3f}°), "
                    f"flat index {closest_info['direction_index']}, "
                    f"angular distance "
                    f"{np.degrees(closest_info['angular_distance_rad']):.4f}°"
                )

        else:
            if args.azim_rad is not None:
                parser.error(
                    "--azim-rad must be used together with --incl-rad or --mu.")
            incl_resolved, azim_resolved, resolution_note = \
                resolve_cli_direction_indices(args.incl, args.azim, angular_grid)

            if angular_grid is not None:
                closest_info = _info_from_indices(
                    incl_resolved, azim_resolved, angular_grid)

        print(
            "Resolved indices:",
            f"incl={incl_resolved}, azim={azim_resolved}",
            f"({resolution_note})",
        )

        output_filed = THDF_read_field_from_hdf5(
            file_name,
            N_frequencies=output_freq_grid["N_frequencies"],
            x_i=args.x,
            y_i=args.y,
            incl_i=incl_resolved,
            azimuth_i=azim_resolved,
            N_azimuth=N_azimuth,
        )

        print("Output field:")
        print("Stokes I:", output_filed["stokes_I"])
        print("Stokes Q:", output_filed["stokes_QI"])
        print("Stokes U:", output_filed["stokes_UI"])
        print("Stokes V:", output_filed["stokes_VI"])

        if args.plot:
            try:
                import matplotlib.pyplot as plt
            except Exception:
                print("matplotlib not available; cannot plot")
            else:
                from wavelength_utils import freq_to_angstrom, vacuum_to_air

                freqs       = output_freq_grid["frequencies"]
                wl_angstrom = freq_to_angstrom(freqs)
                wl_air      = vacuum_to_air(wl_angstrom)

                fig, axs = plt.subplots(2, 2, figsize=(10, 8), sharex=True)
                ax = axs.ravel()

                DIR_COLOR = "C0"

                if closest_info is not None:
                    legend_label = _legend_label(closest_info, is_nearest=angle_search)
                else:
                    legend_label = (
                        f"incl idx={incl_resolved},  azim idx={azim_resolved}"
                    )

                plot_specs = [
                    ("Stokes I",       output_filed["stokes_I"]),
                    ("Stokes Q/I (%)", output_filed["stokes_QI"]),
                    ("Stokes U/I (%)", output_filed["stokes_UI"]),
                    ("Stokes V/I (%)", output_filed["stokes_VI"]),
                ]

                legend_handles = []
                for panel, (ylabel, data) in enumerate(plot_specs):
                    line, = ax[panel].plot(
                        wl_air, data, color=DIR_COLOR,
                        marker=".", label=legend_label)
                    ax[panel].set_ylabel(ylabel)
                    ax[panel].grid(True)
                    if panel == 0:
                        legend_handles.append(line)
                    if panel >= 2:
                        ax[panel].set_xlabel("Wavelength (Å)")

                # --compare: overlay CSV reference profiles
                if args.compare:
                    theta_ref = (closest_info["incl_rad"]
                                 if closest_info is not None
                                 else float(angular_grid["inclination_angles"][incl_resolved])
                                 if angular_grid is not None else None)
                    chi_ref   = (closest_info["azim_rad"]
                                 if closest_info is not None
                                 else float(angular_grid["azimuthal_angles"][azim_resolved])
                                 if angular_grid is not None else None)
                    if theta_ref is None:
                        print("[compare] Cannot resolve angles without angular-grid metadata.")
                    else:
                        cmp = _load_compare_data(theta_ref, chi_ref, args.x, args.y)
                        if cmp is not None:
                            cmp_label = (
                                rf"ref $\theta$={cmp['theta']:.4f}, "
                                rf"$\mu$={cmp['mu']:.4f}"
                            )
                            cmp_keys = ["stokes_I", "stokes_QI",
                                        "stokes_UI", "stokes_VI"]
                            n_csv = len(cmp["stokes_I"])
                            n_h5  = len(wl_air)
                            if n_csv != n_h5:
                                wl_csv = np.linspace(wl_air[0], wl_air[-1], n_csv)
                                cmp_interp = {
                                    k: np.interp(wl_air, wl_csv, cmp[k])
                                    for k in cmp_keys
                                }
                                print(
                                    f"[compare] CSV has {n_csv} pts, HDF5 has {n_h5} pts "
                                    f"— interpolated onto HDF5 wavelength grid."
                                )
                            else:
                                cmp_interp = {k: cmp[k] for k in cmp_keys}
                            for panel, key in enumerate(cmp_keys):
                                ref_line, = ax[panel].plot(
                                    wl_air, cmp_interp[key],
                                    color="darkorange", linestyle="--",
                                    marker="x", markersize=6,
                                    markevery=4, label=cmp_label)
                            legend_handles.append(ref_line)

                ax[0].legend(
                    handles=legend_handles,
                    loc="upper right",
                    fontsize=8,
                    framealpha=0.85,
                    handlelength=1.2,
                )

                abs_path   = os.path.abspath(file_name)
                path_parts = abs_path.split(os.sep)[-3:]
                path_title = "/".join(path_parts)

                # Always show the actual grid azimuth angle in the title
                if closest_info is not None:
                    actual_azim_rad = closest_info["azim_rad"]
                    azim_title = rf"$\chi$={actual_azim_rad:.4f}"
                elif angular_grid is not None:
                    actual_azim_rad = float(
                        angular_grid["azimuthal_angles"][azim_resolved])
                    azim_title = rf"$\chi$={actual_azim_rad:.4f}"
                else:
                    azim_title = f"azim={azim_resolved}"

                title_text = fig.suptitle(
                    f"{path_title}\nStokes vs Wavelength  "
                    f"(x={args.x}, y={args.y}, "
                    f"incl={incl_resolved}, {azim_title})",
                    fontsize=16, fontweight="bold", color="darkblue")
                title_text.set_bbox(
                    dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))
                fig.tight_layout()
                fig.subplots_adjust(top=0.85)

                if args.save:
                    fig.savefig(args.save)
                    print(f"Saved plot to {args.save}")
                else:
                    plt.show()

    # -----------------------------------------------------------------------
    # Branch: --map  (2-D spatial map at a fixed frequency)
    # -----------------------------------------------------------------------
    if args.map:
        try:
            import matplotlib.pyplot as plt
            import matplotlib.colors as mcolors
        except Exception:
            print("matplotlib not available; cannot plot map")
        else:
            from wavelength_utils import freq_to_angstrom, vacuum_to_air

            # Resolve direction from the same CLI options used in the other modes
            if angular_grid is not None:
                if args.mu is not None:
                    map_theta = float(np.arccos(np.clip(args.mu, -1.0, 1.0)))
                    map_mu    = args.mu
                elif args.incl_rad is not None:
                    map_theta = args.incl_rad
                    map_mu    = float(np.cos(map_theta))
                else:
                    # Fall back to --incl index
                    map_theta = float(angular_grid["inclination_angles"][args.incl])
                    map_mu    = float(np.cos(map_theta))

                if args.azim_rad is not None:
                    map_chi = args.azim_rad
                else:
                    map_chi = float(angular_grid["azimuthal_angles"][args.azim])
            else:
                parser.error("--map requires angular-grid metadata in the HDF5 file.")

            freq_idx = args.freq_idx
            print(
                f"[map] reading spatial map: "
                f"mu={map_mu:.4f}, chi={map_chi:.4f} rad, freq_idx={freq_idx}"
            )

            result = read_map_data(file_name, mu=map_mu, chi=map_chi,
                                   freq_idx=freq_idx)

            # Frequency label for the title
            freqs       = output_freq_grid["frequencies"]
            wl_angstrom = freq_to_angstrom(freqs)
            wl_air      = vacuum_to_air(wl_angstrom)
            wl_label    = f"{wl_air[freq_idx]:.4f} Å  (freq_idx={freq_idx})"

            print(
                f"[map] matched direction: "
                f"theta={result['incl_rad']:.4f} rad, "
                f"mu={result['mu']:.4f}, "
                f"chi={result['azim_rad']:.4f} rad, "
                f"angular distance={np.degrees(result['angular_distance_rad']):.4f}°"
            )

            stokes_labels = [
                ("Stokes I",       result["stokes_I"]),
                ("Stokes Q/I (%)", result["stokes_QI"]),
                ("Stokes U/I (%)", result["stokes_UI"]),
                ("Stokes V/I (%)", result["stokes_VI"]),
            ]

            fig, axs = plt.subplots(2, 2, figsize=(10, 8))
            ax = axs.ravel()

            for panel, (label, data) in enumerate(stokes_labels):
                # Use a diverging colormap centred on zero for Q/U/V,
                # and a sequential map for I.
                if panel == 0:
                    cmap = "inferno"
                    norm = None
                else:
                    vmax = np.nanmax(np.abs(data))
                    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
                    cmap = "RdBu_r"

                im = ax[panel].imshow(
                    data.T,           # transpose so x is horizontal, y vertical
                    origin="lower",
                    aspect="equal",
                    cmap=cmap,
                    norm=norm,
                )
                ax[panel].set_title(label, fontsize=10)
                ax[panel].set_xlabel("x index")
                ax[panel].set_ylabel("y index")
                fig.colorbar(im, ax=ax[panel], fraction=0.046, pad=0.04)

            abs_path   = os.path.abspath(file_name)
            path_parts = abs_path.split(os.sep)[-3:]
            path_title = "/".join(path_parts)

            title_text = fig.suptitle(
                f"{path_title}\nSpatial map  "
                rf"$\theta$={result['incl_rad']:.4f},  "
                rf"$\mu$={result['mu']:.4f},  "
                rf"$\chi$={result['azim_rad']:.4f}  |  {wl_label}",
                fontsize=12, fontweight="bold", color="darkblue")
            title_text.set_bbox(
                dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))
            fig.tight_layout()
            fig.subplots_adjust(top=0.85)

            if args.save:
                fig.savefig(args.save)
                print(f"Saved map to {args.save}")
            else:
                plt.show()