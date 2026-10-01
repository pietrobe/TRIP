import h5py
import numpy as np
import argparse
import time
import itertools
from pathlib import Path
from dataclasses import dataclass
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize, TwoSlopeNorm
from scipy.fft import dct, idct
from wavelength_utils import freq_to_angstrom, vacuum_to_air


new_naming = True


beams_directions_names = {
    "beams_directions_group": "/beams_directions",
    "N_directions": "N_directions",
    "chi": "chi",
    "mu": "mu",
}

frequencies_grid_names = {
    "frequencies_grid_group": "/frequencies_grid",
    "N_frequencies": "N_frequencies",
    "frequencies": "frequencies",
}

field_names = {
    "emergent_beams_group": "/emergent_beams",
    "beams_I": "beams_I",
    "beams_QI_pc": "beams_QI_pc",
    "beams_UI_pc": "beams_UI_pc",
    "beams_VI_pc": "beams_VI_pc",
}

# telescope_degradation.py writes the noisy beams either in
# 'emergent_beams_noisy' or in 'emergent_beams_noisy_{noise_name}'.
noisy_beams_prefix = "emergent_beams_noisy"

# telescope_degradation.py --add_slit writes the slit beams either in
# 'emergent_beams_slit' or in 'emergent_beams_slit_{variant}'.
slit_beams_prefix = "emergent_beams_slit"


def resolve_beams_group(h5file, show_noise):
    """
    Return the name of the beams group to read from.

    With show_noise, look for a noisy group ('emergent_beams_noisy', or
    'emergent_beams_noisy_{noise_name}') at the root of the file; if none
    is present, warn and fall back to '/emergent_beams'.
    """
    default_group = field_names["emergent_beams_group"]
    if not show_noise:
        return default_group

    with h5py.File(str(h5file), "r") as f:
        candidates = sorted(
            name for name in f if name.startswith(noisy_beams_prefix)
        )

    if not candidates:
        print(f"Warning: no '{noisy_beams_prefix}' group in '{h5file}', "
              f"using '{default_group}'")
        return default_group

    if noisy_beams_prefix in candidates:
        group_name = noisy_beams_prefix
    else:
        group_name = candidates[0]
        if len(candidates) > 1:
            print(f"Warning: several noisy groups in '{h5file}' "
                  f"({', '.join(candidates)}), using '{group_name}'")

    print(f"Reading noisy beams from '/{group_name}' in '{h5file}'")
    return f"/{group_name}"


def resolve_slit_group(h5file, slit_name):
    """
    Return the name of the slit beams group to read from.

    `slit_name` is the optional argument of --plot_slit: an empty value
    (the bare flag) selects 'emergent_beams_slit', any other value selects
    'emergent_beams_slit_{slit_name}'. If that group is not in the file,
    list the slit groups that are present and exit.
    """
    group_name = slit_beams_prefix
    if slit_name:
        group_name = f"{slit_beams_prefix}_{slit_name}"

    with h5py.File(str(h5file), "r") as f:
        if group_name not in f:
            available = sorted(name for name in f
                               if name.startswith(slit_beams_prefix))
            print(f"Error: no '{group_name}' group in '{h5file}'")
            if available:
                print(f"Available slit groups: {', '.join(available)}")
            else:
                print("No slit group in this file: run "
                      "telescope_degradation.py --add_slit first.")
            raise SystemExit(1)

    print(f"Reading slit beams from '/{group_name}' in '{h5file}'")
    return f"/{group_name}"


def add_noise_label(error_label, show_noise):
    """Mark the plot titles when the noisy beams are shown."""
    if not show_noise:
        return error_label
    return f"{error_label}, noisy" if error_label else "noisy"


def add_slit_label(error_label, slit_group):
    """Mark the plot titles when a slit beams group is shown."""
    if slit_group is None:
        return error_label
    name = slit_group.lstrip("/")
    return f"{error_label}, {name}" if error_label else name


def get_plot_wavelengths(frequencies, no_air=False):
    """Return wavelengths in the medium requested for plotting."""
    freqs = np.asarray(frequencies["frequencies"], dtype=np.float64)
    wavelengths_vacuum = np.asarray(
        freq_to_angstrom(freqs), dtype=np.float64
    )
    if no_air:
        return wavelengths_vacuum
    return np.asarray(vacuum_to_air(wavelengths_vacuum), dtype=np.float64)


def THDF_read_ar_directions_from_hdf5(h5file):
    """
    Read the /beams_directions group and return a dict:
        {
            'N_directions': int,
            'chi' : np.ndarray(dtype=float64),
            'mu' : np.ndarray(dtype=float64),
        }
        the arrays chi and mu are of length N_directions.
    """

    close_when_done = False
    if isinstance(h5file, str):
        f = h5py.File(h5file, "r")
        close_when_done = True
    else:
        f = h5file
    try:
        grp = f[beams_directions_names["beams_directions_group"]]
    except KeyError:
        if close_when_done:
            f.close()
        raise KeyError("Group '/beams_directions' not found")
    N_val = grp[beams_directions_names["N_directions"]][()]

    N = int(N_val.item() if hasattr(N_val, 'item')
            else N_val[0] if hasattr(N_val, '__getitem__') else N_val)

    chi = np.array(grp[beams_directions_names["chi"]], dtype=np.float64)
    mu = np.array(grp[beams_directions_names["mu"]], dtype=np.float64)

    if close_when_done:
        f.close()

    return {
        "N_directions": N,
        "chi": chi,
        "mu": mu,
    }


def THDF_read_ar_frequencies_from_hdf5(h5file):
    """
    Read the /frequencies_grid group and return a dict:
        {
            'frequencies': np.ndarray(dtype=float64),
            'N_frequencies': int,
        }
    """

    close_when_done = False
    if isinstance(h5file, str):
        f = h5py.File(h5file, "r")
        close_when_done = True
    else:
        f = h5file

    try:
        grp = f[frequencies_grid_names["frequencies_grid_group"]]
    except KeyError:
        if close_when_done:
            f.close()
        raise KeyError("Group '/frequencies_grid' not found")

    N_val = grp[frequencies_grid_names["N_frequencies"]][()]
    N = int(N_val.item() if hasattr(N_val, 'item')
            else N_val[0] if hasattr(N_val, '__getitem__') else N_val)

    frequencies = np.array(
        grp[frequencies_grid_names["frequencies"]], dtype=np.float64)
    if close_when_done:
        f.close()
    return {
        "N_frequencies": N,
        "frequencies": frequencies,
    }


def THDF_read_ar_field_directions_grid_from_hdf5(h5file, i, j, mu, chi, beams_group=None):
    """
    Read the /emergent_beams group (or `beams_group`, e.g. the noisy one)
    and return a dict:
      {
        'I': np.ndarray(dtype=float64),
        'QI_pc': np.ndarray(dtype=float64),
        'UI_pc': np.ndarray(dtype=float64),
        'VQ_pc': np.ndarray(dtype=float64),
      }
    """

    beams_group = beams_group or field_names["emergent_beams_group"]

    close_when_done = False
    if isinstance(h5file, str):
        f = h5py.File(h5file, "r")
        close_when_done = True
    else:
        f = h5file

    try:
        grp = f[beams_group]
    except KeyError:
        if close_when_done:
            f.close()
        raise KeyError(f"Group '{beams_group}' not found")

    # make subgroup name beam_I, beam_QI_pc, etc.
    subgroup_I = f"{field_names['beams_I']}"
    subgroup_QI_pc = f"{field_names['beams_QI_pc']}"
    subgroup_UI_pc = f"{field_names['beams_UI_pc']}"
    subgroup_VI_pc = f"{field_names['beams_VI_pc']}"

    dataset_name = f"mu_{mu:1.3f}_chi_{chi:1.3f}"

    I = np.array(grp[subgroup_I][dataset_name][i, j], dtype=np.float64)
    QI_pc = np.array(grp[subgroup_QI_pc][dataset_name][i, j], dtype=np.float64)
    UI_pc = np.array(grp[subgroup_UI_pc][dataset_name][i, j], dtype=np.float64)
    VQ_pc = np.array(grp[subgroup_VI_pc][dataset_name][i, j], dtype=np.float64)

    return {
        "I": I,
        "QI_pc": QI_pc,
        "UI_pc": UI_pc,
        "VI_pc": VQ_pc,
    }


def print_beams_dataset_info(h5file, mu, chi, beams_group=None, label=None):
    """
    Print shape, dtype and group attributes of the beam datasets that are
    about to be read and return the (N_i, N_j, N_frequencies) shape of the
    beams_I dataset. Exit if the group or the datasets are missing.

    The slit groups resample the spatial axes, so their shape differs
    from the one of '/emergent_beams': this is what tells the user which
    (i, j) coordinates are actually available.
    """
    beams_group = beams_group or field_names["emergent_beams_group"]
    dataset_name = f"mu_{mu:1.3f}_chi_{chi:1.3f}"

    close_when_done = False
    if isinstance(h5file, (str, Path)):
        f = h5py.File(str(h5file), "r")
        close_when_done = True
    else:
        f = h5file

    try:
        try:
            grp = f[beams_group]
        except KeyError:
            print(f"Error: group '{beams_group}' not found in '{h5file}'")
            raise SystemExit(1)

        header = f"\nDataset info: '{beams_group}/{dataset_name}' in '{h5file}'"
        if label:
            header += f" [{label}]"
        print(header)

        shape = None
        for key in ("beams_I", "beams_QI_pc", "beams_UI_pc", "beams_VI_pc"):
            name = field_names[key]
            if name not in grp or dataset_name not in grp[name]:
                print(f"Error: dataset '{name}/{dataset_name}' not found in "
                      f"'{beams_group}' of '{h5file}'")
                raise SystemExit(1)
            ds = grp[name][dataset_name]
            print(f"  {name}/{dataset_name}: shape={ds.shape}, dtype={ds.dtype}")
            if shape is None:
                shape = tuple(ds.shape)

        for attr_name, attr_value in grp.attrs.items():
            print(f"  attr {attr_name} = {attr_value}")
    finally:
        if close_when_done:
            f.close()

    print(f"  valid coordinates: i in [0, {shape[0] - 1}], "
          f"j in [0, {shape[1] - 1}], N_frequencies = {shape[2]}")
    return shape


def check_beams_coordinates(shape, i, j, n=None, beams_group=None):
    """
    Exit with an error message when the requested (i, j) point, or the
    n x n region starting at it, does not fit the dataset spatial size.
    """
    region = 1 if n is None else int(n)
    N_i, N_j = shape[0], shape[1]

    if i >= 0 and j >= 0 and i + region <= N_i and j + region <= N_j:
        return

    where = f"(i={i}, j={j})"
    if n is not None:
        where = f"region ({i}:{i + region}, {j}:{j + region})"
    group_label = f" of '{beams_group}'" if beams_group else ""

    if region > N_i or region > N_j:
        limits = (f"a {region} x {region} region does not fit at all")
    else:
        limits = (f"i must be in [0, {N_i - region}] and "
                  f"j in [0, {N_j - region}]")
    print(f"Error: {where} out of bounds{group_label}: the spatial size is "
          f"({N_i}, {N_j}), so {limits}.")
    raise SystemExit(1)


def plot_field_data_multi_direction(directions_field_data, frequencies, i, j,
                                    marker=None, error_label=None,
                                    spectral_pixels=None, mean_region=None,
                                    no_air=False):
    wavelengths = get_plot_wavelengths(frequencies, no_air=no_air)
    sort_idx = np.argsort(wavelengths)
    wavelengths = wavelengths[sort_idx]
    wavelength_medium = "vacuum" if no_air else "air"

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True)
    axs = axes.flatten()
    boundaries = get_spectral_boundaries(
        wavelengths, spectral_pixels, label="profile plot"
    )

    cmap = plt.get_cmap('tab10')
    line_fmt = '-' if marker is None else f'-{marker}'
    original_label_added = False
    for k, direction_data in enumerate(directions_field_data):
        color = cmap(k % 10)
        direction_label = direction_data.get('label')
        if direction_label is None:
            direction_label = (
                f"dir {direction_data['dir_index']} "
                f"(mu={direction_data['mu']:.3f}, chi={direction_data['chi']:.3f})"
            )

        I = np.asarray(direction_data['field_data']['I'])[sort_idx]
        Q = np.asarray(direction_data['field_data']['QI_pc'])[sort_idx]
        U = np.asarray(direction_data['field_data']['UI_pc'])[sort_idx]
        V = np.asarray(direction_data['field_data']['VI_pc'])[sort_idx]

        axs[0].plot(wavelengths, I, line_fmt,
                    color=color, label=direction_label)
        axs[1].plot(wavelengths, Q, line_fmt, color=color)
        axs[2].plot(wavelengths, U, line_fmt, color=color)
        axs[3].plot(wavelengths, V, line_fmt, color=color)

        original_data = direction_data.get('original_field_data')
        if original_data is not None:
            I0 = np.asarray(original_data['I'])[sort_idx]
            Q0 = np.asarray(original_data['QI_pc'])[sort_idx]
            U0 = np.asarray(original_data['UI_pc'])[sort_idx]
            V0 = np.asarray(original_data['VI_pc'])[sort_idx]
            orig_label = 'original profile' if not original_label_added else None
            axs[0].plot(wavelengths, I0, '-', color='red', linewidth=1.2,
                        label=orig_label)
            axs[1].plot(wavelengths, Q0, '-', color='red', linewidth=1.2)
            axs[2].plot(wavelengths, U0, '-', color='red', linewidth=1.2)
            axs[3].plot(wavelengths, V0, '-', color='red', linewidth=1.2)
            original_label_added = True

    # axs[0].set_title('I')
    axs[0].set_ylabel('I')
    axs[0].legend(title='Directions', fontsize='small')

    # axs[1].set_title('Q/I %')
    axs[1].set_ylabel('Q/I %')

    # axs[2].set_title('U/I %')
    axs[2].set_xlabel(f'Wavelength [Angstrom, {wavelength_medium}]')
    axs[2].set_ylabel('U/I %')

    # axs[3].set_title('V/I %')
    axs[3].set_xlabel(f'Wavelength [Angstrom, {wavelength_medium}]')
    axs[3].set_ylabel('V/I %')

    title = f'Field at (i={i}, j={j})'
    if mean_region is not None:
        mi, mj, mn = mean_region
        title += f' mean region ({mi}:{mi+mn}, {mj}:{mj+mn})'
    if error_label:
        title += f' [{error_label}]'
    fig.suptitle(title)
    draw_spectral_boundaries(axs, boundaries)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()


def read_field_data(h5file, i, j, mu, chi, beams_group=None):
    shape = print_beams_dataset_info(h5file, mu, chi, beams_group=beams_group)
    check_beams_coordinates(shape, i, j, beams_group=beams_group)

    start = time.perf_counter()
    field_data = THDF_read_ar_field_directions_grid_from_hdf5(
        h5file,
        i,
        j,
        mu,
        chi,
        beams_group=beams_group,
    )
    elapsed = time.perf_counter() - start
    return field_data, elapsed


def read_field_data_mean(h5file, i, j, n, mu, chi, beams_group=None):
    beams_group = beams_group or field_names["emergent_beams_group"]

    shape = print_beams_dataset_info(h5file, mu, chi, beams_group=beams_group)
    check_beams_coordinates(shape, i, j, n=n, beams_group=beams_group)

    start = time.perf_counter()

    close_when_done = False
    if isinstance(h5file, str):
        f = h5py.File(h5file, "r")
        close_when_done = True
    else:
        f = h5file

    try:
        grp = f[beams_group]
    except KeyError:
        if close_when_done:
            f.close()
        raise KeyError(f"Group '{beams_group}' not found")

    dataset_name = f"mu_{mu:1.3f}_chi_{chi:1.3f}"

    try:
        ds_I = grp[field_names["beams_I"]][dataset_name]
        ds_QI = grp[field_names["beams_QI_pc"]][dataset_name]
        ds_UI = grp[field_names["beams_UI_pc"]][dataset_name]
        ds_VI = grp[field_names["beams_VI_pc"]][dataset_name]
    except KeyError:
        if close_when_done:
            f.close()
        raise KeyError(
            f"Dataset '{dataset_name}' not found for one or more beam components"
        )

    N_i, N_j = ds_I.shape[0], ds_I.shape[1]
    if i < 0 or j < 0 or i + n > N_i or j + n > N_j:
        if close_when_done:
            f.close()
        raise IndexError(
            f"Region ({i}:{i+n}, {j}:{j+n}) out of bounds for dataset "
            f"shape ({N_i}, {N_j}, ...)"
        )

    I = np.mean(np.array(ds_I[i:i+n, j:j+n], dtype=np.float64), axis=(0, 1))
    QI_pc = np.mean(np.array(ds_QI[i:i+n, j:j+n], dtype=np.float64), axis=(0, 1))
    UI_pc = np.mean(np.array(ds_UI[i:i+n, j:j+n], dtype=np.float64), axis=(0, 1))
    VI_pc = np.mean(np.array(ds_VI[i:i+n, j:j+n], dtype=np.float64), axis=(0, 1))

    if close_when_done:
        f.close()

    elapsed = time.perf_counter() - start
    return {
        "I": I,
        "QI_pc": QI_pc,
        "UI_pc": UI_pc,
        "VI_pc": VI_pc,
    }, elapsed


def resolve_h5_input_path(input_path):
    path = Path(input_path)
    if path.is_dir():
        h5_path = path / "emergent_field_abitrary_Omega.h5"
        if not h5_path.is_file():
            raise FileNotFoundError(
                f"Could not find '{h5_path.name}' in directory '{path}'"
            )
        return h5_path

    if not path.is_file():
        raise FileNotFoundError(f"Input path '{path}' does not exist")

    return path


def format_direction_value_for_csv(value):
    return f"{value:0.4f}".replace(".", "")


def find_compare_csv_path(directory, mu, chi, i, j):
    mu_str = format_direction_value_for_csv(mu)
    chi_str = format_direction_value_for_csv(chi)
    expected_name = (
        f"emergent_field_abitrary_Omega.h5_mu{mu_str}_chi{chi_str}_{i}_{j}.csv"
    )
    expected_path = directory / expected_name

    if expected_path.is_file():
        return expected_path

    fallback_pattern = f"emergent_field_abitrary_Omega.h5_mu*_chi*_{i}_{j}.csv"
    matches = sorted(directory.glob(fallback_pattern))
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise FileNotFoundError(
            "Multiple comparison CSV files matched coordinates "
            f"(i={i}, j={j}) in '{directory}'. Expected '{expected_name}'."
        )

    raise FileNotFoundError(
        f"Comparison CSV not found. Expected '{expected_name}' in '{directory}'."
    )


def read_compare_csv_data(csv_path):
    data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
    data = np.atleast_2d(data)

    if data.shape[1] != 4:
        raise ValueError(
            f"Expected 4 columns in '{csv_path}', found {data.shape[1]}"
        )

    return {
        "I": np.asarray(data[:, 0], dtype=np.float64),
        "QI_pc": np.asarray(data[:, 1], dtype=np.float64),
        "UI_pc": np.asarray(data[:, 2], dtype=np.float64),
        "VI_pc": np.asarray(data[:, 3], dtype=np.float64),
    }


def print_wavelength_sampling_info(frequencies, label="plot", lamda_range=None,
                                   no_air=False):
    """Print bandwidth and local nodes-per-Angstrom stats for the plot grid."""
    wavelengths = np.sort(get_plot_wavelengths(frequencies, no_air=no_air))

    if lamda_range is not None:
        lo, delta = lamda_range
        hi = lo + delta
        mask = (wavelengths >= lo) & (wavelengths <= hi)
        wavelengths = wavelengths[mask]

    if wavelengths.size < 2:
        print(f"[{label}] wavelength stats unavailable: need at least 2 nodes")
        return

    dlam = np.diff(wavelengths)
    dlam = dlam[dlam > 0]
    if dlam.size == 0:
        print(f"[{label}] wavelength stats unavailable: non-increasing wavelength grid")
        return

    bandwidth_angstrom = wavelengths[-1] - wavelengths[0]
    nodes_per_angstrom = 1.0 / dlam
    min_nodes_per_angstrom = float(np.min(nodes_per_angstrom))
    max_nodes_per_angstrom = float(np.max(nodes_per_angstrom))

    print(
        f"[{label}] badwidth/bandwidth = {bandwidth_angstrom:.6f} A, "
        f"min nodes/A = {min_nodes_per_angstrom:.6f}, "
        f"max nodes/A = {max_nodes_per_angstrom:.6f}"
    )


def get_spectral_boundaries(wavelengths_air, spectral_pixels, label="plot"):
    """
    Return wavelength boundaries centered at +/- spectral_pixels/2 in index space.
    For uneven grids, boundaries are interpolated in wavelength.
    """
    if spectral_pixels is None:
        return None

    wl = np.asarray(wavelengths_air, dtype=np.float64)
    n = wl.size
    pixels = int(spectral_pixels)

    if pixels <= 0:
        print(f"[{label}] ignoring --spectral_pixels={pixels}: must be > 0")
        return None
    if pixels >= n:
        print(
            f"[{label}] ignoring --spectral_pixels={pixels}: "
            f"grid has {n} points"
        )
        return None

    idx = np.arange(n, dtype=np.float64)
    center = 0.5 * (n - 1)
    half = 0.5 * pixels
    left_idx = center - half
    right_idx = center + half

    left_wl = float(np.interp(left_idx, idx, wl))
    right_wl = float(np.interp(right_idx, idx, wl))
    if left_wl > right_wl:
        left_wl, right_wl = right_wl, left_wl

    print(
        f"[{label}] spectral boundaries for N={pixels}: "
        f"{left_wl:.6f} A to {right_wl:.6f} A"
    )
    return left_wl, right_wl


def draw_spectral_boundaries(axs, boundaries):
    if boundaries is None:
        return
    left_wl, right_wl = boundaries
    for ax in axs:
        xmin, xmax = ax.get_xlim()
        # Highlight the regions outside the selected spectral-pixel window.
        if left_wl > xmin:
            ax.axvspan(xmin, left_wl, color="red", alpha=0.12, zorder=0)
        if right_wl < xmax:
            ax.axvspan(right_wl, xmax, color="red", alpha=0.12, zorder=0)
        ax.axvline(left_wl, color="black", linestyle="--", linewidth=1.1)
        ax.axvline(right_wl, color="black", linestyle="--", linewidth=1.1)


def plot_compare_field_data(h5_data, csv_data, i, j, dir_index, mu, chi, csv_path, marker=None, error_label=None):
    components = [
        ("I", "I"),
        ("Q/I %", "QI_pc"),
        ("U/I %", "UI_pc"),
        ("V/I %", "VI_pc"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    axs = axes.flatten()

    h5_fmt = '-' if marker is None else f'-{marker}'
    csv_fmt = '-' if marker is None else f'-{marker}'

    for ax, (title, key) in zip(axs, components):
        h5_values = np.asarray(h5_data[key], dtype=np.float64)
        csv_values = np.asarray(csv_data[key], dtype=np.float64)

        if h5_values.shape != csv_values.shape:
            raise ValueError(
                f"Shape mismatch for {key}: HDF5 {h5_values.shape} vs CSV {csv_values.shape}"
            )

        sample_idx = np.arange(h5_values.size)
        ax.plot(sample_idx, h5_values, h5_fmt, label="HDF5")
        ax.plot(sample_idx, csv_values, csv_fmt, label="CSV")
        ax.set_title(title)
        ax.set_xlabel("Sample")
        ax.set_ylabel(title)
        ax.grid(True, alpha=0.3)

    axs[0].legend()

    title = f"Compare dir={dir_index}, i={i}, j={j}, mu={mu:.3f}, chi={chi:.3f}\n{csv_path.name}"
    if error_label:
        title += f' [{error_label}]'
    fig.suptitle(title)
    plt.show()


def compare_field_data(input_path, i, j, dir_index, marker=None, diff_path=None, relerr=False, adjrelerr=False,
                       show_noise=False, slit_name=None, mean=None,
                       no_air=False):
    h5_path = resolve_h5_input_path(input_path)
    directions = THDF_read_ar_directions_from_hdf5(str(h5_path))
    frequencies = THDF_read_ar_frequencies_from_hdf5(str(h5_path))

    if dir_index < 0 or dir_index >= directions["N_directions"]:
        raise IndexError(
            f"dir_index={dir_index} out of range [0, {directions['N_directions'] - 1}]"
        )

    mu = directions["mu"][dir_index]
    chi = directions["chi"][dir_index]
    if slit_name is not None:
        beams_group = resolve_slit_group(h5_path, slit_name)
    else:
        beams_group = resolve_beams_group(h5_path, show_noise)
    if mean is not None:
        h5_data, elapsed = read_field_data_mean(
            str(h5_path), i, j, mean[2], mu, chi,
            beams_group=beams_group)
    else:
        h5_data, elapsed = read_field_data(str(h5_path), i, j, mu, chi,
                                           beams_group=beams_group)
    if diff_path is not None:
        if slit_name is not None:
            diff_beams_group = resolve_slit_group(diff_path, slit_name)
        else:
            diff_beams_group = resolve_beams_group(diff_path, show_noise)
        if mean is not None:
            h5_data2, elapsed2 = read_field_data_mean(
                str(diff_path), i, j, mean[2], mu, chi,
                beams_group=diff_beams_group)
        else:
            h5_data2, elapsed2 = read_field_data(
                str(diff_path), i, j, mu, chi,
                beams_group=diff_beams_group)
        if adjrelerr:
            print_adj_relative_error(h5_data, h5_data2,
                                     f"dir={dir_index}, i={i}, j={j}")
            h5_data = adj_relative_error_field_data(h5_data, h5_data2)
        elif relerr:
            print_relative_error(h5_data, h5_data2,
                                 f"dir={dir_index}, i={i}, j={j}")
            h5_data = relative_error_field_data(h5_data, h5_data2)
        else:
            h5_data = subtract_field_data(h5_data, h5_data2)
        elapsed = elapsed + elapsed2

    csv_path = find_compare_csv_path(h5_path.parent, mu, chi, i, j)
    csv_data = read_compare_csv_data(csv_path)

    print_directions_table(directions, selected_indices=[dir_index])
    print(
        f"Read HDF5 direction {dir_index} (mu={mu:.3f}, chi={chi:.3f}) "
        f"at (i={i}, j={j}) in {elapsed:.6f} seconds"
    )
    print(f"Loaded comparison CSV: {csv_path}")

    if adjrelerr:
        error_label = "adj. rel. err."
    elif relerr:
        error_label = "rel. err."
    else:
        error_label = None
    error_label = add_noise_label(error_label, show_noise)
    if slit_name is not None:
        error_label = add_slit_label(error_label, beams_group)
    print_wavelength_sampling_info(
        frequencies, label="compare plot", no_air=no_air
    )
    plot_compare_field_data(h5_data, csv_data, i, j,
                            dir_index, mu, chi, csv_path, marker=marker,
                            error_label=error_label)


def read_map_data(h5file, mu, chi, freq_idx, beams_group=None):
    beams_group = beams_group or field_names["emergent_beams_group"]

    print_beams_dataset_info(h5file, mu, chi, beams_group=beams_group)

    close_when_done = False
    if isinstance(h5file, str):
        f = h5py.File(h5file, "r")
        close_when_done = True
    else:
        f = h5file

    try:
        grp = f[beams_group]
    except KeyError:
        if close_when_done:
            f.close()
        raise KeyError(f"Group '{beams_group}' not found")

    dataset_name = f"mu_{mu:1.3f}_chi_{chi:1.3f}"

    try:
        ds_I = grp[field_names["beams_I"]][dataset_name]
        ds_QI = grp[field_names["beams_QI_pc"]][dataset_name]
        ds_UI = grp[field_names["beams_UI_pc"]][dataset_name]
        ds_VI = grp[field_names["beams_VI_pc"]][dataset_name]
    except KeyError:
        if close_when_done:
            f.close()
        raise KeyError(
            f"Dataset '{dataset_name}' not found for one or more beam components"
        )

    n_frequencies = ds_I.shape[2]
    if freq_idx is not None and (freq_idx < 0 or freq_idx >= n_frequencies):
        if close_when_done:
            f.close()
        raise IndexError(
            f"freq_idx={freq_idx} out of range [0, {n_frequencies - 1}]"
        )

    if freq_idx is None:
        freq_idx = slice(None)
    else:
        freq_idx = int(freq_idx)

    map_data = {
        "I": np.array(ds_I[:, :, freq_idx], dtype=np.float64),
        "QI_pc": np.array(ds_QI[:, :, freq_idx], dtype=np.float64),
        "UI_pc": np.array(ds_UI[:, :, freq_idx], dtype=np.float64),
        "VI_pc": np.array(ds_VI[:, :, freq_idx], dtype=np.float64),
    }

    if close_when_done:
        f.close()

    return map_data


def read_slit_data(h5file, mu, chi, slit_axis, slit_index, beams_group=None):
    beams_group = beams_group or field_names["emergent_beams_group"]

    shape = print_beams_dataset_info(h5file, mu, chi, beams_group=beams_group)
    n_axis = shape[0] if slit_axis == "i" else shape[1]
    if slit_index < 0 or slit_index >= n_axis:
        print(f"Error: --slit {slit_axis} index {slit_index} out of bounds: "
              f"the spatial size is ({shape[0]}, {shape[1]}), so {slit_axis} "
              f"must be in [0, {n_axis - 1}].")
        raise SystemExit(1)

    close_when_done = False
    if isinstance(h5file, str):
        f = h5py.File(h5file, "r")
        close_when_done = True
    else:
        f = h5file

    try:
        grp = f[beams_group]
    except KeyError:
        if close_when_done:
            f.close()
        raise KeyError(f"Group '{beams_group}' not found")

    dataset_name = f"mu_{mu:1.3f}_chi_{chi:1.3f}"

    try:
        ds_I = grp[field_names["beams_I"]][dataset_name]
        ds_QI = grp[field_names["beams_QI_pc"]][dataset_name]
        ds_UI = grp[field_names["beams_UI_pc"]][dataset_name]
        ds_VI = grp[field_names["beams_VI_pc"]][dataset_name]
    except KeyError:
        if close_when_done:
            f.close()
        raise KeyError(
            f"Dataset '{dataset_name}' not found for one or more beam components"
        )

    if slit_axis == "i":
        slit_data = {
            "I": np.array(ds_I[slit_index, :, :], dtype=np.float64),
            "QI_pc": np.array(ds_QI[slit_index, :, :], dtype=np.float64),
            "UI_pc": np.array(ds_UI[slit_index, :, :], dtype=np.float64),
            "VI_pc": np.array(ds_VI[slit_index, :, :], dtype=np.float64),
        }
        spatial_label = "j"
    else:
        slit_data = {
            "I": np.array(ds_I[:, slit_index, :], dtype=np.float64),
            "QI_pc": np.array(ds_QI[:, slit_index, :], dtype=np.float64),
            "UI_pc": np.array(ds_UI[:, slit_index, :], dtype=np.float64),
            "VI_pc": np.array(ds_VI[:, slit_index, :], dtype=np.float64),
        }
        spatial_label = "i"

    if close_when_done:
        f.close()

    return slit_data, spatial_label


def plot_map_data(map_data, mu, chi, freq_value=None, freq_idx=None, log_scale=False, cmap=False, error_label=None, diff_plot=False):
    components = [
        ("I", map_data["I"], "gray"),
        ("Q/I %", map_data["QI_pc"], "RdBu_r"),
        ("U/I %", map_data["UI_pc"], "RdBu_r"),
        ("V/I %", map_data["VI_pc"], "RdBu_r"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    axs = axes.flatten()

    for ax, (title, values, cmap_name) in zip(axs, components):
        if title == "I":
            if diff_plot:
                cmap_name = "RdBu_r"
                vmax = np.nanmax(np.abs(values))
                if vmax > 0:
                    norm = TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)
                else:
                    norm = None
            else:
                norm = LogNorm() if log_scale else None
        else:
            vmax = np.nanmax(np.abs(values))
            if vmax > 0:
                norm = TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax)
            else:
                norm = None
        image = ax.imshow(values, origin="lower", aspect="auto", cmap=cmap_name, norm=norm)
        ax.set_title(title)
        ax.set_xlabel("j")
        ax.set_ylabel("i")
        if cmap:
            ax.set_box_aspect(mu)
        fig.colorbar(image, ax=ax)

    title_parts = [f"mu={mu:.3f}", f"chi={chi:.3f}"]
    if freq_idx is not None:
        title_parts.append(f"freq_idx={freq_idx}")
    if freq_value is not None:
        title_parts.append(f"freq={freq_value:.6e}")
    title = "Map slices at " + ", ".join(title_parts)
    if error_label:
        title += f' [{error_label}]'
    fig.suptitle(title)
    plt.show()


def plot_slit_data(slit_data, spatial_label, slit_index, mu, chi, frequencies,
                   log_scale=False, lamda_range=None, cslit=False, scslit=None,
                   csv_path=None, error_label=None, spectral_pixels=None,
                   no_air=False):
    components = [
        ("I", slit_data["I"], "gray"),
        ("Q/I %", slit_data["QI_pc"], "RdBu_r"),
        ("U/I %", slit_data["UI_pc"], "RdBu_r"),
        ("V/I %", slit_data["VI_pc"], "RdBu_r"),
    ]

    wavelengths = get_plot_wavelengths(frequencies, no_air=no_air)
    sort_idx = np.argsort(wavelengths)
    wavelengths = wavelengths[sort_idx]
    wavelength_medium = "vacuum" if no_air else "air"

    if lamda_range is not None:
        lo, delta = lamda_range
        hi = lo + delta
        mask = (wavelengths >= lo) & (wavelengths <= hi)
        keep = np.where(mask)[0]
        if keep.size == 0:
            raise ValueError(
                f"No wavelengths in range [{lo}, {hi}] Å")
        wavelengths = wavelengths[keep]
        sort_idx = sort_idx[keep] if keep.size > 0 else sort_idx

    n_freq = len(wavelengths)
    wl_edges = np.empty(n_freq + 1)
    wl_edges[0] = wavelengths[0] - (wavelengths[1] - wavelengths[0]) / 2
    wl_edges[1:-1] = (wavelengths[1:] + wavelengths[:-1]) / 2
    wl_edges[-1] = wavelengths[-1] + (wavelengths[-1] - wavelengths[-2]) / 2

    key_map = {"I": "I", "Q/I %": "QI_pc", "U/I %": "UI_pc", "V/I %": "VI_pc"}
    sorted_arrays = {k: slit_data[v][:, sort_idx] for k, v in key_map.items()}

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    axs = axes.flatten()
    boundaries = get_spectral_boundaries(
        wavelengths, spectral_pixels, label="slit plot"
    )

    for ax, (title, _values, cmap_name) in zip(axs, components):
        sorted_values = sorted_arrays[title]
        if title == "I":
            if log_scale:
                pos = sorted_values[sorted_values > 0]
                vmin = max(np.nanmin(pos) if pos.size > 0 else 1e-300, 1e-300)
                vmax = np.nanmax(sorted_values)
                norm = LogNorm(vmin=vmin, vmax=vmax)
            else:
                norm = None
        else:
            vmax_d = np.nanmax(np.abs(sorted_values))
            if vmax_d > 0:
                if scslit is not None:
                    vmax_d = vmax_d * scslit / 100
                norm = TwoSlopeNorm(vcenter=0, vmin=-vmax_d, vmax=vmax_d)
            else:
                norm = None
        n_spatial = sorted_values.shape[0]
        y_edges = np.arange(n_spatial + 1) - 0.5
        ylabel = f"Spatial index ({spatial_label})"
        image = ax.pcolormesh(wl_edges, y_edges, sorted_values,
                              cmap=cmap_name, norm=norm)
        ax.set_title(title)
        ax.set_xlabel(rf"Wavelength [$\AA$, {wavelength_medium}]")
        ax.set_ylabel(ylabel)
        if cslit:
            ax.set_box_aspect(mu)
        fig.colorbar(image, ax=ax)

    title = f"Slit at {spatial_label}={slit_index}, mu={mu:.3f}, chi={chi:.3f}"
    if error_label:
        title += f' [{error_label}]'
    fig.suptitle(title)
    draw_spectral_boundaries(axs, boundaries)

    if csv_path is not None:
        if csv_path.endswith(("/", "\\")):
            base = csv_path
            sep = ""
        elif "/" in csv_path or "\\" in csv_path:
            base = csv_path
            sep = "_"
        else:
            base = f"/tmp/{csv_path}_"
            sep = ""
        csv_pairs = [("I", "I"), ("QI", "QI_pc"), ("UI", "UI_pc"), ("VI", "VI_pc")]
        for cname, ckey in csv_pairs:
            out = Path(f"{base}{sep}{cname}.csv")
            out.parent.mkdir(parents=True, exist_ok=True)
            arr = slit_data[ckey][:, sort_idx]
            header = f"wavelength_{wavelength_medium}_angstrom," + ",".join(
                f"{spatial_label}_{i}" for i in range(arr.shape[0])
            )
            table = np.column_stack([wavelengths, arr.T])
            np.savetxt(out, table, delimiter=",", header=header, comments="")
            print(f"Saved: {out.resolve()}")

    plt.show()


def print_directions_table(directions, selected_indices=None, title="Beam Directions:"):
    n = directions["N_directions"]
    chi = np.asarray(directions["chi"], dtype=np.float64)
    mu = np.asarray(directions["mu"], dtype=np.float64)
    selected = set(selected_indices or [])

    print(title)
    print(f"Number of directions: {n}")

    idx_w = max(5, len(str(n - 1)))
    sel_w = 3
    val_w = 12

    header = f"| {'sel':>{sel_w}} | {'idx':>{idx_w}} | {'mu':>{val_w}} | {'chi':>{val_w}} |"
    sep = f"+-{'-' * sel_w}-+-{'-' * idx_w}-+-{'-' * val_w}-+-{'-' * val_w}-+"
    print(sep)
    print(header)
    print(sep)
    bold_start = "\033[1m"
    bold_end = "\033[0m"
    for idx, (mu_val, chi_val) in enumerate(zip(mu, chi)):
        marker = "*" if idx in selected else ""
        row = f"| {marker:>{sel_w}} | {idx:>{idx_w}d} | {mu_val:>{val_w}.6f} | {chi_val:>{val_w}.6f} |"
        if idx in selected:
            row = f"{bold_start}{row}{bold_end}"
        print(row)
    print(sep)
    if selected:
        print("* selected direction index (bold row)")


def compress(input_arr, target_size):
    x = np.asarray(input_arr)
    if target_size == 0:
        return np.empty(0, dtype=x.dtype)
    if x.size == 0:
        raise ValueError("compress: empty input")

    N = x.size
    M = int(target_size)
    if M == N:
        return x.copy()

    xd = x.astype(np.float64, copy=False)
    ratio = N / M
    edges = np.arange(M + 1) * ratio
    edges = np.minimum(edges, float(N))

    cs = np.concatenate(([0.0], np.cumsum(xd)))
    idx = edges.astype(np.int64)
    frac = edges - idx
    safe = np.minimum(idx, N - 1)
    F = cs[idx] + frac * xd[safe]

    out = (F[1:] - F[:-1]) / np.diff(edges)
    return out.astype(x.dtype, copy=False)


def decompress(comp, original_size):
    c = np.asarray(comp)
    if original_size == 0:
        return np.empty(0, dtype=c.dtype)
    if c.size == 0:
        raise ValueError("decompress: empty input")

    M = c.size
    N = int(original_size)
    if M == N:
        return c.copy()
    if M == 1:
        return np.full(N, c[0], dtype=c.dtype)

    cd = c.astype(np.float64, copy=False)
    pos = (np.arange(N) + 0.5) * (M / N) - 0.5
    xp = np.arange(M)
    out = np.interp(pos, xp, cd)
    return out.astype(c.dtype, copy=False)


def compress_fft(v, M):
    v = np.asarray(v, dtype=np.float64)
    N = v.size
    if M == 0:
        return np.empty(0)
    if N == 0:
        raise ValueError("compress_fft: empty input")
    if M == N:
        return v.copy()

    X = np.fft.rfft(v)
    K = M // 2 + 1

    if M < N:
        Y = X[:K].copy()
        if M % 2 == 0:
            Y[-1] = Y[-1].real
    else:
        Y = np.zeros(K, dtype=complex)
        Y[:X.size] = X

    return np.fft.irfft(Y, M) * (M / N)


def decompress_fft(c, N):
    c = np.asarray(c, dtype=np.float64)
    M = c.size
    if N == 0:
        return np.empty(0)
    if M == 0:
        raise ValueError("decompress_fft: empty input")
    if M == N:
        return c.copy()

    X = np.fft.rfft(c)
    K = N // 2 + 1

    if N > M:
        Y = np.zeros(K, dtype=complex)
        Y[:X.size] = X
    else:
        Y = X[:K].copy()
        if N % 2 == 0:
            Y[-1] = Y[-1].real

    return np.fft.irfft(Y, N) * (N / M)


@dataclass
class DCTPayload:
    indices: np.ndarray
    values: np.ndarray
    original_size: int

    @property
    def size(self):
        return self.indices.size

    def storage_floats(self):
        return 2 * self.indices.size


def compress_dct_topk(v, K):
    v = np.asarray(v, dtype=np.float64)
    N = v.size
    if K >= N:
        coefs = dct(v, type=2, norm='ortho')
        return DCTPayload(np.arange(N), coefs, N)
    coefs = dct(v, type=2, norm='ortho')
    keep = np.argpartition(np.abs(coefs), N - K)[-K:]
    keep.sort()
    return DCTPayload(keep.astype(np.int32), coefs[keep], N)


def decompress_dct_topk(payload):
    N = payload.original_size
    full = np.zeros(N, dtype=np.float64)
    full[payload.indices] = payload.values
    return idct(full, type=2, norm='ortho')


def plot_compress_test(directions_field_data, frequencies, i, j, target_size,
                       marker=None, use_fft=False, error_label=None,
                       spectral_pixels=None, no_air=False):
    components = [
        ("I", "I"),
        ("Q/I %", "QI_pc"),
        ("U/I %", "UI_pc"),
        ("V/I %", "VI_pc"),
    ]

    wavelengths = get_plot_wavelengths(frequencies, no_air=no_air)
    sort_idx = np.argsort(wavelengths)
    wavelengths = wavelengths[sort_idx]
    wavelength_medium = "vacuum" if no_air else "air"

    _compress = compress_fft if use_fft else compress
    _decompress = decompress_fft if use_fft else decompress

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True)
    axs = axes.flatten()
    boundaries = get_spectral_boundaries(
        wavelengths, spectral_pixels, label="compress plot"
    )

    cmap = plt.get_cmap('tab10')
    line_fmt = '-' if marker is None else f'-{marker}'
    recon_marker = marker if marker is not None else '.'
    dash_fmt = f'--{recon_marker}'

    bold = "\033[1m"
    reset = "\033[0m"

    def print_table(rows, key_width=10, val_width=18):
        sep = "+" + "=" * (key_width + 2) + "+" + "=" * (val_width + 2) + "+" + "=" * (val_width + 2) + "+" + "=" * (val_width + 2) + "+"
        hdr = f"| {'':>{key_width}} | {'RMS':>{val_width}} | {'NRMS':>{val_width}} | {'max_rel%':>{val_width}} |"
        print(sep)
        print(hdr)
        dash = sep.replace("=", "-")
        print(dash)
        for key, rms, nrms, max_rel in rows:
            print(f"| {key:>{key_width}} | {rms:>{val_width}.6e} | {nrms:>{val_width}.6e} | {max_rel:>{val_width}.4f} |")
        print(sep)

    for k, direction_data in enumerate(directions_field_data):
        color = cmap(k % 10)
        direction_label = (
            f"dir {direction_data['dir_index']} "
            f"(mu={direction_data['mu']:.3f}, chi={direction_data['chi']:.3f})"
        )

        rows = []
        for ax, (title, key) in zip(axs, components):
            arr = np.asarray(direction_data['field_data'][key])[sort_idx]
            comp = _compress(arr, target_size)
            recon = _decompress(comp, arr.size)
            diff = recon - arr
            rms = np.sqrt(np.mean(diff ** 2))
            norm = np.max(np.abs(arr))
            nrms = rms / norm if norm != 0 else 0.0
            max_rel = np.max(np.abs(diff)) / norm * 100 if norm != 0 else 0.0
            rows.append((key, rms, nrms, max_rel))
            ax.plot(wavelengths, arr, line_fmt,
                    color=color, label=f"{direction_label} orig")
            ax.plot(wavelengths, recon, dash_fmt,
                    color=color, label=f"{direction_label} c{target_size}")

        print(f"\n{bold}{direction_label}{reset}")
        print_table(rows)

    axs[0].set_ylabel('I')
    axs[0].legend(title='Directions', fontsize='small')

    axs[1].set_ylabel('Q/I %')

    axs[2].set_xlabel(f'Wavelength [Angstrom, {wavelength_medium}]')
    axs[2].set_ylabel('U/I %')

    axs[3].set_xlabel(f'Wavelength [Angstrom, {wavelength_medium}]')
    axs[3].set_ylabel('V/I %')

    title = f'Compression test (target={target_size}) at (i={i}, j={j})'
    if error_label:
        title += f' [{error_label}]'
    fig.suptitle(title)
    draw_spectral_boundaries(axs, boundaries)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()


def plot_compress_test_dct(directions_field_data, frequencies, i, j, K,
                           marker=None, error_label=None, spectral_pixels=None,
                           no_air=False):
    components = [
        ("I", "I"),
        ("Q/I %", "QI_pc"),
        ("U/I %", "UI_pc"),
        ("V/I %", "VI_pc"),
    ]

    wavelengths = get_plot_wavelengths(frequencies, no_air=no_air)
    sort_idx = np.argsort(wavelengths)
    wavelengths = wavelengths[sort_idx]
    wavelength_medium = "vacuum" if no_air else "air"

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True)
    axs = axes.flatten()
    boundaries = get_spectral_boundaries(
        wavelengths, spectral_pixels, label="compress DCT plot"
    )

    cmap = plt.get_cmap('tab10')
    line_fmt = '-' if marker is None else f'-{marker}'
    recon_marker = marker if marker is not None else '.'
    dash_fmt = f'--{recon_marker}'

    bold = "\033[1m"
    reset = "\033[0m"

    def print_table(rows, key_width=10, val_width=18):
        sep = "+" + "=" * (key_width + 2) + "+" + "=" * (val_width + 2) + "+" + "=" * (val_width + 2) + "+" + "=" * (val_width + 2) + "+"
        hdr = f"| {'':>{key_width}} | {'RMS':>{val_width}} | {'NRMS':>{val_width}} | {'max_rel%':>{val_width}} |"
        print(sep)
        print(hdr)
        dash = sep.replace("=", "-")
        print(dash)
        for key, rms, nrms, max_rel in rows:
            print(f"| {key:>{key_width}} | {rms:>{val_width}.6e} | {nrms:>{val_width}.6e} | {max_rel:>{val_width}.4f} |")
        print(sep)

    for k, direction_data in enumerate(directions_field_data):
        color = cmap(k % 10)
        direction_label = (
            f"dir {direction_data['dir_index']} "
            f"(mu={direction_data['mu']:.3f}, chi={direction_data['chi']:.3f})"
        )

        rows = []
        for ax, (title, key) in zip(axs, components):
            arr = np.asarray(direction_data['field_data'][key])[sort_idx]
            payload = compress_dct_topk(arr, K)
            recon = decompress_dct_topk(payload)
            diff = recon - arr
            rms = np.sqrt(np.mean(diff ** 2))
            norm = np.max(np.abs(arr))
            nrms = rms / norm if norm != 0 else 0.0
            max_rel = np.max(np.abs(diff)) / norm * 100 if norm != 0 else 0.0
            rows.append((key, rms, nrms, max_rel))
            ax.plot(wavelengths, arr, line_fmt,
                    color=color, label=f"{direction_label} orig")
            ax.plot(wavelengths, recon, dash_fmt,
                    color=color, label=f"{direction_label} K={K}")

        print(f"\n{bold}{direction_label}{reset}")
        print_table(rows)

    axs[0].set_ylabel('I')
    axs[0].legend(title='Directions', fontsize='small')

    axs[1].set_ylabel('Q/I %')

    axs[2].set_xlabel(f'Wavelength [Angstrom, {wavelength_medium}]')
    axs[2].set_ylabel('U/I %')

    axs[3].set_xlabel(f'Wavelength [Angstrom, {wavelength_medium}]')
    axs[3].set_ylabel('V/I %')

    title = f'DCT top-K compression test (K={K}) at (i={i}, j={j})'
    if error_label:
        title += f' [{error_label}]'
    fig.suptitle(title)
    draw_spectral_boundaries(axs, boundaries)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()


def subtract_field_data(data_a, data_b):
    return {
        "I": data_a["I"] - data_b["I"],
        "QI_pc": data_a["QI_pc"] - data_b["QI_pc"],
        "UI_pc": data_a["UI_pc"] - data_b["UI_pc"],
        "VI_pc": data_a["VI_pc"] - data_b["VI_pc"],
    }


def print_relative_error(data_a, data_b, label=""):
    components = [("I", "I"), ("Q/I %", "QI_pc"), ("U/I %", "UI_pc"), ("V/I %", "VI_pc")]
    print(f"\nRelative error {label}:")
    for name, key in components:
        a = np.asarray(data_a[key], dtype=np.float64)
        b = np.asarray(data_b[key], dtype=np.float64)
        denom = np.abs(a)
        mask = denom > 0
        if np.any(mask):
            rel_err = np.abs(a[mask] - b[mask]) / denom[mask]
            mean_rel = np.mean(rel_err) * 100
            max_rel = np.max(rel_err) * 100
            print(f"  {name:>8s}: mean_rel_err = {mean_rel:.4f}%, max_rel_err = {max_rel:.4f}%")
        else:
            print(f"  {name:>8s}: all zeros, cannot compute relative error")


def print_adj_relative_error(data_a, data_b, label=""):
    print("Adjusted relative error definition: |v0 - v1| / (max(|v0|) + |v0|)")
    components = [("I", "I"), ("Q/I %", "QI_pc"), ("U/I %", "UI_pc"), ("V/I %", "VI_pc")]
    print(f"\nAdjusted relative error {label}:")
    for name, key in components:
        a = np.asarray(data_a[key], dtype=np.float64)
        b = np.asarray(data_b[key], dtype=np.float64)
        max_abs = np.max(np.abs(a))
        if max_abs == 0:
            print(f"  {name:>8s}: all zeros in reference, cannot compute adjusted relative error")
            continue
        denom = max_abs + np.abs(a)
        rel_err = np.abs(a - b) / denom
        mean_rel = np.mean(rel_err) * 100
        max_rel = np.max(rel_err) * 100
        print(f"  {name:>8s}: mean_adj_rel_err = {mean_rel:.4f}%, max_adj_rel_err = {max_rel:.4f}%")


def relative_error_field_data(data_a, data_b):
    result = {}
    for key in ("I", "QI_pc", "UI_pc", "VI_pc"):
        a = np.asarray(data_a[key], dtype=np.float64)
        b = np.asarray(data_b[key], dtype=np.float64)
        denom = np.abs(a)
        mask = denom > 0
        rel_err = np.zeros_like(a)
        rel_err[mask] = np.abs(a[mask] - b[mask]) / denom[mask]
        result[key] = rel_err
    return result


def adj_relative_error_field_data(data_a, data_b):
    result = {}
    for key in ("I", "QI_pc", "UI_pc", "VI_pc"):
        a = np.asarray(data_a[key], dtype=np.float64)
        b = np.asarray(data_b[key], dtype=np.float64)
        max_abs = np.max(np.abs(a))
        if max_abs == 0:
            result[key] = np.zeros_like(a)
        else:
            result[key] = np.abs(a - b) / (max_abs + np.abs(a))
    return result


def find_closest_direction(directions, target_mu):
    """
    Find the direction whose chi is closest to 0 and, among those, whose mu
    is closest to target_mu.

    The search is lexicographic: minimize |chi| first, then |mu - target_mu|.
    Each file may have a different mu grid, so this must be applied per file.
    """
    mu = np.asarray(directions["mu"], dtype=np.float64)
    chi = np.asarray(directions["chi"], dtype=np.float64)

    min_chi = np.min(np.abs(chi))
    candidates = np.where(np.abs(chi) == min_chi)[0]
    if candidates.size == 0:
        candidates = np.arange(directions["N_directions"])

    best_idx = int(candidates[np.argmin(np.abs(mu[candidates] - target_mu))])

    return {
        "dir_index": best_idx,
        "mu": float(mu[best_idx]),
        "chi": float(chi[best_idx]),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Read AR beam directions and frequencies from HDF5 file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    g_input = parser.add_argument_group(
        "Input", "HDF5 files to read and which beam group to take.")
    g_input.add_argument("h5file", type=str, nargs="?", default=None,
                         help="Path to the HDF5 file or to a directory containing it.")
    g_input.add_argument("--file_diff", type=str, default=None,
                         help="Path to a second HDF5 file. The difference h5file - file_diff is "
                              "computed and used for all subsequent operations.")
    g_input.add_argument("--file_comp", type=str, default=None,
                         help="Path to a second HDF5 file. The profiles from both files are "
                              "plotted overlaid (no difference is computed).")
    g_input.add_argument("--show_noise", action="store_true",
                         help="Plot the noisy beams ('emergent_beams_noisy' group, written by "
                              "telescope_degradation.py --add_noise) instead of 'emergent_beams'. "
                              "Falls back to 'emergent_beams' if the noisy group is absent.")
    g_input.add_argument("--plot_slit", type=str, nargs="?", const="",
                         default=None, metavar="NAME",
                         help="Plot the slit beams ('emergent_beams_slit' group, written by "
                              "telescope_degradation.py --add_slit) instead of 'emergent_beams'. "
                              "With an argument (e.g. --plot_slit noisy) the group "
                              "'emergent_beams_slit_{NAME}' is read instead. All the other "
                              "options apply unchanged to the slit data. Exits if the group "
                              "is absent. Mutually exclusive with --show_noise.")

    g_space = parser.add_argument_group(
        "Spatial selection", "Which point (or region) of the spatial grid to profile.")
    g_space.add_argument("-i", "--icoord", type=int, default=3,
                         help="Spatial grid index i (default: 3).")
    g_space.add_argument("-j", "--jcoord", type=int, default=3,
                         help="Spatial grid index j (default: 3).")
    g_space.add_argument("--ij", type=int, nargs=2, metavar=("I", "J"),
                         help="Set both i and j coordinates (overrides -i/-j).")
    g_space.add_argument("--mean", type=int, nargs=3, metavar=("I", "J", "N"),
                         help="Average profiles in an N x N region starting at (I, J). "
                              "Alternative to --ij; mutually exclusive with it.")

    g_dir = parser.add_argument_group(
        "Direction selection", "Which beam direction (mu, chi) to plot.")
    g_dir.add_argument("-d", "--dir-index", type=int, nargs="+", default=[4],
                       dest="dir_indices",
                       help="One or more direction indices for mu and chi (default: 4).")
    g_dir.add_argument("--mu", type=float, default=None,
                       help="Select the direction with chi closest to 0 and mu "
                            "closest to this value, instead of -d. With --file_diff "
                            "the direction is looked up independently in each file, "
                            "so their mu tables may differ; the difference of the "
                            "two selected directions is plotted.")

    g_mode = parser.add_argument_group(
        "Plot modes", "Alternatives to the default Stokes-profile plot.")
    g_mode.add_argument("-m", "--map", type=int, nargs=2, metavar=("DIR_INDEX", "FREQ_INDEX"),
                        help="Read and plot 2D maps for one direction index and one frequency index.")
    g_mode.add_argument("--map_mu", type=float, nargs=2, metavar=("MU", "FREQ_INDEX"),
                        help="Like --map but selects the direction with chi closest to 0 "
                             "and mu closest to the given value.")
    g_mode.add_argument("--slit", type=str, nargs=2, metavar=("AXIS", "INDEX"),
                        help="Plot I/Q/U/V colormaps vs frequency along axis (i or j) at given index.")
    g_mode.add_argument("--compare", type=int, nargs=3, metavar=("I", "J", "DIR_INDEX"),
                        help="Compare HDF5 data with CSV for coordinates i, j and one direction index.")

    g_compress = parser.add_argument_group(
        "Compression tests", "Reduce each Stokes profile and overlay the reconstruction.")
    g_compress.add_argument("--compress", type=int, metavar="N",
                            help="Compress/decompress test (linear): reduce each Stokes profile to N elements and plot overlay.")
    g_compress.add_argument("--compressfft", type=int, metavar="N",
                            help="Compress/decompress test (FFT): reduce each Stokes profile to N elements and plot overlay.")
    g_compress.add_argument("--compressdct", type=int, metavar="K",
                            help="Compress/decompress test (DCT top-K): keep K DCT coefficients and plot overlay.")

    g_err = parser.add_argument_group(
        "Error metrics", "Only meaningful together with --file_diff.")
    g_err.add_argument("--relerr", action="store_true",
                       help="Report relative error |v0 - v1| / v0 when --file_diff is used.")
    g_err.add_argument("--adjrelerr", action="store_true",
                       help="Report adjusted relative error |v0 - v1| / (max(|v0|) + |v0|) when --file_diff is used.")

    g_appear = parser.add_argument_group(
        "Plot appearance", "Scaling, axes and styling of the produced figures.")
    g_appear.add_argument("--marker", type=str, default=None,
                          help="Marker style for plots (default: no marker).")
    g_appear.add_argument("--log", action="store_true",
                          help="Plot I map in logarithmic scale.")
    g_appear.add_argument("--no_air", action="store_true",
                          help="Use vacuum wavelengths instead of converting to air wavelengths.")
    g_appear.add_argument("--lamdarange", type=float, nargs=2, metavar=("MIN", "DELTA"),
                          help="Wavelength range [Å, air] for --slit plot: MIN to MIN+DELTA "
                               "(use vacuum wavelengths with --no_air).")
    g_appear.add_argument("--cslit", action="store_true",
                          help="Compress y-axis proportionally to mu (direction cosine): y/x = mu/1.")
    g_appear.add_argument("--cmap", action="store_true",
                          help="Compress map y-axis proportionally to mu: y/x = mu/1.")
    g_appear.add_argument("--scslit", type=float, metavar="PC",
                          help="Saturate slit colormap at PC%% of maximum (values above are clipped).")
    g_appear.add_argument("--spectral_pixels", type=int, metavar="N", default=None,
                          help="Draw boundaries at center +/- N/2 pixels in wavelength plots. "
                               "Ignored when N <= 0 or N >= number of plotted points.")

    g_out = parser.add_argument_group("Output", "Saving data to disk.")
    g_out.add_argument("--csv", type=str, nargs="?", const="/tmp/", default=None, metavar="PREFIX",
                       help="Save slit data to 4 CSV files. If PREFIX has no '/', it's a name under /tmp/. "
                            "Examples: --csv slits -> /tmp/slits_I.csv; --csv /mypath/slits -> /mypath/slits_I.csv")

    args = parser.parse_args()

    if args.ij is not None:
        args.icoord, args.jcoord = args.ij

    if args.mean is not None and args.ij is not None:
        parser.error("--mean and --ij are mutually exclusive")

    if args.plot_slit is not None and args.show_noise:
        parser.error("--plot_slit and --show_noise are mutually exclusive: "
                     "use '--plot_slit noisy' to read "
                     "'emergent_beams_slit_noisy'")

    if args.mean is not None:
        args.icoord, args.jcoord = args.mean[0], args.mean[1]

    if args.h5file is None:
        parser.print_help()
        raise SystemExit(0)

    h5_path = resolve_h5_input_path(args.h5file)
    diff_h5_path = resolve_h5_input_path(args.file_diff) if args.file_diff else None
    comp_h5_path = resolve_h5_input_path(args.file_comp) if args.file_comp else None

    if diff_h5_path is not None:
        print(f"Computing difference: {h5_path} - {diff_h5_path}")

    if args.compare is not None:
        compare_i, compare_j, compare_dir_index = args.compare
        compare_field_data(args.h5file, compare_i,
                           compare_j, compare_dir_index, marker=args.marker,
                           diff_path=diff_h5_path, relerr=args.relerr,
                           adjrelerr=args.adjrelerr, show_noise=args.show_noise,
                           slit_name=args.plot_slit, mean=args.mean,
                           no_air=args.no_air)
        raise SystemExit(0)

    error_label = None
    if args.adjrelerr:
        error_label = "adj. rel. err."
    elif args.relerr:
        error_label = "rel. err."
    error_label = add_noise_label(error_label, args.show_noise)

    if args.plot_slit is not None:
        beams_group = resolve_slit_group(h5_path, args.plot_slit)
        diff_beams_group = (resolve_slit_group(diff_h5_path, args.plot_slit)
                            if diff_h5_path is not None else None)
        error_label = add_slit_label(error_label, beams_group)
    else:
        beams_group = resolve_beams_group(h5_path, args.show_noise)
        diff_beams_group = (resolve_beams_group(diff_h5_path, args.show_noise)
                            if diff_h5_path is not None else None)

    directions = THDF_read_ar_directions_from_hdf5(str(h5_path))
    frequencies = THDF_read_ar_frequencies_from_hdf5(str(h5_path))

    if args.map is not None or args.map_mu is not None:
        if args.map is not None:
            map_dir_index, map_freq_index = args.map
        else:
            map_mu, map_freq_index = args.map_mu
            map_freq_index = int(map_freq_index)
            map_spec = find_closest_direction(directions, map_mu)
            map_dir_index = map_spec['dir_index']
            print(f"--map_mu {map_mu}: '{h5_path}' -> dir {map_dir_index} "
                  f"(mu={map_spec['mu']:.3f}, chi={map_spec['chi']:.3f})")

        diff_directions = (THDF_read_ar_directions_from_hdf5(str(diff_h5_path))
                           if diff_h5_path is not None else None)
        diff_spec = None
        if diff_directions is not None and args.map_mu is not None:
            map_mu, _ = args.map_mu
            diff_spec = find_closest_direction(diff_directions, map_mu)
            print(f"--map_mu {map_mu}: '{diff_h5_path}' -> dir {diff_spec['dir_index']} "
                  f"(mu={diff_spec['mu']:.3f}, chi={diff_spec['chi']:.3f})")

        if map_dir_index < 0 or map_dir_index >= directions['N_directions']:
            parser.error(f"--map direction index {map_dir_index} out of range "
                         f"[0, {directions['N_directions'] - 1}].")

        if map_freq_index < 0 or map_freq_index >= frequencies['N_frequencies']:
            parser.error(f"--map frequency index {map_freq_index} out of range "
                         f"[0, {frequencies['N_frequencies'] - 1}].")

        print_directions_table(directions, selected_indices=[map_dir_index],
                               title=f"Beam Directions: {h5_path.name} (source)")
        if diff_directions is not None:
            diff_selected = ([diff_spec['dir_index']] if diff_spec is not None
                             else [map_dir_index])
            print_directions_table(diff_directions, selected_indices=diff_selected,
                                   title=f"Beam Directions: {diff_h5_path.name} (diff)")

        mu = directions['mu'][map_dir_index]
        chi = directions['chi'][map_dir_index]
        map_data = read_map_data(str(h5_path), mu, chi, map_freq_index,
                                 beams_group=beams_group)
        if diff_h5_path is not None:
            if diff_spec is not None:
                mu2 = diff_spec['mu']
                chi2 = diff_spec['chi']
            else:
                mu2, chi2 = mu, chi
            map_data2 = read_map_data(str(diff_h5_path), mu2, chi2, map_freq_index,
                                      beams_group=diff_beams_group)
            if args.adjrelerr:
                print_adj_relative_error(map_data, map_data2,
                                         f"dir={map_dir_index}")
                map_data = adj_relative_error_field_data(map_data, map_data2)
            elif args.relerr:
                print_relative_error(map_data, map_data2,
                                     f"dir={map_dir_index}")
                map_data = relative_error_field_data(map_data, map_data2)
            else:
                map_data = subtract_field_data(map_data, map_data2)
        freq_value = frequencies['frequencies'][map_freq_index]

        print(f"Read map for direction {map_dir_index} "
              f"(mu={mu:.3f}, chi={chi:.3f}) at frequency index {map_freq_index} "
              f"(freq={freq_value:.6e})")
        if diff_spec is not None:
            print(f"Diff map for direction {diff_spec['dir_index']} "
                  f"(mu={mu2:.3f}, chi={chi2:.3f})")

        print_wavelength_sampling_info(
            frequencies, label="map plot", no_air=args.no_air
        )

        plot_map_data(
            map_data,
            mu,
            chi,
            freq_value=freq_value,
            freq_idx=map_freq_index,
            log_scale=args.log,
            cmap=args.cmap,
            error_label=error_label,
            diff_plot=diff_h5_path is not None,
        )

        if args.csv is not None:
            csv_path = args.csv
            if csv_path.endswith(("/", "\\")):
                base, sep = csv_path, ""
            elif "/" in csv_path or "\\" in csv_path:
                base, sep = csv_path, "_"
            else:
                base, sep = f"/tmp/{csv_path}_", ""
            csv_pairs = [("I", "I"), ("QI", "QI_pc"), ("UI", "UI_pc"), ("VI", "VI_pc")]
            for cname, ckey in csv_pairs:
                out = Path(f"{base}{sep}{cname}.csv")
                out.parent.mkdir(parents=True, exist_ok=True)
                arr = map_data[ckey]
                np.savetxt(out, arr, delimiter=",", fmt="%.8e")
                print(f"Saved: {out.resolve()}")

        raise SystemExit(0)

    if args.slit is not None:
        slit_axis = args.slit[0]
        if slit_axis not in ("i", "j"):
            parser.error("--slit axis must be 'i' or 'j'")
        slit_index = int(args.slit[1])

        dir_index = args.dir_indices[0]
        if dir_index < 0 or dir_index >= directions['N_directions']:
            parser.error(f"--slit direction index {dir_index} out of range "
                         f"[0, {directions['N_directions'] - 1}].")

        print_directions_table(directions, selected_indices=[dir_index])

        mu = directions['mu'][dir_index]
        chi = directions['chi'][dir_index]
        slit_data, spatial_label = read_slit_data(
            str(h5_path), mu, chi, slit_axis, slit_index,
            beams_group=beams_group)
        if diff_h5_path is not None:
            slit_data2, _ = read_slit_data(
                str(diff_h5_path), mu, chi, slit_axis, slit_index,
                beams_group=diff_beams_group)
            if args.adjrelerr:
                print_adj_relative_error(slit_data, slit_data2,
                                         f"dir={dir_index}, slit={slit_axis}={slit_index}")
                slit_data = adj_relative_error_field_data(slit_data, slit_data2)
            elif args.relerr:
                print_relative_error(slit_data, slit_data2,
                                     f"dir={dir_index}, slit={slit_axis}={slit_index}")
                slit_data = relative_error_field_data(slit_data, slit_data2)
            else:
                slit_data = subtract_field_data(slit_data, slit_data2)

        print(f"Slit along '{spatial_label}' at {slit_axis}={slit_index}, "
              f"direction {dir_index} (mu={mu:.3f}, chi={chi:.3f}), "
              f"shape={slit_data['I'].shape}")

        print_wavelength_sampling_info(
            frequencies, label="slit plot", lamda_range=args.lamdarange,
            no_air=args.no_air
        )

        plot_slit_data(slit_data, spatial_label, slit_index, mu, chi,
                       frequencies, log_scale=args.log,
                       lamda_range=args.lamdarange, cslit=args.cslit,
                       scslit=args.scslit, csv_path=args.csv,
                       error_label=error_label,
                       spectral_pixels=args.spectral_pixels,
                       no_air=args.no_air)
        raise SystemExit(0)

    if args.compress is not None:
        target_size = args.compress
        for dir_index in args.dir_indices:
            if dir_index < 0 or dir_index >= directions['N_directions']:
                parser.error(f"--dir-index {dir_index} out of range "
                             f"[0, {directions['N_directions'] - 1}].")

        print_directions_table(directions, selected_indices=args.dir_indices)

        isv = args.icoord
        jsv = args.jcoord
        directions_field_data = []
        for dir_index in args.dir_indices:
            mu = directions['mu'][dir_index]
            chi = directions['chi'][dir_index]
            if args.mean is not None:
                field_data, elapsed = read_field_data_mean(
                    str(h5_path), isv, jsv, args.mean[2], mu, chi,
                    beams_group=beams_group)
            else:
                field_data, elapsed = read_field_data(str(h5_path), isv, jsv, mu, chi,
                                                      beams_group=beams_group)
            if diff_h5_path is not None:
                if args.mean is not None:
                    field_data2, elapsed2 = read_field_data_mean(
                        str(diff_h5_path), isv, jsv, args.mean[2], mu, chi,
                        beams_group=diff_beams_group)
                else:
                    field_data2, elapsed2 = read_field_data(
                        str(diff_h5_path), isv, jsv, mu, chi,
                        beams_group=diff_beams_group)
                if args.adjrelerr:
                    print_adj_relative_error(field_data, field_data2,
                                             f"dir={dir_index}, i={isv}, j={jsv}")
                    field_data = adj_relative_error_field_data(field_data, field_data2)
                elif args.relerr:
                    print_relative_error(field_data, field_data2,
                                         f"dir={dir_index}, i={isv}, j={jsv}")
                    field_data = relative_error_field_data(field_data, field_data2)
                else:
                    field_data = subtract_field_data(field_data, field_data2)
                elapsed = elapsed + elapsed2
            print(
                f"Read direction {dir_index} (mu={mu:.3f}, chi={chi:.3f}) in {elapsed:.6f} seconds")
            directions_field_data.append({
                "dir_index": dir_index,
                "mu": mu,
                "chi": chi,
                "field_data": field_data,
            })

        print_wavelength_sampling_info(
            frequencies, label="compress plot", no_air=args.no_air
        )
        plot_compress_test(directions_field_data, frequencies, isv, jsv,
                           target_size, marker=args.marker,
                           error_label=error_label,
                           spectral_pixels=args.spectral_pixels,
                           no_air=args.no_air)
        raise SystemExit(0)

    if args.compressfft is not None:
        target_size = args.compressfft
        for dir_index in args.dir_indices:
            if dir_index < 0 or dir_index >= directions['N_directions']:
                parser.error(f"--dir-index {dir_index} out of range "
                             f"[0, {directions['N_directions'] - 1}].")

        print_directions_table(directions, selected_indices=args.dir_indices)

        isv = args.icoord
        jsv = args.jcoord
        directions_field_data = []
        for dir_index in args.dir_indices:
            mu = directions['mu'][dir_index]
            chi = directions['chi'][dir_index]
            if args.mean is not None:
                field_data, elapsed = read_field_data_mean(
                    str(h5_path), isv, jsv, args.mean[2], mu, chi,
                    beams_group=beams_group)
            else:
                field_data, elapsed = read_field_data(str(h5_path), isv, jsv, mu, chi,
                                                      beams_group=beams_group)
            if diff_h5_path is not None:
                if args.mean is not None:
                    field_data2, elapsed2 = read_field_data_mean(
                        str(diff_h5_path), isv, jsv, args.mean[2], mu, chi,
                        beams_group=diff_beams_group)
                else:
                    field_data2, elapsed2 = read_field_data(
                        str(diff_h5_path), isv, jsv, mu, chi,
                        beams_group=diff_beams_group)
                if args.adjrelerr:
                    print_adj_relative_error(field_data, field_data2,
                                             f"dir={dir_index}, i={isv}, j={jsv}")
                    field_data = adj_relative_error_field_data(field_data, field_data2)
                elif args.relerr:
                    print_relative_error(field_data, field_data2,
                                         f"dir={dir_index}, i={isv}, j={jsv}")
                    field_data = relative_error_field_data(field_data, field_data2)
                else:
                    field_data = subtract_field_data(field_data, field_data2)
                elapsed = elapsed + elapsed2
            print(
                f"Read direction {dir_index} (mu={mu:.3f}, chi={chi:.3f}) in {elapsed:.6f} seconds")
            directions_field_data.append({
                "dir_index": dir_index,
                "mu": mu,
                "chi": chi,
                "field_data": field_data,
            })

        print_wavelength_sampling_info(
            frequencies, label="compress FFT plot", no_air=args.no_air
        )
        plot_compress_test(directions_field_data, frequencies, isv, jsv,
                           target_size, marker=args.marker, use_fft=True,
                           error_label=error_label,
                           spectral_pixels=args.spectral_pixels,
                           no_air=args.no_air)
        raise SystemExit(0)

    if args.compressdct is not None:
        K = args.compressdct
        for dir_index in args.dir_indices:
            if dir_index < 0 or dir_index >= directions['N_directions']:
                parser.error(f"--dir-index {dir_index} out of range "
                             f"[0, {directions['N_directions'] - 1}].")

        print_directions_table(directions, selected_indices=args.dir_indices)

        isv = args.icoord
        jsv = args.jcoord
        directions_field_data = []
        for dir_index in args.dir_indices:
            mu = directions['mu'][dir_index]
            chi = directions['chi'][dir_index]
            if args.mean is not None:
                field_data, elapsed = read_field_data_mean(
                    str(h5_path), isv, jsv, args.mean[2], mu, chi,
                    beams_group=beams_group)
            else:
                field_data, elapsed = read_field_data(str(h5_path), isv, jsv, mu, chi,
                                                      beams_group=beams_group)
            if diff_h5_path is not None:
                if args.mean is not None:
                    field_data2, elapsed2 = read_field_data_mean(
                        str(diff_h5_path), isv, jsv, args.mean[2], mu, chi,
                        beams_group=diff_beams_group)
                else:
                    field_data2, elapsed2 = read_field_data(
                        str(diff_h5_path), isv, jsv, mu, chi,
                        beams_group=diff_beams_group)
                if args.adjrelerr:
                    print_adj_relative_error(field_data, field_data2,
                                             f"dir={dir_index}, i={isv}, j={jsv}")
                    field_data = adj_relative_error_field_data(field_data, field_data2)
                elif args.relerr:
                    print_relative_error(field_data, field_data2,
                                         f"dir={dir_index}, i={isv}, j={jsv}")
                    field_data = relative_error_field_data(field_data, field_data2)
                else:
                    field_data = subtract_field_data(field_data, field_data2)
                elapsed = elapsed + elapsed2
            print(
                f"Read direction {dir_index} (mu={mu:.3f}, chi={chi:.3f}) in {elapsed:.6f} seconds")
            directions_field_data.append({
                "dir_index": dir_index,
                "mu": mu,
                "chi": chi,
                "field_data": field_data,
            })

        print_wavelength_sampling_info(
            frequencies, label="compress DCT plot", no_air=args.no_air
        )
        plot_compress_test_dct(directions_field_data, frequencies, isv, jsv,
                               K, marker=args.marker,
                               error_label=error_label,
                               spectral_pixels=args.spectral_pixels,
                               no_air=args.no_air)
        raise SystemExit(0)

    if comp_h5_path is not None:
        if args.plot_slit is not None:
            comp_beams_group = resolve_slit_group(comp_h5_path, args.plot_slit)
        else:
            comp_beams_group = resolve_beams_group(comp_h5_path, args.show_noise)
        comp_directions = THDF_read_ar_directions_from_hdf5(str(comp_h5_path))

        if args.mu is not None:
            main_spec = find_closest_direction(directions, args.mu)
            comp_spec = find_closest_direction(comp_directions, args.mu)
            print(f"--mu {args.mu}: '{h5_path}' -> "
                  f"dir {main_spec['dir_index']} (mu={main_spec['mu']:.3f}, "
                  f"chi={main_spec['chi']:.3f}); '{comp_h5_path}' -> "
                  f"dir {comp_spec['dir_index']} (mu={comp_spec['mu']:.3f}, "
                  f"chi={comp_spec['chi']:.3f})")
            selected_indices = [main_spec['dir_index']]
        else:
            main_spec = None
            comp_spec = None
            selected_indices = args.dir_indices

        for dir_index in selected_indices:
            if dir_index < 0 or dir_index >= directions['N_directions']:
                parser.error(f"--dir-index {dir_index} out of range "
                             f"[0, {directions['N_directions'] - 1}].")

        comp_selected = ([comp_spec['dir_index']] if comp_spec is not None
                         else selected_indices)
        print_directions_table(directions, selected_indices=selected_indices,
                               title=f"Beam Directions: {h5_path.name} (source)")
        print_directions_table(comp_directions, selected_indices=comp_selected,
                               title=f"Beam Directions: {comp_h5_path.name} (comp)")

        isv = args.icoord
        jsv = args.jcoord

        directions_field_data = []
        for dir_index in selected_indices:
            if main_spec is not None:
                mu = main_spec['mu']
                chi = main_spec['chi']
            else:
                mu = directions['mu'][dir_index]
                chi = directions['chi'][dir_index]
            if args.mean is not None:
                field_data, elapsed = read_field_data_mean(
                    str(h5_path), isv, jsv, args.mean[2], mu, chi,
                    beams_group=beams_group)
            else:
                field_data, elapsed = read_field_data(str(h5_path), isv, jsv, mu, chi,
                                                      beams_group=beams_group)

            if comp_spec is not None:
                mu2 = comp_spec['mu']
                chi2 = comp_spec['chi']
                comp_dir_index = comp_spec['dir_index']
            else:
                mu2, chi2 = mu, chi
                comp_dir_index = dir_index
            if args.mean is not None:
                comp_field_data, elapsed2 = read_field_data_mean(
                    str(comp_h5_path), isv, jsv, args.mean[2], mu2, chi2,
                    beams_group=comp_beams_group)
            else:
                comp_field_data, elapsed2 = read_field_data(
                    str(comp_h5_path), isv, jsv, mu2, chi2,
                    beams_group=comp_beams_group)

            print(f"Read direction {dir_index} (mu={mu:.3f}, chi={chi:.3f}) "
                  f"vs comp direction {comp_dir_index} (mu={mu2:.3f}, "
                  f"chi={chi2:.3f}) in {elapsed + elapsed2:.6f} seconds")

            directions_field_data.append({
                "dir_index": dir_index,
                "mu": mu,
                "chi": chi,
                "field_data": field_data,
                "original_field_data": None,
                "label": f"source dir {dir_index} (mu={mu:.3f}, chi={chi:.3f})",
            })
            directions_field_data.append({
                "dir_index": comp_dir_index,
                "mu": mu2,
                "chi": chi2,
                "field_data": comp_field_data,
                "original_field_data": None,
                "label": f"comp dir {comp_dir_index} "
                         f"(mu={mu2:.3f}, chi={chi2:.3f})",
            })

        print_wavelength_sampling_info(
            frequencies, label="profile plot", no_air=args.no_air
        )
        plot_field_data_multi_direction(
            directions_field_data, frequencies, isv, jsv, marker=args.marker,
            error_label=error_label,
            spectral_pixels=args.spectral_pixels,
            mean_region=tuple(args.mean) if args.mean is not None else None,
            no_air=args.no_air)
        raise SystemExit(0)

    if args.mu is not None:
        main_spec = find_closest_direction(directions, args.mu)
        diff_spec = None
        diff_directions = None
        if diff_h5_path is not None:
            diff_directions = THDF_read_ar_directions_from_hdf5(str(diff_h5_path))
            diff_spec = find_closest_direction(diff_directions, args.mu)
            print(f"--mu {args.mu}: '{h5_path}' -> dir {main_spec['dir_index']} "
                  f"(mu={main_spec['mu']:.3f}, chi={main_spec['chi']:.3f}); "
                  f"'{diff_h5_path}' -> dir {diff_spec['dir_index']} "
                  f"(mu={diff_spec['mu']:.3f}, chi={diff_spec['chi']:.3f})")
        else:
            print(f"--mu {args.mu}: '{h5_path}' -> dir {main_spec['dir_index']} "
                  f"(mu={main_spec['mu']:.3f}, chi={main_spec['chi']:.3f})")
        selected_indices = [main_spec['dir_index']]
    else:
        main_spec = None
        diff_spec = None
        diff_directions = (THDF_read_ar_directions_from_hdf5(str(diff_h5_path))
                           if diff_h5_path is not None else None)
        selected_indices = args.dir_indices

    for dir_index in selected_indices:
        if dir_index < 0 or dir_index >= directions['N_directions']:
            parser.error(f"--dir-index {dir_index} out of range "
                         f"[0, {directions['N_directions'] - 1}].")

    print_directions_table(directions, selected_indices=selected_indices,
                           title=f"Beam Directions: {h5_path.name} (source)")
    if diff_h5_path is not None:
        diff_selected = ([diff_spec['dir_index']] if diff_spec is not None
                         else selected_indices)
        print_directions_table(diff_directions, selected_indices=diff_selected,
                               title=f"Beam Directions: {diff_h5_path.name} (diff)")

    # print("\nFrequencies:")
    # print(f"Number of frequencies: {frequencies['N_frequencies']}")
    # print(f"Frequencies: {frequencies['frequencies']}")

    isv = args.icoord
    jsv = args.jcoord

    directions_field_data = []
    for dir_index in selected_indices:
        if main_spec is not None:
            mu = main_spec['mu']
            chi = main_spec['chi']
        else:
            mu = directions['mu'][dir_index]
            chi = directions['chi'][dir_index]

        if args.mean is not None:
            field_data, elapsed = read_field_data_mean(
                str(h5_path), isv, jsv, args.mean[2], mu, chi,
                beams_group=beams_group)
        else:
            field_data, elapsed = read_field_data(str(h5_path), isv, jsv, mu, chi,
                                                  beams_group=beams_group)
        original_field_data = None
        # With --show_noise, overlay the original non-noisy profile in red.
        if args.show_noise and diff_h5_path is None and not args.relerr and not args.adjrelerr:
            original_field_data, _ = read_field_data(
                str(h5_path), isv, jsv, mu, chi,
                beams_group=field_names["emergent_beams_group"]
            )

        diff_dir_index = dir_index
        if diff_h5_path is not None:
            if diff_spec is not None:
                mu2 = diff_spec['mu']
                chi2 = diff_spec['chi']
                diff_dir_index = diff_spec['dir_index']
            else:
                mu2, chi2 = mu, chi
            if args.mean is not None:
                field_data2, elapsed2 = read_field_data_mean(
                    str(diff_h5_path), isv, jsv, args.mean[2], mu2, chi2,
                    beams_group=diff_beams_group)
            else:
                field_data2, elapsed2 = read_field_data(
                    str(diff_h5_path), isv, jsv, mu2, chi2,
                    beams_group=diff_beams_group)
            if args.adjrelerr:
                print_adj_relative_error(field_data, field_data2,
                                         f"dir={dir_index}, i={isv}, j={jsv}")
                field_data = adj_relative_error_field_data(field_data, field_data2)
            elif args.relerr:
                print_relative_error(field_data, field_data2,
                                     f"dir={dir_index}, i={isv}, j={jsv}")
                field_data = relative_error_field_data(field_data, field_data2)
            else:
                field_data = subtract_field_data(field_data, field_data2)
            elapsed = elapsed + elapsed2

        if diff_spec is not None:
            print(f"Read direction {dir_index} (mu={mu:.3f}, chi={chi:.3f}) - "
                  f"diff direction {diff_dir_index} (mu={mu2:.3f}, chi={chi2:.3f}) "
                  f"in {elapsed:.6f} seconds")
        else:
            print(f"Read direction {dir_index} (mu={mu:.3f}, chi={chi:.3f}) "
                  f"in {elapsed:.6f} seconds")
        directions_field_data.append({
            "dir_index": dir_index,
            "mu": mu,
            "chi": chi,
            "field_data": field_data,
            "original_field_data": original_field_data,
        })

    print_wavelength_sampling_info(
        frequencies, label="profile plot", no_air=args.no_air
    )
    plot_field_data_multi_direction(
        directions_field_data, frequencies, isv, jsv, marker=args.marker,
        error_label=error_label,
        spectral_pixels=args.spectral_pixels,
        mean_region=tuple(args.mean) if args.mean is not None else None,
        no_air=args.no_air)


if __name__ == "__main__":
    main()
