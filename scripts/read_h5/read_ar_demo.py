import h5py
import numpy as np
import argparse
import time
import itertools
from pathlib import Path
import matplotlib.pyplot as plt

####################################################################


def speed_of_light():
    return 2.99792458e10

####################################################################


def cm_to_angstrom(cm):
    return cm * 1e8

####################################################################


def angstrom_to_cm(cm):
    return cm * 1e-8

####################################################################


def freq_to_angstrom(freq):
    return cm_to_angstrom(speed_of_light() / freq)


####################################################################
def vacuum_to_air(wave, to_air_limit=200.0):
    """
    https://github.com/ITA-Solar/rh/blob/master/idl/vacuumtoair.pro
    """
    wave2 = wave * wave

    fact = 1.0 + 2.735182e-4 + (1.314182e0 + 2.76249e+4 / wave2) / wave2
    fact = fact * (wave > to_air_limit) + 1.0 * (wave < to_air_limit)

    return wave / fact


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


def THDF_read_ar_field_directions_grid_from_hdf5(h5file, i, j, mu, chi):
    """
    Read the /emergent_beams group and return a dict:
      {
        'I': np.ndarray(dtype=float64),
        'QI_pc': np.ndarray(dtype=float64),
        'UI_pc': np.ndarray(dtype=float64),
        'VQ_pc': np.ndarray(dtype=float64),
      }
    """

    close_when_done = False
    if isinstance(h5file, str):
        f = h5py.File(h5file, "r")
        close_when_done = True
    else:
        f = h5file

    try:
        grp = f[field_names["emergent_beams_group"]]
    except KeyError:
        if close_when_done:
            f.close()
        raise KeyError("Group '/emergent_beams' not found")

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


def plot_field_data_multi_direction(directions_field_data, frequencies, i, j):
    freqs = np.asarray(frequencies['frequencies'], dtype=np.float64)
    wavelengths_vacuum = freq_to_angstrom(freqs)
    wavelengths_air = vacuum_to_air(wavelengths_vacuum)
    sort_idx = np.argsort(wavelengths_air)
    wavelengths_air = wavelengths_air[sort_idx]

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True)
    axs = axes.flatten()

    cmap = plt.get_cmap('tab10')
    for k, direction_data in enumerate(directions_field_data):
        color = cmap(k % 10)
        direction_label = (
            f"dir {direction_data['dir_index']} "
            f"(mu={direction_data['mu']:.3f}, chi={direction_data['chi']:.3f})"
        )

        I = np.asarray(direction_data['field_data']['I'])[sort_idx]
        Q = np.asarray(direction_data['field_data']['QI_pc'])[sort_idx]
        U = np.asarray(direction_data['field_data']['UI_pc'])[sort_idx]
        V = np.asarray(direction_data['field_data']['VI_pc'])[sort_idx]

        axs[0].plot(wavelengths_air, I, '-o', color=color, label=direction_label)
        axs[1].plot(wavelengths_air, Q, '-o', color=color)
        axs[2].plot(wavelengths_air, U, '-o', color=color)
        axs[3].plot(wavelengths_air, V, '-o', color=color)

    # axs[0].set_title('I')
    axs[0].set_ylabel('I')
    axs[0].legend(title='Directions', fontsize='small')

    # axs[1].set_title('Q/I %')
    axs[1].set_ylabel('Q/I %')

    # axs[2].set_title('U/I %')
    axs[2].set_xlabel('Wavelength [Angstrom, air]')
    axs[2].set_ylabel('U/I %')

    # axs[3].set_title('V/I %')
    axs[3].set_xlabel('Wavelength [Angstrom, air]')
    axs[3].set_ylabel('V/I %')

    fig.suptitle(f'Field at (i={i}, j={j})')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()


def read_field_data(h5file, i, j, mu, chi):
    start = time.perf_counter()
    field_data = THDF_read_ar_field_directions_grid_from_hdf5(
        h5file,
        i,
        j,
        mu,
        chi,
    )
    elapsed = time.perf_counter() - start
    return field_data, elapsed


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


def plot_compare_field_data(h5_data, csv_data, i, j, dir_index, mu, chi, csv_path):
    components = [
        ("I", "I"),
        ("Q/I %", "QI_pc"),
        ("U/I %", "UI_pc"),
        ("V/I %", "VI_pc"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    axs = axes.flatten()

    for ax, (title, key) in zip(axs, components):
        h5_values = np.asarray(h5_data[key], dtype=np.float64)
        csv_values = np.asarray(csv_data[key], dtype=np.float64)

        if h5_values.shape != csv_values.shape:
            raise ValueError(
                f"Shape mismatch for {key}: HDF5 {h5_values.shape} vs CSV {csv_values.shape}"
            )

        sample_idx = np.arange(h5_values.size)
        ax.plot(sample_idx, h5_values, "-", label="HDF5", marker="o")
        ax.plot(sample_idx, csv_values, "-", label="CSV", marker="x")
        ax.set_title(title)
        ax.set_xlabel("Sample")
        ax.set_ylabel(title)
        ax.grid(True, alpha=0.3)

    axs[0].legend()

    fig.suptitle(
        f"Compare dir={dir_index}, i={i}, j={j}, mu={mu:.3f}, chi={chi:.3f}\n{csv_path.name}"
    )
    plt.show()


def compare_field_data(input_path, i, j, dir_index):
    h5_path = resolve_h5_input_path(input_path)
    directions = THDF_read_ar_directions_from_hdf5(str(h5_path))

    if dir_index < 0 or dir_index >= directions["N_directions"]:
        raise IndexError(
            f"dir_index={dir_index} out of range [0, {directions['N_directions'] - 1}]"
        )

    mu = directions["mu"][dir_index]
    chi = directions["chi"][dir_index]
    h5_data, elapsed = read_field_data(str(h5_path), i, j, mu, chi)

    csv_path = find_compare_csv_path(h5_path.parent, mu, chi, i, j)
    csv_data = read_compare_csv_data(csv_path)

    print_directions_table(directions, selected_indices=[dir_index])
    print(
        f"Read HDF5 direction {dir_index} (mu={mu:.3f}, chi={chi:.3f}) "
        f"at (i={i}, j={j}) in {elapsed:.6f} seconds"
    )
    print(f"Loaded comparison CSV: {csv_path}")

    plot_compare_field_data(h5_data, csv_data, i, j,
                            dir_index, mu, chi, csv_path)


def read_map_data(h5file, mu, chi, freq_idx):
    close_when_done = False
    if isinstance(h5file, str):
        f = h5py.File(h5file, "r")
        close_when_done = True
    else:
        f = h5file

    try:
        grp = f[field_names["emergent_beams_group"]]
    except KeyError:
        if close_when_done:
            f.close()
        raise KeyError("Group '/emergent_beams' not found")

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
    if freq_idx < 0 or freq_idx >= n_frequencies:
        if close_when_done:
            f.close()
        raise IndexError(
            f"freq_idx={freq_idx} out of range [0, {n_frequencies - 1}]"
        )

    map_data = {
        "I": np.array(ds_I[:, :, freq_idx], dtype=np.float64),
        "QI_pc": np.array(ds_QI[:, :, freq_idx], dtype=np.float64),
        "UI_pc": np.array(ds_UI[:, :, freq_idx], dtype=np.float64),
        "VI_pc": np.array(ds_VI[:, :, freq_idx], dtype=np.float64),
    }

    if close_when_done:
        f.close()

    return map_data


def plot_map_data(map_data, mu, chi, freq_value=None, freq_idx=None):
    components = [
        ("I", map_data["I"]),
        ("Q/I %", map_data["QI_pc"]),
        ("U/I %", map_data["UI_pc"]),
        ("V/I %", map_data["VI_pc"]),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    axs = axes.flatten()

    for ax, (title, values) in zip(axs, components):
        image = ax.imshow(values, origin="lower", aspect="auto")
        ax.set_title(title)
        ax.set_xlabel("j")
        ax.set_ylabel("i")
        fig.colorbar(image, ax=ax)

    title_parts = [f"mu={mu:.3f}", f"chi={chi:.3f}"]
    if freq_idx is not None:
        title_parts.append(f"freq_idx={freq_idx}")
    if freq_value is not None:
        title_parts.append(f"freq={freq_value:.6e}")
    fig.suptitle("Map slices at " + ", ".join(title_parts))
    plt.show()


def print_directions_table(directions, selected_indices=None):
    n = directions["N_directions"]
    chi = np.asarray(directions["chi"], dtype=np.float64)
    mu = np.asarray(directions["mu"], dtype=np.float64)
    selected = set(selected_indices or [])

    print("Beam Directions:")
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


def main():
    parser = argparse.ArgumentParser(
        description="Read AR beam directions and frequencies from HDF5 file."
    )
    parser.add_argument("h5file", type=str,
                        help="Path to the HDF5 file or to a directory containing it.")
    parser.add_argument("-i", "--icoord", type=int, default=3,
                        help="Spatial grid index i (default: 3).")
    parser.add_argument("-j", "--jcoord", type=int, default=3,
                        help="Spatial grid index j (default: 3).")
    parser.add_argument("-d", "--dir-index", type=int, nargs="+", default=[4],
                        dest="dir_indices",
                        help="One or more direction indices for mu and chi (default: 4).")
    parser.add_argument("-m", "--map", type=int, nargs=2, metavar=("DIR_INDEX", "FREQ_INDEX"),
                        help="Read and plot 2D maps for one direction index and one frequency index.")
    parser.add_argument("--compare", type=int, nargs=3, metavar=("I", "J", "DIR_INDEX"),
                        help="Compare HDF5 data with CSV for coordinates i, j and one direction index.")
    args = parser.parse_args()

    if args.compare is not None:
        compare_i, compare_j, compare_dir_index = args.compare
        compare_field_data(args.h5file, compare_i,
                           compare_j, compare_dir_index)
        raise SystemExit(0)

    h5_path = resolve_h5_input_path(args.h5file)

    directions = THDF_read_ar_directions_from_hdf5(str(h5_path))
    frequencies = THDF_read_ar_frequencies_from_hdf5(str(h5_path))

    if args.map is not None:
        map_dir_index, map_freq_index = args.map

        if map_dir_index < 0 or map_dir_index >= directions['N_directions']:
            parser.error(f"--map direction index {map_dir_index} out of range "
                         f"[0, {directions['N_directions'] - 1}].")

        if map_freq_index < 0 or map_freq_index >= frequencies['N_frequencies']:
            parser.error(f"--map frequency index {map_freq_index} out of range "
                         f"[0, {frequencies['N_frequencies'] - 1}].")

        print_directions_table(directions, selected_indices=[map_dir_index])

        mu = directions['mu'][map_dir_index]
        chi = directions['chi'][map_dir_index]
        map_data = read_map_data(str(h5_path), mu, chi, map_freq_index)
        freq_value = frequencies['frequencies'][map_freq_index]

        print(f"Read map for direction {map_dir_index} "
              f"(mu={mu:.3f}, chi={chi:.3f}) at frequency index {map_freq_index} "
              f"(freq={freq_value:.6e})")

        plot_map_data(
            map_data,
            mu,
            chi,
            freq_value=freq_value,
            freq_idx=map_freq_index,
        )
        raise SystemExit(0)

    for dir_index in args.dir_indices:
        if dir_index < 0 or dir_index >= directions['N_directions']:
            parser.error(f"--dir-index {dir_index} out of range "
                         f"[0, {directions['N_directions'] - 1}].")

    print_directions_table(directions, selected_indices=args.dir_indices)

    # print("\nFrequencies:")
    # print(f"Number of frequencies: {frequencies['N_frequencies']}")
    # print(f"Frequencies: {frequencies['frequencies']}")

    isv = args.icoord
    jsv = args.jcoord

    directions_field_data = []
    for dir_index in args.dir_indices:
        mu = directions['mu'][dir_index]
        chi = directions['chi'][dir_index]
        field_data, elapsed = read_field_data(str(h5_path), isv, jsv, mu, chi)
        print(
            f"Read direction {dir_index} (mu={mu:.3f}, chi={chi:.3f}) in {elapsed:.6f} seconds")
        directions_field_data.append({
            "dir_index": dir_index,
            "mu": mu,
            "chi": chi,
            "field_data": field_data,
        })

    plot_field_data_multi_direction(
        directions_field_data, frequencies, isv, jsv)


if __name__ == "__main__":
    main()
