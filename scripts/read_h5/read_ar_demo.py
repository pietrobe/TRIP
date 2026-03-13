import h5py
import numpy as np
import argparse
import time
import itertools
import matplotlib.pyplot as plt


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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Read AR beam directions and frequencies from HDF5 file."
    )
    parser.add_argument("h5file", type=str,
                        help="Path to the HDF5 file.")
    parser.add_argument("-i", "--icoord", type=int, default=3,
                        help="Spatial grid index i (default: 3).")
    parser.add_argument("-j", "--jcoord", type=int, default=3,
                        help="Spatial grid index j (default: 3).")
    parser.add_argument("-d", "--dir-index", type=int, default=4,
                        dest="dir_index",
                        help="Index into the directions arrays for mu and chi (default: 4).")
    args = parser.parse_args()

    directions = THDF_read_ar_directions_from_hdf5(args.h5file)
    frequencies = THDF_read_ar_frequencies_from_hdf5(args.h5file)

    print("Beam Directions:")
    print(f"Number of directions: {directions['N_directions']}")
    print(f"Chi angles: {directions['chi']}")
    print(f"Mu angles: {directions['mu']}")

    print("\nFrequencies:")
    print(f"Number of frequencies: {frequencies['N_frequencies']}")
    print(f"Frequencies: {frequencies['frequencies']}")

    isv = args.icoord
    jsv = args.jcoord
    if args.dir_index < 0 or args.dir_index >= directions['N_directions']:
        parser.error(f"--dir-index {args.dir_index} out of range "
                     f"[0, {directions['N_directions'] - 1}].")
    mu = directions['mu'][args.dir_index]
    chi = directions['chi'][args.dir_index]

    start = time.perf_counter()
    
    data_field_list = []
    
    for (i, j) in itertools.product(range(5), range(5)):
    
        field_data_tmp = THDF_read_ar_field_directions_grid_from_hdf5(args.h5file,
                                                              i, j,
                                                              mu, chi)
        data_field_list.append(field_data_tmp)
        
        
        if (i == isv) and (j == jsv):
            field_data = field_data_tmp
            
    elapsed = time.perf_counter() - start
    print(f"Elapsed: {elapsed:.6f} seconds")
    
    print(field_data)

    print(f"\nField data at (i={isv}, j={jsv}), mu={mu}, chi={chi}:")
    print(f"I: {field_data['I']}")
    print(f"QI_pc: {field_data['QI_pc']}")
    print(f"UI_pc: {field_data['UI_pc']}")
    print(f"VI_pc: {field_data['VI_pc']}")
    
    print(f"Elapsed: {elapsed:.6f} seconds")
    # Plot I, QI_pc, UI_pc, VI_pc vs frequencies in 2x2 tiled panels
    freqs = frequencies['frequencies']
    I = np.asarray(field_data['I'])
    Q = np.asarray(field_data['QI_pc'])
    U = np.asarray(field_data['UI_pc'])
    V = np.asarray(field_data['VI_pc'])

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True)
    axs = axes.flatten()

    axs[0].plot(freqs, I, '-o', color='C0')
    axs[0].set_title('I')
    axs[0].set_ylabel('I')

    axs[1].plot(freqs, Q, '-o', color='C1')
    axs[1].set_title('QI_pc')
    axs[1].set_ylabel('QI_pc')

    axs[2].plot(freqs, U, '-o', color='C2')
    axs[2].set_title('UI_pc')
    axs[2].set_xlabel('Frequency')
    axs[2].set_ylabel('UI_pc')

    axs[3].plot(freqs, V, '-o', color='C3')
    axs[3].set_title('VI_pc')
    axs[3].set_xlabel('Frequency')
    axs[3].set_ylabel('VI_pc')

    fig.suptitle(f'Field at (i={isv}, j={jsv}), mu={mu:.3f}, chi={chi:.3f}')
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()