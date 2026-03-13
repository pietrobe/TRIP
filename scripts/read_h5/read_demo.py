import h5py
import numpy as np
import argparse
import os


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
        grp = f["/emergent_field"]
    except KeyError:
        if close_when_done:
            f.close()
        raise KeyError("Group '/emergent_field' not found")

    dsI = grp["emergent_I"]
    dsQ = grp["emergent_QI_pc"]
    dsU = grp["emergent_UI_pc"]
    dsV = grp["emergent_VI_pc"]

    # infer N_frequencies if not provided
    if N_frequencies is None:
        N_frequencies = dsI.shape[4]

    # The HDF5 shape is (x, y, incl, azim, freq), so use incl_i and azimuth_i directly
    stokes_I = np.array(dsI[x_i, y_i, incl_i, azimuth_i, :],
                        dtype=np.float64).reshape((N_frequencies,))
    stokes_Q = np.array(dsQ[x_i, y_i, incl_i, azimuth_i, :],
                        dtype=np.float64).reshape((N_frequencies,))
    stokes_U = np.array(dsU[x_i, y_i, incl_i, azimuth_i, :],
                        dtype=np.float64).reshape((N_frequencies,))
    stokes_V = np.array(dsV[x_i, y_i, incl_i, azimuth_i, :],
                        dtype=np.float64).reshape((N_frequencies,))

    if close_when_done:
        f.close()

    return {
        "index_i": x_i,
        "index_j": y_i,
        "index_incl": incl_i,
        "index_azimuth": azimuth_i,
        "stokes_I": stokes_I,
        "stokes_QI": stokes_Q,
        "stokes_UI": stokes_U,
        "stokes_VI": stokes_V,
    }


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
        N = int(N_val.item() if hasattr(N_val, 'item') else N_val[0] if hasattr(N_val, '__getitem__') else N_val)
    else:
        # prefer inclination_angles if present
        if "inclination_angles" in grp:
            N = int(grp["inclination_angles"].shape[0])
        else:
            if close_when_done:
                f.close()
            raise KeyError(
                "Neither 'N_directions' nor 'inclination_angles' found in group")

    N_incl_val = grp["N_inclination_angles"][()]
    N_inclination_angles = int(N_incl_val.item() if hasattr(N_incl_val, 'item') else N_incl_val[0] if hasattr(N_incl_val, '__getitem__') else N_incl_val)
    N_azim_val = grp["N_azimuthal_angles"][()]
    N_azimuthal_angles = int(N_azim_val.item() if hasattr(N_azim_val, 'item') else N_azim_val[0] if hasattr(N_azim_val, '__getitem__') else N_azim_val)

    inclination_angles = np.array(grp["inclination_angles"], dtype=np.float64)
    azimuthal_angles = np.array(grp["azimuthal_angles"], dtype=np.float64)
    incl_indices = np.array(grp["inclinations_indices"], dtype=np.int32)
    azim_indices = np.array(grp["azimuthal_indices"], dtype=np.int32)

    if close_when_done:
        f.close()

    return {
        "N_directions": N,
        "inclination_angles": inclination_angles,
        "azimuthal_angles": azimuthal_angles,
        "inclinations_indices": incl_indices,
        "azimuthal_indices": azim_indices,
        "N_inclination_angles": N_inclination_angles,
        "N_azimuthal_angles": N_azimuthal_angles,
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

    # Read N_frequencies if present, otherwise infer from the frequencies dataset length
    if "N_frequencies" in grp:
        N = int(grp["N_frequencies"][()])
    else:
        N = int(grp["frequencies"].shape[0])

    frequencies = np.array(grp["frequencies"], dtype=np.float64)

    if close_when_done:
        f.close()

    return {"N_frequencies": N, "frequencies": frequencies}


if __name__ == "__main__":
    
    # Example usage: accept optional filename from command line
    parser = argparse.ArgumentParser(description="Read HDF5 demo outputs")
    parser.add_argument("file", nargs="?", default="output_field_mpi.h5",
                        help="HDF5 file to read (default: output_field_mpi.h5)")
    parser.add_argument("--plot", action="store_true",
                        help="Plot Stokes I,Q,U,V vs frequency")
    parser.add_argument("--save", help="If set, save plot to this filename (e.g. plot.png)")
    parser.add_argument("--x", type=int, default=1, help="x index (default 1)")
    parser.add_argument("--y", type=int, default=2, help="y index (default 2)")
    parser.add_argument("--incl", type=int, default=2, help="inclination index (default 2)")
    parser.add_argument("--azim", type=int, default=2, help="azimuth index (default 2)")
    args = parser.parse_args()
    file_name = args.file

    output_freq_grid = THDF_read_frequencies_grid_from_hdf5(file_name)
    print("Number of frequencies:", output_freq_grid["N_frequencies"])
    print("Frequencies (Hz):", output_freq_grid["frequencies"])

    angular_grid = THDF_read_angular_grid_from_hdf5(file_name)

    print("Number of directions:", angular_grid["N_directions"])
    print("Inclination angles (deg):", angular_grid["inclination_angles"])
    print("Azimuthal angles (deg):", angular_grid["azimuthal_angles"])
    print("Inclination indices:", angular_grid["inclinations_indices"])
    print("Azimuthal indices:", angular_grid["azimuthal_indices"])
    print("N_inclination_angles: ", angular_grid["N_inclination_angles"])
    print("N_azimuthal_angles: ", angular_grid["N_azimuthal_angles"])
    print("N_directions:", angular_grid["N_directions"])

    # Use the value from the angular grid
    N_azimuth = angular_grid["N_azimuthal_angles"]

    output_filed = THDF_read_field_from_hdf5(
        file_name,
        N_frequencies=output_freq_grid["N_frequencies"],
        x_i=args.x,
        y_i=args.y,
        incl_i=args.incl,
        azimuth_i=args.azim,
        N_azimuth=N_azimuth,
    )

    print("Output field at (0,0) for incl=0, azimuth=0:")
    print("Stokes I:", output_filed["stokes_I"])
    print("Stokes Q:", output_filed["stokes_QI"])
    print("Stokes U:", output_filed["stokes_UI"])
    print("Stokes V:", output_filed["stokes_VI"])
    # Plot if requested (2x2 tiled plot: I,Q,U,V)
    if args.plot:
        try:
            import matplotlib.pyplot as plt
        except Exception:
            print("matplotlib not available; cannot plot")
        else:
            freqs = output_freq_grid["frequencies"]
            I = output_filed["stokes_I"]
            Q = output_filed["stokes_QI"]
            U = output_filed["stokes_UI"]
            V = output_filed["stokes_VI"]

            fig, axs = plt.subplots(2, 2, figsize=(10, 8), sharex=True)
            ax = axs.ravel()

            ax[0].plot(freqs, I, color="C0", marker="o")
            ax[0].set_title("I")
            ax[0].set_ylabel("Stokes I")
            ax[0].grid(True)

            ax[1].plot(freqs, Q, color="C1", marker="o")
            ax[1].set_title("Q")
            ax[1].set_ylabel("Stokes Q")
            ax[1].grid(True)

            ax[2].plot(freqs, U, color="C2", marker="o")
            ax[2].set_title("U")
            ax[2].set_xlabel("Frequency (Hz)")
            ax[2].set_ylabel("Stokes U")
            ax[2].grid(True)

            ax[3].plot(freqs, V, color="C3", marker="o")
            ax[3].set_title("V")
            ax[3].set_xlabel("Frequency (Hz)")
            ax[3].set_ylabel("Stokes V")
            ax[3].grid(True)

            # Extract last 3 path components from input file
            abs_path = os.path.abspath(file_name)
            path_parts = abs_path.split(os.sep)[-3:]
            path_title = "/".join(path_parts)
            
            title_text = fig.suptitle(f"{path_title}\nStokes vs Frequency (x={args.x}, y={args.y}, incl={args.incl}, azim={args.azim})", 
                                      fontsize=16, fontweight='bold', color='darkblue')
            title_text.set_bbox(dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
            fig.tight_layout()
            fig.subplots_adjust(top=0.85)

            if args.save:
                fig.savefig(args.save)
                print(f"Saved plot to {args.save}")
            else:
                plt.show()
    
    
