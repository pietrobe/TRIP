import h5py
import numpy as np
import argparse
import os


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
        N = int(grp["N_frequencies"][()])
    else:
        N = int(grp["frequencies"].shape[0])

    frequencies = np.array(grp["frequencies"], dtype=np.float64)

    if close_when_done:
        f.close()

    return {"N_frequencies": N, "frequencies": frequencies}


def THDF_read_KQ_table_from_hdf5(h5file):
    """
    Read the KQ table from an HDF5 file.

    `h5file` may be a path (str) or an open h5py.File.
    """
    close_when_done = False
    if isinstance(h5file, str):
        f = h5py.File(h5file, "r")
        close_when_done = True
    else:
        f = h5file

    try:
        grp = f["/KQ_table"]
    except KeyError:
        try:
            grp = f["/KQ_matrix"]
        except KeyError:
            if close_when_done:
                f.close()
            raise KeyError("Groups '/KQ_table' and '/KQ_matrix' not found")

    def get_dataset(name):
        if name in grp:
            data = np.array(grp[name], dtype=np.int32)
        else:
            raise KeyError("Dataset '/KQ_table/%s' not found" % name)
        return data

    KQ_table = get_dataset("KQ_table")
    KQ_compressed_real = get_dataset("KQ_compressed_real")
    KQ_compressed_imag = get_dataset("KQ_compressed_imag")

    if close_when_done:
        f.close()

    return (KQ_table, KQ_compressed_real, KQ_compressed_imag)


def THDF_read_JKQ_CRD_from_hdf5(h5file, i, j, k):
    """
    Read JKQ CRD data from an HDF5 file at spatial index (i, j, k).

    CRD data does not use normalization multipliers or frequency dimension.
    The datasets are stored as:
      /JKQ/JKQ_re[N_x][N_y][N_z][N_JKQ_real]
      /JKQ/JKQ_im[N_x][N_y][N_z][N_JKQ_imag]

    `h5file` may be a path (str) or an open h5py.File.
    """
    close_when_done = False
    if isinstance(h5file, str):
        f = h5py.File(h5file, "r")
        close_when_done = True
    else:
        f = h5file

    try:
        grp = f["/JKQ"]
    except KeyError:
        if close_when_done:
            f.close()
        raise KeyError("Group '/JKQ' not found")

    def get_JKQ(name):
        if name in grp:
            data = np.array(grp[name][i, j, k, :])
        else:
            raise KeyError("Dataset '/JKQ/%s' not found" % name)
        return data

    JKQ_re = get_JKQ("JKQ_re")
    JKQ_im = get_JKQ("JKQ_im")

    if close_when_done:
        f.close()

    return (JKQ_re, JKQ_im)


def THDF_read_JKQ_CRD_full(h5file):
    """
    Read the full JKQ CRD datasets from an HDF5 file.

    Returns the full 4D arrays for real and imaginary parts.
    """
    close_when_done = False
    if isinstance(h5file, str):
        f = h5py.File(h5file, "r")
        close_when_done = True
    else:
        f = h5file

    try:
        grp = f["/JKQ"]
    except KeyError:
        if close_when_done:
            f.close()
        raise KeyError("Group '/JKQ' not found")

    JKQ_re = np.array(grp["JKQ_re"])
    JKQ_im = np.array(grp["JKQ_im"])

    if close_when_done:
        f.close()

    return (JKQ_re, JKQ_im)


def THDF_assemble_JKQ_CRD(KQ_table, KQ_compressed_real, KQ_compressed_imag,
                          JKQ_re, JKQ_im):
    """
    Assemble JKQ CRD data into (K, Q) keyed dictionaries.

    CRD data has no frequency dimension. JKQ_re and JKQ_im are 1D arrays
    with one value per KQ component, indexed by the full KQ_table (not compressed).

    For negative Q, symmetry relations are applied:
      Real part:  J(K, -Q) = (-1)^Q * J(K, Q)
      Imag part:  J(K, -Q) = (-1)^(Q+1) * J(K, Q)
      For Q=0:    Imag part is zero.
    """
    KQ_real = dict()
    KQ_imag = dict()

    for idx, KQ_itr in enumerate(KQ_table):
        K = int(KQ_itr[0])
        Q = int(KQ_itr[1])

        val_re = float(JKQ_re[idx])
        val_im = float(JKQ_im[idx])

        KQ_real[(K, Q)] = val_re
        KQ_imag[(K, Q)] = val_im

        if Q < 0:
            KQ_real[(K, -Q)] = (-1)**Q * val_re
            KQ_imag[(K, -Q)] = (-1)**(Q+1) * val_im

    for idx, KQ_itr in enumerate(KQ_table):
        K = int(KQ_itr[0])
        Q = int(KQ_itr[1])
        if Q == 0:
            KQ_imag[(K, Q)] = 0.0

    return (KQ_real, KQ_imag)


def plot_KQ_CRD(KQ_real, KQ_imag, frequencies=None, file_name=None, i=0, j=0, k=0):
    import matplotlib.pyplot as plt

    keys = list(KQ_real.keys())
    n_plots = len(keys)

    fig, axes = plt.subplots(3, 3, figsize=(15, 12))
    axes = axes.flatten()

    for idx, key in enumerate(keys):
        if idx >= 9:
            break
        K = key[0]
        Q = key[1]
        ax = axes[idx]
        real_val = KQ_real[key]
        imag_val = KQ_imag[key]
        if isinstance(real_val, np.ndarray) and real_val.ndim > 0:
            if frequencies is not None:
                ax.plot(frequencies, real_val, label='Real Part', color='blue')
                ax.plot(frequencies, imag_val, label='Imaginary Part', color='red')
                ax.set_xlabel('Frequency (Hz)')
            else:
                ax.plot(real_val, label='Real Part', color='blue')
                ax.plot(imag_val, label='Imaginary Part', color='red')
                ax.set_xlabel('Index')
        else:
            ax.bar(0, real_val, color='blue', alpha=0.7)
            ax.bar(1, imag_val, color='red', alpha=0.7)
            ax.set_xticks([0, 1])
            ax.set_xticklabels(['Real', 'Imag'])

            va_re = 'bottom' if real_val >= 0 else 'top'
            va_im = 'bottom' if imag_val >= 0 else 'top'
            ax.text(0, real_val, f"{real_val:.3e}",
                    ha='center', va=va_re, fontsize=8, fontweight='bold')
            ax.text(1, imag_val, f"{imag_val:.3e}",
                    ha='center', va=va_im, fontsize=8, fontweight='bold')

        ax.set_title(f"K={K}, Q={Q}")
        ax.set_ylabel('Value')
        if isinstance(real_val, np.ndarray) and real_val.ndim > 0:
            ax.legend()
        ax.grid()

    for idx in range(n_plots, 9):
        axes[idx].axis('off')

    abs_path = os.path.abspath(file_name)
    path_parts = abs_path.split(os.sep)[-3:]
    path_title = "/".join(path_parts)

    title_text = fig.suptitle(f"{path_title}\nJKQ CRD Real and Imaginary Parts (i={i}, j={j}, k={k})",
                              fontsize=16, fontweight='bold', color='darkblue')
    title_text.set_bbox(dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Read JKQ CRD data from HDF5 file")
    parser.add_argument("--file", "-f", type=str, default="JKQ_CRD_field.h5",
                        help="HDF5 file name (default: JKQ_CRD_field.h5)")

    parser.add_argument("--i", "-i", type=int, default=2,
                        help="Index i for JKQ CRD (default: 2)")
    parser.add_argument("--j", "-j", type=int, default=2,
                        help="Index j for JKQ CRD (default: 2)")
    parser.add_argument("--k", "-k", type=int, default=0,
                        help="Index k for JKQ CRD (default: 0)")

    parser.add_argument("--plot", "-p", action="store_true",
                        help="Show matplotlib plot")

    args = parser.parse_args()

    file_name = args.file
    i = args.i
    j = args.j
    k = args.k

    (KQ, KQ_compressed_real, KQ_compressed_imag) = THDF_read_KQ_table_from_hdf5(file_name)
    (JKQ_re, JKQ_im) = THDF_read_JKQ_CRD_from_hdf5(file_name, i, j, k)
    (KQ_real, KQ_imag) = THDF_assemble_JKQ_CRD(KQ,
                                               KQ_compressed_real,
                                               KQ_compressed_imag,
                                               JKQ_re,
                                               JKQ_im)

    abs_path = os.path.abspath(file_name)
    path_parts = abs_path.split(os.sep)[-3:]
    path_title = "/".join(path_parts)

    print(f"File:    {path_title}")
    print(f"Grid:    N_x=8, N_y=8, N_z=143, N_frequencies=96")
    print(f"Point:   i={i}, j={j}, k={k}")
    print()
    print(f"{'K':>3}  {'Q':>3}  {'Re(J)':>18}  {'Im(J)':>18}")
    print(f"{'-'*3:>3}  {'-'*3:>3}  {'-'*18}  {'-'*18}")

    for key in sorted(KQ_real.keys(), key=lambda x: (x[0], x[1])):
        K = key[0]
        Q = key[1]
        print(f"{K:>3}  {Q:>3}  {KQ_real[key]:>18.6e}  {KQ_imag[key]:>18.6e}")

    if args.plot:
        plot_KQ_CRD(KQ_real, KQ_imag, file_name=file_name, i=i, j=j, k=k)
