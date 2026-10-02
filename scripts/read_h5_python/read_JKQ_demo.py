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

    # Read N_frequencies if present, otherwise infer from the frequencies dataset length
    if "N_frequencies" in grp:
        N = int(grp["N_frequencies"][()])
    else:
        N = int(grp["frequencies"].shape[0])

    frequencies = np.array(grp["frequencies"], dtype=np.float64)

    if close_when_done:
        f.close()

    return {"N_frequencies": N, "frequencies": frequencies}


def THDF_read_KQ_grid_from_hdf5(h5file):
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

    def get_KQ(name):
        if name in grp:
            KQ = np.array(grp[name], dtype=np.int32)
        else:
            raise KeyError("Dataset '/KQ_table/%s' not found" % name)
        return KQ

    KQ = get_KQ("KQ_table")
    KQ_compressed_real = get_KQ("KQ_compressed_real")
    KQ_compressed_imag = get_KQ("KQ_compressed_imag")

    if close_when_done:
        f.close()

    return (KQ, KQ_compressed_real, KQ_compressed_imag)


def THDF_read_JKQ_from_hdf5(h5file, i, j, k):
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
            JKQ = np.array(grp[name][i, j, k, :, :])
        else:
            raise KeyError("Dataset '/JKQ/%s' not found" % name)
        return JKQ

    def get_JKQ_norm_mult(name):
        if name in grp:
            JKQ = np.array(grp[name][i, j, k, :])
        else:
            raise KeyError("Dataset '/JKQ/%s' not found" % name)
        return JKQ

    JKQ_re = get_JKQ("JKQ_re")
    JKQ_im = get_JKQ("JKQ_im")

    JKQ_norm_re = get_JKQ_norm_mult("JKQ_norm_re")
    JKQ_norm_im = get_JKQ_norm_mult("JKQ_norm_im")

    if close_when_done:
        f.close()

    return (JKQ_re, JKQ_im, JKQ_norm_re, JKQ_norm_im)


def THDF_read_JKQ_denormalize(KQ,
                              KQ_compressed_real,
                              KQ_compressed_imag,
                              JKQ_re,
                              JKQ_im,
                              JKQ_norm_re,
                              JKQ_norm_im):
    KQ_real = dict()
    KQ_imag = dict()

    for _, KQ_itr in enumerate(KQ):
        K = int(KQ_itr[0])
        Q = int(KQ_itr[1])
        KQ_real[(K, Q)] = None
        KQ_imag[(K, Q)] = None
          
    # print(KQ_real) 
    # print(KQ_imag)

    for idx, KQ_re in enumerate(KQ_compressed_real):
        K = int(KQ_re[0])
        Q = int(KQ_re[1])
        print(KQ_re)
        
        norm_mult = JKQ_norm_re[idx]
        KQ_re_denorm = np.array(JKQ_re[idx, :], dtype=np.float64) * norm_mult
        KQ_real[(K, Q)] = KQ_re_denorm
        if Q < 0:
            KQ_real[(K, -Q)] = (-1)**Q * KQ_re_denorm

    for idx, KQ_im in enumerate(KQ_compressed_imag):
        K = int(KQ_im[0])
        Q = int(KQ_im[1])
        print(KQ_im)
        
        norm_mult = JKQ_norm_im[idx]
        KQ_im_denorm = np.array(JKQ_im[idx, :], dtype=np.float64) * norm_mult
        KQ_imag[(K, Q)] = KQ_im_denorm
        
        if Q < 0:
            KQ_imag[(K, -Q)] = (-1)**(Q+1) * KQ_im_denorm
            
    for idx, KQ_itr in enumerate(KQ):
        K = int(KQ_itr[0])
        Q = int(KQ_itr[1])
        if Q == 0:
            KQ_imag[(K, Q)] = np.zeros_like(KQ_real[(K, Q)], dtype=np.float64)

    return (KQ_real, KQ_imag)

def plot_KQ(KQ_real, KQ_imag, frequencies=None, file_name=None, i=0, j=0, k=0):
    import matplotlib.pyplot as plt

    keys = list(KQ_real.keys())
    n_plots = len(keys)
    
    fig, axes = plt.subplots(3, 3, figsize=(15, 12))
    axes = axes.flatten()
    
    for idx, key in enumerate(keys):
        if idx >= 9:  # Only plot first 9
            break
        K = key[0]
        Q = key[1]
        ax = axes[idx]
        if frequencies is not None:
            ax.plot(frequencies, KQ_real[key], label='Real Part', color='blue')
            ax.plot(frequencies, KQ_imag[key], label='Imaginary Part', color='red')
            ax.set_xlabel('Frequency (Hz)')
        else:
            ax.plot(KQ_real[key], label='Real Part')
            ax.plot(KQ_imag[key], label='Imaginary Part')
            ax.set_xlabel('Index')
        ax.set_title(f"K={K}, Q={Q}")
        ax.set_ylabel('Value')
        ax.legend()
        ax.grid()
    
    # Hide unused subplots
    for idx in range(n_plots, 9):
        axes[idx].axis('off')
    
    # Extract last 3 path components from input file
    abs_path = os.path.abspath(file_name)
    path_parts = abs_path.split(os.sep)[-3:]
    path_title = "/".join(path_parts)
    
    title_text = fig.suptitle(f"{path_title}\nKQ Real and Imaginary Parts (i={i}, j={j}, k={k})", 
                              fontsize=16, fontweight='bold', color='darkblue')
    title_text.set_bbox(dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Read JKQ data from HDF5 file")
    parser.add_argument("--file", "-f", type=str, default="output_JKQ_mpi.h5",
                        help="HDF5 file name (default: output_JKQ_mpi.h5)")
    
    parser.add_argument("--i", "-i", type=int, default=2,
                        help="Index i for JKQ (default: 2)")
    parser.add_argument("--j", "-j", type=int, default=2,
                        help="Index j for JKQ (default: 2)")
    parser.add_argument("--k", "-k", type=int, default=0,
                        help="Index k for JKQ (default: 0)")
    
    args = parser.parse_args()
    
    file_name = args.file
    i = args.i
    j = args.j
    k = args.k

    output_freq_grid = THDF_read_frequencies_grid_from_hdf5(file_name)
    print("Number of frequencies:", output_freq_grid["N_frequencies"])
    print("Frequencies (Hz):", output_freq_grid["frequencies"])

    (KQ, KQ_compressed_real, KQ_compressed_imag) = THDF_read_KQ_grid_from_hdf5(file_name)

    print("KQ: \n", KQ)
    print("KQ: \n", KQ_compressed_real)
    print("KQ: \n", KQ_compressed_imag)

    (JKQ_re, JKQ_im,
     JKQ_norm_re, JKQ_norm_im) = THDF_read_JKQ_from_hdf5(file_name, i, j, k)

    print("JKQ_re: ", JKQ_re)
    print("JKQ_im: ", JKQ_im)

    print("JKQ_norm_re: ", JKQ_norm_re)
    print("JKQ_norm_im: ", JKQ_norm_im)

    (KQ_real, KQ_imag) = THDF_read_JKQ_denormalize(KQ,
                                                   KQ_compressed_real,
                                                   KQ_compressed_imag,
                                                   JKQ_re,
                                                   JKQ_im,
                                                   JKQ_norm_re,
                                                   JKQ_norm_im)
    
    print ("KQ_real: ", KQ_real)
    print ("KQ_imag: ", KQ_imag)
    
    plot_KQ(KQ_real, KQ_imag, output_freq_grid["frequencies"], file_name, i, j, k)
