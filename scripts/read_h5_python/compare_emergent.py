import argparse
import os
import matplotlib.pyplot as plt
import numpy as np

# Import helper functions from read_emergent.py
from read_emergent import (
    THDF_read_angular_grid_from_hdf5,
    THDF_read_field_from_hdf5,
    THDF_read_frequencies_grid_from_hdf5,
    resolve_cli_direction_indices,
)

CONFIGS = [ 
    # {
    #     "name": "CaI_B0_TwoTerm_routines",
    #     # List of HDF5 file paths to compare
    #     "files": [
    #         "/home/melanie-tonarelli/Scrivania/SolarCode/data/output_CaI_v2/FAL-C-B0.PRD_TWOTERM/emergent_field_angular_grid.h5",
    #         "/home/melanie-tonarelli/Scrivania/SolarCode/data/output_CaI_v2/FAL-C-B0.CRD_TWOTERM_limit/emergent_field_angular_grid.h5",
    #         "/home/melanie-tonarelli/Scrivania/SolarCode/data/output_CaI_v2/FAL-C-B0.PRD_AA_TWOTERM/emergent_field_angular_grid.h5",
    #     ],
    #     "labels": [
    #         "PRD two-term", 
    #         "CRD two-term", 
    #         "PRDAA two-term",
    #     ],
    #     "save": "CaI_B0_TwoTerm_routines.png",  # Base output filename or None
    #     # "output_dir": "plots/output",      # Relative or absolute folder path (optional)
    #     "x": 0,
    #     "y": 0,
    #     "incl_list": [0, 1, 2, 4],  # List of inclination indices to iterate over
    #     "azim": 1,
    #     "xlim": [4225, 4229],
    # },
    # {
    #     "name": "CaI_B20_TwoTerm_routines",
    #     # List of HDF5 file paths to compare
    #     "files": [
    #         "/home/melanie-tonarelli/Scrivania/SolarCode/data/output_CaI_v2/FAL-C-B20.PRD_TWOTERM/emergent_field_angular_grid.h5",
    #         "/home/melanie-tonarelli/Scrivania/SolarCode/data/output_CaI_v2/FAL-C-B20.CRD_TWOTERM_limit/emergent_field_angular_grid.h5",
    #         "/home/melanie-tonarelli/Scrivania/SolarCode/data/output_CaI_v2/FAL-C-B20.PRD_AA_TWOTERM/emergent_field_angular_grid.h5",
    #     ],
    #     "labels": [
    #         "PRD two-term", 
    #         "CRD two-term", 
    #         "PRDAA two-term",
    #     ],
    #     "save": "CaI_B20_TwoTerm_routines.png",  # Base output filename or None
    #     # "output_dir": "plots/output",      # Relative or absolute folder path (optional)
    #     "x": 0,
    #     "y": 0,
    #     "incl_list": [0, 1, 2, 4],  # List of inclination indices to iterate over
    #     "azim": 1,
    #     "xlim": [4225, 4229],
    # },
    {
        "name": "CaI_unmag_comp",
        # List of HDF5 file paths to compare
        "files": [
            # "/home/melanie-tonarelli/Scrivania/SolarCode/data/output_CaI/FAL-C.CRD_TWOTERM/emergent_field_angular_grid.h5",
            "/home/melanie-tonarelli/Scrivania/SolarCode/data/output_CaI_v3/FAL-C.PRD_AA_GB/emergent_field_angular_grid.h5",
            "/home/melanie-tonarelli/Scrivania/SolarCode/data/output_CaI_v3/FAL-C-B0.PRD_AA_TWOTERM/emergent_field_angular_grid.h5",
            "/home/melanie-tonarelli/Scrivania/SolarCode/data/output_CaI_v3/FAL-C-B0.PRD_AA_TWOTERM_unmag/emergent_field_angular_grid.h5",
        ],
        "labels": [
            # "PRD two-term", 
            # "CRD two-term (no limit)", 
            "PRDAA GB",
            "PRDAA (with B=0)",
            "PRDAA unmag"
        ],
        "save": "CaI_unmag_comp.png",  # Base output filename or None
        # "output_dir": "plots/output",      # Relative or absolute folder path (optional)
        "x": 0,
        "y": 0,
        "incl_list": [0, 1, 2, 4],  # List of inclination indices to iterate over
        "azim": 1,
        "xlim": [4225, 4229],
    },
]

LINESTYLES = ["-", "--", "-.", ":", (0, (3, 1, 1, 1)), (0, (5, 1))]


def create_single_plot(incl_index, params, output_dir=None):
    files = params["files"]
    labels = params["labels"]
    x_idx = params["x"]
    y_idx = params["y"]
    xlim = params["xlim"]
    azim_idx = params["azim"]

    fig, axs = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    ax = axs.ravel()

    plot_keys = ["stokes_I", "stokes_QI", "stokes_UI", "stokes_VI"]
    ylabels = ["Stokes I", "Stokes Q/I (%)", "Stokes U/I (%)", "Stokes V/I (%)"]

    colors = plt.cm.tab10(np.linspace(0, 1, max(10, len(files))))
    legend_handles = []

    extracted_mu = None
    extracted_theta = None

    for idx, file_path in enumerate(files):
        if not os.path.exists(file_path):
            print(f"Warning: File not found: {file_path}. Skipping...")
            continue

        freq_grid = THDF_read_frequencies_grid_from_hdf5(file_path)
        freqs = freq_grid["frequencies"]

        try:
            from wavelength_utils import freq_to_angstrom, vacuum_to_air
            wl_air = vacuum_to_air(freq_to_angstrom(freqs))
        except ImportError:
            wl_air = freqs

        try:
            angular_grid = THDF_read_angular_grid_from_hdf5(file_path)
            n_azim = angular_grid["N_azimuthal_angles"]
        except KeyError:
            angular_grid = None
            n_azim = None

        incl_i, azim_i, _ = resolve_cli_direction_indices(
            incl_index, azim_idx, angular_grid
        )
        if angular_grid is not None and "inclination_angles" in angular_grid:
            extracted_theta = float(angular_grid["inclination_angles"][incl_i])
            extracted_mu = float(np.cos(extracted_theta))

        field = THDF_read_field_from_hdf5(
            file_path,
            N_frequencies=freq_grid["N_frequencies"],
            x_i=x_idx,
            y_i=y_idx,
            incl_i=incl_i,
            azimuth_i=azim_i,
            N_azimuth=n_azim,
        )

        file_label = labels[idx] if labels else os.path.basename(file_path)
        color = colors[idx % len(colors)]
        linestyle = LINESTYLES[idx % len(LINESTYLES)]

        for panel, key in enumerate(plot_keys):
            line = ax[panel].plot(
                wl_air,
                field[key],
                color=color,
                label=file_label,
                linestyle=linestyle,
                linewidth=1.5,
            )[0]
            if panel == 0:
                legend_handles.append(line)

    for panel, ylabel in enumerate(ylabels):
        ax[panel].set_ylabel(ylabel)
        ax[panel].grid(True, linestyle="--", alpha=0.6)
        if xlim is not None:
            ax[panel].set_xlim(xlim[0], xlim[1])
        if panel >= 2:
            ax[panel].set_xlabel("Wavelength (Å)" if "wl_air" in locals() else "Frequency / Axis Index")

    ax[0].legend(handles=legend_handles, loc="lower right", fontsize=8, framealpha=0.9)

    title_text = f"Stokes vs Wavelength (x={x_idx}, y={y_idx}, incl_idx={incl_index})"
    if extracted_mu is not None:
        title_text += rf"  |  $\mu = {extracted_mu:.4f}$ ($\theta = {extracted_theta:.4f}$ rad)"

    fig.suptitle(title_text, fontsize=13, fontweight="bold", color="darkblue")
    fig.tight_layout()
    fig.subplots_adjust(top=0.88)

    if params["save"]:
        base_name, ext = os.path.splitext(params["save"])
        filename = f"{base_name}_incl_{incl_index}{ext}"

        # 1. Determine folder (CLI argument > config output_dir > current folder)
        raw_target_dir = output_dir or params.get("output_dir", ".")

        # 2. Convert to an absolute path & expand '~' if used (e.g., ~/Desktop)
        target_dir = os.path.abspath(os.path.expanduser(raw_target_dir))

        # 3. Create folder (and intermediate folders) if missing
        os.makedirs(target_dir, exist_ok=True)

        # 4. Save file
        out_path = os.path.join(target_dir, filename)
        fig.savefig(out_path, dpi=300)
        print(f"Saved: {out_path}")
        plt.close(fig)
    else:
        plt.show()


def run_single_config(cfg, output_dir=None):
    params = cfg.copy()

    if "incl_list" in params and params["incl_list"]:
        incl_list = params["incl_list"]
    else:
        incl_list = [params.get("incl", 0)]

    for incl_idx in incl_list:
        print(f"\n--- Generating plot for inclination index = {incl_idx} ---")
        create_single_plot(incl_index=incl_idx, params=params, output_dir=output_dir)


def main():
    parser = argparse.ArgumentParser(description="Plot emergent field Stokes parameters.")
    parser.add_argument(
        "-o", "--output-dir",
        type=str,
        default=None,
        help="Directory path where PNG plots will be saved (relative or absolute)."
    )
    args = parser.parse_args()

    print(f"Found {len(CONFIGS)} configuration sets to process.\n")
    for idx, cfg in enumerate(CONFIGS, 1):
        cfg_name = cfg.get("name", f"Set_{idx}")
        print("=" * 45)
        print(f" Running [{idx}/{len(CONFIGS)}]: {cfg_name}")
        print("=" * 45)
        run_single_config(cfg, output_dir=args.output_dir)
        print()


if __name__ == "__main__":
    main()