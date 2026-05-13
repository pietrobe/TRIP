import os
import shutil
import glob
import re
import argparse

def get_id_from_file(file_path, prefix_pattern):
    """
    Helper function to find a specific pattern in a file and extract 
     the numeric ID following it.
    """
    try:
        with open(file_path, 'r') as f:
            for line in f:
                # Search for the prefix (case-insensitive) and capture the digits after it
                match = re.search(rf'{prefix_pattern}\s*(\d+)', line, re.IGNORECASE)
                if match:
                    return match.group(1).strip()
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
    return None

def process_slurm_files(search_directory):
    # 1. Find all files starting with "slurm-" in the given directory
    pattern = os.path.join(search_directory, "slurm-*")
    slurm_files = glob.glob(pattern)

    if not slurm_files:
        print(f"No files matching 'slurm-*' found in {search_directory}")
        return

    for file_path in slurm_files:
        filename = os.path.basename(file_path)
        
        # Skip files already processed (renamed to slurm_c-*)
        if filename.startswith("slurm_c-"):
            continue

        output_dir = None
        slurm_job_id = None
        
        # 2. Extract Output path and SLURM_JOB_ID from the current slurm file
        try:
            with open(file_path, 'r') as f:
                for line in f:
                    # Find the Output path
                    if 'Output path:' in line:
                        path_match = re.search(r'Output path:\s*"([^"]+)"', line)
                        if path_match:
                            output_dir = path_match.group(1)
                    
                    # Find the SLURM_JOB_ID inside the slurm file
                    if 'SLURM_JOB_ID:' in line:
                        id_match = re.search(r'SLURM_JOB_ID:\s*(\d+)', line)
                        if id_match:
                            slurm_job_id = id_match.group(1)
                            
                    # Break early if both pieces of info are found
                    if output_dir and slurm_job_id:
                        break
        except Exception as e:
            print(f"Error processing {filename}: {e}")
            continue

        # 3. Validation and ID Comparison
        if output_dir and slurm_job_id:
            info_file_path = os.path.join(output_dir, "info.txt")
            
            if not os.path.exists(info_file_path):
                print(f"Skipping {filename}: 'info.txt' not found in {output_dir}")
                continue

            # Extract ID from info.txt (Expected line: "SLURM Job ID: 12345")
            info_job_id = get_id_from_file(info_file_path, "SLURM Job ID:")

            if info_job_id and info_job_id == slurm_job_id:
                try:
                    # 4. Copy the file to the output path
                    dest_path = os.path.join(output_dir, filename)
                    shutil.copy2(file_path, dest_path)
                    print(f"Match found! Copied {filename} to {output_dir}")

                    # 5. Rename the original file in the search directory
                    new_filename = filename.replace("slurm-", "slurm_c-", 1)
                    new_file_path = os.path.join(search_directory, new_filename)
                    os.rename(file_path, new_file_path)
                    print(f"Renamed {filename} -> {new_filename}")

                except Exception as e:
                    print(f"Failed to copy/rename {filename}: {e}")
            else:
                print(f"ID mismatch for {filename}: Slurm file({slurm_job_id}) vs info.txt({info_job_id})")
        else:
            print(f"Missing info in {filename}: Path found? {bool(output_dir)}, ID found? {bool(slurm_job_id)}")

def parse_args():
    parser = argparse.ArgumentParser(
        description="Copy matching slurm output files into their job output directories."
    )
    parser.add_argument(
        "-d",
        "--slurm-dir",
        default="./",
        help="Directory containing the input slurm-* files. Defaults to the current directory.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    process_slurm_files(args.slurm_dir)
