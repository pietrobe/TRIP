#!/bin/bash -l

#SBATCH --ntasks=512
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --account=ehpc238
#SBATCH --job-name="TRIP_PRD_3D"
#SBATCH --qos=gp_ehpc
#SBATCH --exclusive

echo -e "\n\033[1;34m========== Running Test 1 (CRD with B = 0) ==========\033[0m"
srun ./solar_3D_tests --yaml_config ../tests/config_test_5x5_CRD_B0.yml 

echo -e "\n\033[1;34m========== Running Test 2 (CRD with uniform B) ==========\033[0m"
srun ./solar_3D_tests --yaml_config ../tests/config_test_5x5_CRD_B.yml
