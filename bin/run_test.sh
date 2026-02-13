#!/bin/bash -l

#SBATCH --ntasks=2048
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --account=ehpc238
#SBATCH --job-name="TRIP_PRD_3D"
#SBATCH --qos=gp_ehpc
#SBATCH --exclusive

export APP_PATH=${HOME}/git/TRIP/build
export TESTS_DIR=${HOME}/git/TRIP/tests


echo -e "\n\033[1;34m========== Running Test 1 (CRD with B = 0) ==========\033[0m"
srun ${APP_PATH}/solar_3D_tests --yaml_config ${TESTS_DIR}/config_test_5x5_CRD_B0.yml 

echo -e "\n\033[1;34m========== Running Test 2 (CRD with uniform B) ==========\033[0m"
srun ${APP_PATH}/solar_3D_tests --yaml_config ${TESTS_DIR}/tests/config_test_5x5_CRD_B.yml

echo -e "\n\033[1;34m========== Running Test 3 (PRD with uniform B) ==========\033[0m"
srun ${APP_PATH}/solar_3D_tests --yaml_config ${TESTS_DIR}/tests/config_test_5x5_PRD_B.yml
