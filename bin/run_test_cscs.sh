#!/bin/bash -l

#SBATCH --ntasks=512
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --account=u2
#SBATCH --job-name="TRIP_PRD_3D"
# #SBATCH --qos=gp_ehpc
#SBATCH --partition=debug
#SBATCH --exclusive
#SBATCH --ntasks-per-socket=46

export APP_PATH=${HOME}/git/TRIP/build
export TESTS_DIR=${HOME}/git/TRIP/tests
export BIN_DIR=${HOME}/git/TRIP/bin

export SCRIPT_DIR=${HOME}/git/TRIP/bin

export MPICH_GPU_SUPPORT_ENABLED=1

echo -e "\n\033[1;34m========== Running Test 1 (CRD with B = 0) ==========\033[0m"
srun --cpu-bind=socket ${SCRIPT_DIR}/mps-wrapper.sh ${APP_PATH}/solar_3D_tests --yaml_config ${TESTS_DIR}/config_test_5x5_CRD_B0.yml 
cd ${TESTS_DIR}/test_output
rm test_profiles*.csv
cd ${BIN_DIR}

echo -e "\n\033[1;34m========== Running Test 2 (CRD with uniform B) ==========\033[0m"
srun --cpu-bind=socket ${SCRIPT_DIR}/mps-wrapper.sh ${APP_PATH}/solar_3D_tests --yaml_config ${TESTS_DIR}/config_test_5x5_CRD_B.yml
cd ${TESTS_DIR}/test_output
rm test_profiles*.csv
cd ${BIN_DIR}

echo -e "\n\033[1;34m========== Running Test 3 (PRD with B = 0) ==========\033[0m"
srun --cpu-bind=socket ${SCRIPT_DIR}/mps-wrapper.sh ${APP_PATH}/solar_3D_tests --yaml_config ${TESTS_DIR}/config_test_5x5_PRD_B0.yml
cd ${TESTS_DIR}/test_output
rm test_profiles*.csv
cd ${BIN_DIR}
