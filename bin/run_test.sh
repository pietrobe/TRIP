#!/bin/bash -l

#SBATCH --ntasks=2048
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --account=ehpc597
#SBATCH --job-name="TRIP_PRD_3D"
#SBATCH --qos=gp_ehpc
#SBATCH --exclusive

export APP_PATH=${HOME}/git/TRIP/build
export TESTS_DIR=${HOME}/git/TRIP/tests

echo -e "\n\033[1;34m========== Running Test 1 (CRD with B = 0) ==========\033[0m"
srun ${APP_PATH}/solar_3D_tests --yaml_config ${TESTS_DIR}/config_test_5x5_CRD_B0.yml 
rm ${TESTS_DIR}/test_output/test_profiles*.csv

echo -e "\n\033[1;34m========== Running Test 2 (CRD with uniform B) ==========\033[0m"
srun ${APP_PATH}/solar_3D_tests --yaml_config ${TESTS_DIR}/config_test_5x5_CRD_B.yml
rm ${TESTS_DIR}/test_output/test_profiles*.csv

echo -e "\n\033[1;34m========== Running Test 3 (CRD-JKQ with uniform B and no eps_csc) ==========\033[0m"
srun ${APP_PATH}/solar_3D_tests --yaml_config ${TESTS_DIR}/config_test_5x5_CRD_B_JKQ.yml
rm ${TESTS_DIR}/test_output/test_profiles*.csv

echo -e "\n\033[1;34m========== Running Test 4 (PRD with uniform B) ==========\033[0m"
if [ "$PETSC_PRECISION" = "single" ]; then
    echo -e "\033[1;33mWARNING: ksp_rtol might be too strict for single precision. PETSc precision: $PETSC_PRECISION\033[0m"
fi
srun ${APP_PATH}/solar_3D_tests --yaml_config ${TESTS_DIR}/config_test_5x5_PRD_B.yml
rm ${TESTS_DIR}/test_output/test_profiles*.csv
