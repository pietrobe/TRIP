#!/bin/bash -l

### to be defined, 12288 is tested !!!!!!!!!
##### #SBATCH --ntasks=16384
#SBATCH --ntasks=12288

## to be defined !!!!!!!!! 
#SBATCH --time=10:12:27

#### #SBATCH --cpus-per-task=1
#SBATCH --account=iac90
#SBATCH --job-name="TRIP_PRD_3D"
#SBATCH --qos=gp_resa

#SBATCH --mail-type=ALL
#SBATCH --mail-user=my.email

#### #SBATCH --constraint=highmem

###### #SBATCH --nodes=1


export CRD="--CRD"
export MAX_ITER=170

### MareNostrum 5
export MAIN_DATA_DIR=/gpfs/projects/iac90/

export SERIES_DIR=DataSet_TRIP_PORTA

### Set the horizontal grid size
export GRID_SIZE=64x64

### MareNostrum5
export INPUT_DIR=${MAIN_DATA_DIR}/input/${GRID_SIZE}
export OUTPUT_DIR=${MAIN_DATA_DIR}/output/main/${GRID_SIZE}   ### ATTENTION: test output directory

export RTOL=1e-8

######################################################
########### Jobs for the first serie of runs #########

## 3 TODO
## make for thest also with the classic RII grid (fast)
## For 64 x 64 grid size in the presence of Magnetic and bulk velocity fields
# export INPUT_CONFIG=AR_385_Cut_64x64_mirrorxy-CRD_I_B_V.conf

## 3.1 NOT TODO (non serve ormai il FAST e' sufficientemente verificato)
## make for thest also with the classic RII grid (fast)
## For 64 x 64 grid size in the presence of Magnetic and bulk velocity fields
## export INPUT_CONFIG=AR_385_Cut_64x64_mirrorxy-CRD_I_B_V.conf

#################

## 1 TODO
## For 64 x 64 grid size in the absence of Magnetic and bulk velocity fields
# export INPUT_CONFIG=AR_385_Cut_64x64_mirrorxy-CRD_I_B0_V0.conf

#################

## 2 TODO
## For 64 x 64 grid size in the presence of Magnetic and absence of bulk velocity fields
### export INPUT_CONFIG=AR_385_Cut_64x64_mirrorxy-CRD_I_B_V0.conf

#################

## 4 TODO
## For 64 x 64 grid size in the absence of Magnetic and presence of bulk velocity fields
# export INPUT_CONFIG=AR_385_Cut_64x64_mirrorxy-CRD_I_B0_V.conf
######################################################

echo "INPUT_DIR:     $INPUT_DIR"
echo "OUTPUT_DIR:    $OUTPUT_DIR"
echo "INPUT_CONFIG:  $INPUT_CONFIG"
echo "CRD:           $CRD"
echo "RTOL:          $RTOL"
echo ""
echo "Command to run: "
echo "/home/usi/usi441290/git/solar_3d/build/solar_3D $CRD --input_dir $INPUT_DIR  --problem_input_config $INPUT_CONFIG --output_dir $OUTPUT_DIR -ksp_type fgmres -ksp_gmres_restart 30 -ksp_max_it 22 -ksp_monitor -ksp_view -ksp_rtol $RTOL"

echo ""
echo "Running Solar 3D ........ "
echo ""
echo " ---------------------------- "
echo ""

srun /home/usi/usi441290/git/solar_3d/build/solar_3D $CRD --input_dir $INPUT_DIR  --problem_input_config $INPUT_CONFIG --output_dir $OUTPUT_DIR -ksp_type fgmres -ksp_gmres_restart 30 -ksp_max_it $MAX_ITER -ksp_monitor -ksp_view -ksp_rtol $RTOL




