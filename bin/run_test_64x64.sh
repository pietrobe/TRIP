#!/bin/bash -l

#SBATCH --ntasks=12288
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --account=iac90
#SBATCH --job-name="PRD_3D"
#SBATCH --qos=gp_resa
### #SBATCH --nodes=110
#SBATCH --exclusive    
#SBATCH --constraint=highmem
## #SBATCH --mem=251G

### add the joj id to the output file name
#### #SBATCH --output=/home/usi/usi441290/ms5_jobs_log/PRD_3D_64x64_%j.out
### #SBATCH  --error=/home/usi/usi441290/ms5_jobs_log/PRD_3D_64x64_%j.err

#SBATCH --mail-type=all
#SBATCH --mail-user=simone.riva@usi.ch    

#### Per ricevere email di notifica ci si deve registrare in https://userportal.bsc.es/

srun solar_3D -ksp_monitor -ksp_monitor_true_residual  -ksp_view -ksp_rtol 1e-4 -ksp_max_it 1
 
