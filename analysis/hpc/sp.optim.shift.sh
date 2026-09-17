#!/usr/bin/bash

# ---- EDIT THESE THREE FOR YOUR CLUSTER ----
my_container=${HOME}/containers
my_scripts=${HOME}/workspace/MPA_Timeseries
my_indices=${HOME}/data/Modskurt   # must match config.R dir_data
while IFS=$'\t' read P1 P2

do
JOB=`sbatch << EOF
#!/usr/bin/bash
#SBATCH -p normal
#SBATCH --job-name=optim.shift.test
#SBATCH --output=optim.shift_%A_%a.ou
#SBATCH --error=optim.shift_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --cpus-per-task=1

#echo "apptainer exec --bind ${HOME}:${HOME} ${my_container}/stats.sif Rscript ${my_scripts}/sp_optimloc_ecoreg_shift.R ${P1} ${P2}" 
apptainer exec --bind ${HOME}:${HOME} ${my_container}/stats.sif Rscript ${my_scripts}/sp_optimloc_ecoreg_shift.R ${P1} ${P2} 

EOF
`
echo "JobID = ${JOB} for parameters: ${P1} ${P2}" 
#submitted on `date`"
#echo SLURM_JOB_NODELIST is $SLURM_JOB_NODELIST

done < ${my_indices}/id.optim.abund.relocated.txt
