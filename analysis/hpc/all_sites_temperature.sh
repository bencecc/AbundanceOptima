#!/usr/bin/bash

my_container=/home/lisandro/containers/
my_scripts=/home/lisandro/workspace/MPA_Timeseries
my_indices=/home/lisandro/Lavori/MPA_timeseries/Modskurt

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

#echo "apptainer exec --bind /home/lisandro:/home/lisandro ${my_container}/stats.sif Rscript ${my_scripts}/all_sites_temperature.R ${P1} ${P2}" 
apptainer exec --bind /home/lisandro:/home/lisandro ${my_container}/stats.sif Rscript ${my_scripts}/all_sites_temperature.R ${P1} ${P2} 

EOF
`
echo "JobID = ${JOB} for parameters: ${P1} ${P2}" 
#submitted on `date`"
#echo SLURM_JOB_NODELIST is $SLURM_JOB_NODELIST

done < ${my_indices}/id.all.sites.temperature.txt
