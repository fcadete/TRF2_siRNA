#!/bin/bash
#SBATCH --job-name=fastQC
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --image=docker:argrosso/htstools:0.2.1
#SBATCH --workdir=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/TRF2_siRNA


for f in encode_files/*bigBed; do

srun --exclusive -n1 -N1 shifter bigBedToBed $f ${f//bigBed/bed} &

done

wait


echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID
