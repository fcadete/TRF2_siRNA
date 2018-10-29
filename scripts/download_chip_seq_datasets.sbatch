#!/bin/bash
#SBATCH --job-name=download_datasets
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --image=docker:argrosso/htspreprocessing:0.1.2
#SBATCH --workdir=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/TRF2_siRNA/


while read line; do
  array=($line)
  SRA_ID=${array[1]}
  srun --exclusive -n1 -N1 shifter fastq-dump -O chip_seq_raw_data $SRA_ID &
done < chip_seq_raw_data/datasets.txt
wait

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

