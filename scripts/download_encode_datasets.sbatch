#!/bin/bash
#SBATCH --job-name=download_datasets
#SBATCH --time=6:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --image=docker:argrosso/htspreprocessing:0.1.2
#SBATCH --workdir=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/TRF2_siRNA/

grep --perl-regexp --only-matching  'http\S*' encode_hela_chipseq_metadata_GRCh38_bigbedsandwigs.tsv > links.txt

mkdir encode_files

srun bash -c "cd encode_files; xargs -n 1 curl -O -L < ../links.txt"

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

