#!/bin/bash
#SBATCH --job-name=fastQC
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --image=docker:argrosso/htspreprocessing:0.1.1
#SBATCH --workdir=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/TRF2_siRNA

mkdir filtered_fastqc_reports

for f in filtered_fastq/*fq.gz; do

filename=$(basename "${f%.*}")

mkdir filtered_fastqc_reports/$filename

srun --exclusive -n1 -N1 shifter bash -c "zcat $f | fastqc --outdir=filtered_fastqc_reports/$filename/ stdin" &

done

wait


echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID
