#!/bin/bash
#SBATCH --job-name=rRNA-bowtie2
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=32G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --image=docker:argrosso/htstools:0.2.1
#SBATCH --workdir=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/TRF2_siRNA/


srun -n1 -N1 --exclusive shifter macs2 callpeak -t chip_seq_raw_data/SRR831105.bam -c chip_seq_raw_data/SRR831099.bam \
                                                --outdir chip_seq_raw_data/NR2F2


srun -n1 -N1 --exclusive shifter macs2 callpeak -t chip_seq_raw_data/SRR831106.bam -c chip_seq_raw_data/SRR831099.bam \
                                                --outdir chip_seq_raw_data/NR2C2


echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID


