#!/bin/bash
#SBATCH --job-name=tophat
#SBATCH --time=60:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=40
#SBATCH --image=docker:araposo/tophat:latest
#SBATCH --workdir=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/TRF2_siRNA/

mkdir tophat_on_hg19_output

srun -n1 -N1 --exclusive shifter tophat --GTF=transcriptome_reference/hg38_gencode24.gtf \
                                        --transcriptome-index=transcriptome_reference/transcriptome-index \
                                        /mnt/nfs/lobo/IMM-NFS/genomes/hg38/Sequence/Bowtie2Index/genome

for f in filtered_fastq/*1_filtered_paired.fq; do

filename=$(basename "${f%.*}").fq

srun -n1 -N1 --exclusive shifter tophat -p 40 -o tophat_on_hg19_output/${filename//_1_filtered_paired.fq/} \
                   --GTF=transcriptome_reference/hg38_gencode24.gtf \
                   --transcriptome-index=transcriptome_reference/transcriptome-index \
                   /mnt/nfs/lobo/IMM-NFS/genomes/hg38/Sequence/Bowtie2Index/genome \
                   filtered_fastq/$filename \
                   filtered_fastq/${filename//_1_filtered/_2_filtered} &
done
wait


echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

