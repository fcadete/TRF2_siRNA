#!/bin/bash
#SBATCH --job-name=tophat
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=1
#SBATCH --image=docker:argrosso/htstools:0.2.1
#SBATCH --workdir=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/TRF2_siRNA/

mkdir tophat_extracted_fastq

for f in filtered_fastq/*1_filtered_paired.fq; do

filename=$(basename "${f%.*}").fq
sample=${filename//_1_filtered_paired.fq/}

echo $sample

srun -n1 -N1 --exclusive shifter bash -c "samtools sort -n -o tophat_on_hg19_output/$sample/accepted_hits.qsort tophat_on_hg19_output/$sample/accepted_hits.bam
                                          bedtools bamtofastq -i tophat_on_hg19_output/$sample/accepted_hits.qsort \
                                                              -fq tophat_extracted_fastq/${sample}_1.fq \
                                                              -fq2 tophat_extracted_fastq/${sample}_2.fq" &
done
wait


echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

