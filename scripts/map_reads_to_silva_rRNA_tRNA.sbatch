#!/bin/bash
#SBATCH --job-name=rRNA-bowtie2
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=32
#SBATCH --image=docker:araposo/tophat:latest
#SBATCH --workdir=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/TRF2_siRNA/


srun -n1 -N1 --exclusive shifter bowtie2-build rRNA_tRNA_reference/hg38_rRNA_tRNA.fa rRNA_tRNA_reference/hg38_rRNA_tRNA

OUTDIR=bowtie2_to_rRNA_tRNA

mkdir $OUTDIR

for f in filtered_fastq/*1_filtered_paired.fq; do

filename=$(basename "${f%.*}").fq

srun -n1 -N1 --exclusive shifter bash -c "bowtie2 --very-sensitive-local --threads 32 \
                                         -x rRNA_tRNA_reference/hg38_rRNA_tRNA -1 $f -2 ${f//_1_filtered/_2_filtered} \
                                         --met-file $OUTDIR/${filename//.fq/_bowtie_metrics.txt} \
                                         --al-conc-gz $OUTDIR/${filename//.fq/maps_to_rRNA_tRNA} \
                                         --un-conc-gz $OUTDIR/${filename//.fq/not_rRNA_tRNA} \
                                         --al-gz $OUTDIR/${filename//.fq/singlets_maps_to_rRNA_tRNA} \
                                         --un-gz $OUTDIR/${filename//.fq/singlets_not_rRNA_tRNA} | samtools view -Sb -@ 32 - > $OUTDIR/${filename//.fq/_mapped_to_rRNA_tRNA.bam}" &
                                                 
done
wait

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID


