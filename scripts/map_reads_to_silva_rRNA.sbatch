#!/bin/bash
#SBATCH --job-name=rRNA-bowtie2
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=32
#SBATCH --image=docker:araposo/tophat:latest
#SBATCH --workdir=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/TRF2_siRNA/

cat arb-silva.de_2018-10-18_id563046_tax_silva.fasta arb-silva.de_2018-10-18_id563048_tax_silva.fasta > silva_Hsapiens_rRNA.fa

sed -i 's/U/T/g' silva_Hsapiens_rRNA.fa

srun -n1 -N1 --exclusive shifter bowtie2-build rRNA_reference/silva_Hsapiens_rRNA.fa rRNA_reference/silva_Hsapiens_rRNA

mkdir bowtie2_to_rRNA

for f in filtered_fastq/*1_filtered_paired.fq.gz; do

filename=$(basename "${f%.*}").gz

srun -n1 -N1 --exclusive shifter bash -c "bowtie2 --very-sensitive-local --threads 32 \
                                         -x rRNA_reference/silva_Hsapiens_rRNA -1 $f -2 ${f//_1_filtered/_2_filtered} \
                                         --met-file bowtie2_to_rRNA/${filename//.fq.gz/_bowtie_metrics.txt} \
                                         --al-conc-gz bowtie2_to_rRNA/${filename//.fq.gz/maps_to_rRNA} \
                                         --un-conc-gz bowtie2_to_rRNA/${filename//.fq.gz/not_rRNA} \
                                         --al-gz bowtie2_to_rRNA/${filename//.fq.gz/singlets_maps_to_rRNA} \
                                         --un-gz bowtie2_to_rRNA/${filename//.fq.gz/singlets_not_rRNA} | samtools view -Sb -@ 32 - > bowtie2_to_rRNA/${filename//.fq.gz/_mapped_to_rRNA.bam}" &
                                                 
done
wait

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID


