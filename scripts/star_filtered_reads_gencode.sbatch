#!/bin/bash
#SBATCH --job-name=star
#SBATCH --time=60:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=40
#SBATCH --image=docker:argrosso/star:2.6.0b
#SBATCH --workdir=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/TRF2_siRNA/


OUTPUTDIR=star_on_hg38_output_gencode
mkdir $OUTPUTDIR

REFERENCEDIR=/mnt/nfs/lobo/IMM-NFS/genomes/hg38/Sequence/StarIndex/

LNCRNASANNOT=transcriptome_reference/gencode.v29.annotation.gtf

for f in filtered_fastq/*1_filtered_paired.fq; do

filename=$(basename "${f%.*}").fq

srun -n1 -N1 --exclusive shifter /software/STAR-2.6.0b/bin/Linux_x86_64/STAR --runThreadN 40 \
                   --outFileNamePrefix $OUTPUTDIR/${filename//_1_filtered_paired.fq/} \
                   --sjdbGTFfile $LNCRNASANNOT \
                   --genomeDir $REFERENCEDIR \
                   --readFilesIn filtered_fastq/$filename \
                   filtered_fastq/${filename//_1_filtered/_2_filtered} \
                  --outFilterType BySJout \
                  --outFilterMultimapNmax 20 \
                  --alignSJoverhangMin 8 \
                  --alignSJDBoverhangMin 1 \
                  --outFilterMismatchNmax 999 \
                  --outFilterMismatchNoverReadLmax 0.04 \
                  --alignIntronMin 20 \
                  --alignIntronMax 1000000 \
                  --alignMatesGapMax 1000000 \
                  --quantMode GeneCounts &
done
wait


echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

