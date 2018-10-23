#!/bin/bash
#SBATCH --job-name=rRNA-sortmerna
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --ntasks=18
#SBATCH --cpus-per-task=1
#SBATCH --workdir=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/TRF2_siRNA/

mkdir sortmerna_to_rRNA

RRNADATABASES=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/TRF2_siRNA/sortmerna-2.1b/rRNA_databases
RRNAINDEXES=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/TRF2_siRNA/sortmerna-2.1b/index

#srun -n1 -N1 --exclusive ./sortmerna-2.1b/indexdb_rna --ref $RRNADATABASES/silva-bac-16s-id90.fasta,$RRNAINDEXES/silva-bac-16s-db:\
#$RRNADATABASES/silva-bac-23s-id98.fasta,$RRNAINDEXES/silva-bac-23s-db:\
#$RRNADATABASES/silva-arc-16s-id95.fasta,$RRNAINDEXES/silva-arc-16s-db:\
#$RRNADATABASES/silva-arc-23s-id98.fasta,$RRNAINDEXES/silva-arc-23s-db:\
#$RRNADATABASES/silva-euk-18s-id95.fasta,$RRNAINDEXES/silva-euk-18s-db:\
#$RRNADATABASES/silva-euk-28s-id98.fasta,$RRNAINDEXES/silva-euk-28s:\
#$RRNADATABASES/rfam-5s-database-id98.fasta,$RRNAINDEXES/rfam-5s-db:\
#$RRNADATABASES/rfam-5.8s-database-id98.fasta,$RRNAINDEXES/rfam-5.8s-db

for f in filtered_fastq/*1_filtered_paired.fq; do

filename=$(basename "${f%.*}")

#srun -n1 -N1 --exclusive bash -c "./sortmerna-2.1b/scripts/merge-paired-reads.sh $f ${f//_1_filtered/_2_filtered} sortmerna_to_rRNA/${filename//_1_filtered/_merged_filtered}
#                                  ./sortmerna-2.1b/sortmerna --ref $RRNADATABASES/silva-bac-16s-id90.fasta,$RRNAINDEXES/silva-bac-16s-db:\
#$RRNADATABASES/silva-bac-23s-id98.fasta,$RRNAINDEXES/silva-bac-23s-db:\
#$RRNADATABASES/silva-arc-16s-id95.fasta,$RRNAINDEXES/silva-arc-16s-db:\
#$RRNADATABASES/silva-arc-23s-id98.fasta,$RRNAINDEXES/silva-arc-23s-db:\
#$RRNADATABASES/silva-euk-18s-id95.fasta,$RRNAINDEXES/silva-euk-18s-db:\
#$RRNADATABASES/silva-euk-28s-id98.fasta,$RRNAINDEXES/silva-euk-28s:\
#$RRNADATABASES/rfam-5s-database-id98.fasta,$RRNAINDEXES/rfam-5s-db:\
#$RRNADATABASES/rfam-5.8s-database-id98.fasta,$RRNAINDEXES/rfam-5.8s-db \
#--reads sortmerna_to_rRNA/${filename//_1_filtered/_merged_filtered} \
#--aligned sortmerna_to_rRNA/${filename//_1_filtered/_rRNA_merged_filtered} \
#--other sortmerna_to_rRNA/${filename//_1_filtered/_non_rRNA_merged_filtered} \
#--fastx \
#--log \
#--paired_in \
#-a 32
srun -n1 -N1 --exclusive                                  ./sortmerna-2.1b/scripts/unmerge-paired-reads.sh sortmerna_to_rRNA/${filename//_1_filtered/_non_rRNA_merged_filtered}.fastq \
                                                                                    sortmerna_to_rRNA/${filename//_1_filtered/_non_rRNA_1_filtered}.fq \
                                                                                    sortmerna_to_rRNA/${filename//_1_filtered/_non_rRNA_2_filtered}.fq &
#                                   rm sortmerna_to_rRNA/${filename//_1_filtered/_merged_filtered}
#                                   rm sortmerna_to_rRNA/${filename//_1_filtered/_non_rRNA_merged_filtered.fastq}" &
done
wait

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID


