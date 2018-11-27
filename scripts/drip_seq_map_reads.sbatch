#!/bin/bash
#SBATCH --job-name=rRNA-bowtie2
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=32
#SBATCH --image=docker:araposo/tophat:latest
#SBATCH --workdir=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/TRF2_siRNA/


while read line; do

array=($line)
SRA_ID=${array[1]}

srun -n1 -N1 --exclusive shifter --image=docker:araposo/tophat:latest bash -c "bowtie2 --very-sensitive-local --threads 32 \
                                         -x /mnt/nfs/lobo/IMM-NFS/genomes/hg38/Sequence/Bowtie2Index/genome -1 drip_seq_raw_data/${SRA_ID}_1.fastq -2 drip_seq_raw_data/${SRA_ID}_2.fastq \
                                         --met-file drip_seq_raw_data/${SRA_ID}_bowtie_metrics.txt \
                                         | samtools view -Sb -@ 32 - > drip_seq_raw_data/${SRA_ID}.bam" &
                                                 
done < drip_seq_raw_data/datasets.txt
wait

for f in drip_seq_raw_data/*bam; do

filename=$(basename "${f%.*}")

srun -n1 -N1 --exclusive shifter --image=docker:argrosso/htstools:0.2.1 macs2 filterdup -i $f --keep-dup=1 \
                                                                              -o drip_seq_raw_data/$filename.bed &

done
wait


echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID


