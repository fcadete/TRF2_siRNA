#!/bin/bash
#SBATCH --job-name=rRNA-bowtie2
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=32G
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --image=docker:argrosso/htstools:0.2.1
#SBATCH --workdir=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/TRF2_siRNA/


srun -n1 -N1 --exclusive shifter macs2 callpeak -t drip_seq_raw_data/SRR5164421.bam -c drip_seq_raw_data/SRR5164422.bam \
                                                -B --SPMR --outdir drip_seq_raw_data/DRIP_seq_1 -n DRIP_seq_1 &


srun -n1 -N1 --exclusive shifter macs2 callpeak -t drip_seq_raw_data/SRR5683211.bam -c drip_seq_raw_data/SRR5683212.bam \
                                                -B --SPMR --outdir drip_seq_raw_data/DRIP_seq_2 -n DRIP_seq_2 &

wait

srun -n1 -N1 --exclusive shifter macs2 bdgcmp -t drip_seq_raw_data/DRIP_seq_1/DRIP_seq_1_treat_pileup.bdg \
                                              -c drip_seq_raw_data/DRIP_seq_1/DRIP_seq_1_control_lambda.bdg \
                                              -o drip_seq_raw_data/DRIP_seq_1/DRIP_seq_1_FE.bdg -m FE &

srun -n1 -N1 --exclusive shifter macs2 bdgcmp -t drip_seq_raw_data/DRIP_seq_1/DRIP_seq_1_treat_pileup.bdg \
                                              -c drip_seq_raw_data/DRIP_seq_1/DRIP_seq_1_control_lambda.bdg \
                                              -o drip_seq_raw_data/DRIP_seq_1/DRIP_seq_1_logLR.bdg -m logLR -p 0.00001 &

srun -n1 -N1 --exclusive shifter macs2 bdgcmp -t drip_seq_raw_data/DRIP_seq_2/DRIP_seq_2_treat_pileup.bdg \
                                              -c drip_seq_raw_data/DRIP_seq_2/DRIP_seq_2_control_lambda.bdg \
                                              -o drip_seq_raw_data/DRIP_seq_2/DRIP_seq_2_FE.bdg -m FE &

srun -n1 -N1 --exclusive shifter macs2 bdgcmp -t drip_seq_raw_data/DRIP_seq_2/DRIP_seq_2_treat_pileup.bdg \
                                              -c drip_seq_raw_data/DRIP_seq_2/DRIP_seq_2_control_lambda.bdg \
                                              -o drip_seq_raw_data/DRIP_seq_2/DRIP_seq_2_logLR.bdg -m logLR -p 0.00001 &

wait

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID


