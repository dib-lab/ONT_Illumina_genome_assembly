#!/bin/bash -l
#SBATCH -D /pylon5/mc5phkp/ljcohen/kfish_pilon/Axen
#SBATCH -J Axen_bwa
#SBATCH -o /pylon5/mc5phkp/ljcohen/kfish_pilon/Axen/bwa-Axen%j.o
#SBATCH -e /pylon5/mc5phkp/ljcohen/kfish_pilon/Axen/bwa-Axen%j.o
#SBATCH -t 24:00:00
#SBATCH -p LM
#SBATCH --mem=500GB
#SBATCH -c 32
module load samtools/1.9 bwa/0.7.13 
bwa index Axen.ctg.fa
bwa mem -t 32 Axen.ctg.fa Axen_L6L7_1.qc.fq.gz Axen_L6L7_2.qc.fq.gz | samtools view -Sb - | samtools sort - -o Axen_pilon.sorted.bam
