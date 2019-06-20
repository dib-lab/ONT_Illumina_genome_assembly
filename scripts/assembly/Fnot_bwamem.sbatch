#!/bin/bash -l
#SBATCH -D /pylon5/mc5phkp/ljcohen/kfish_pilon/Fnot
#SBATCH -J Fnot_bwa
#SBATCH -o /pylon5/mc5phkp/ljcohen/kfish_pilon/Fnot/bwa-Fnot%j.o
#SBATCH -e /pylon5/mc5phkp/ljcohen/kfish_pilon/Fnot/bwa-Fnot%j.o
#SBATCH -t 24:00:00
#SBATCH -p LM
#SBATCH --mem=500GB
#SBATCH -c 32
module load samtools/1.9 bwa/0.7.13 
bwa index Fnot.ctg.fa
bwa mem -t 32 Fnot.ctg.fa Fnot_L6L7_1.qc.fq.gz Fnot_L6L7_2.qc.fq.gz | samtools view -Sb - | samtools sort - -o Fnot_pilon.sorted.bam
