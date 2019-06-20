#!/bin/bash -l
#SBATCH -D /pylon5/mc5phkp/ljcohen/kfish_pilon/Foli
#SBATCH -J Foli_bwa
#SBATCH -o /pylon5/mc5phkp/ljcohen/kfish_pilon/Foli/bwa-Foli%j.o
#SBATCH -e /pylon5/mc5phkp/ljcohen/kfish_pilon/Foli/bwa-Foli%j.o
#SBATCH -t 24:00:00
#SBATCH -p LM
#SBATCH --mem=500GB
#SBATCH -c 32
module load samtools/1.9 bwa/0.7.13 
bwa index Foli.ctg.fa
bwa mem -t 32 Foli.ctg.fa Killi_Fish_S19_L001_R1_001.qc.fq.gz Killi_Fish_S19_L001_R2_001.qc.fq.gz | samtools view -Sb - | samtools sort - -o Foli_pilon.sorted.bam
