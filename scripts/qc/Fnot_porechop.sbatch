#!/bin/bash -l
#SBATCH -D /pylon5/mc5phkp/ljcohen/kfish_ONT
#SBATCH -J Fnot-porechop
#SBATCH -o /pylon5/mc5phkp/ljcohen/kfish_ONT/porec-Fnot%j.o
#SBATCH -e /pylon5/mc5phkp/ljcohen/kfish_ONT/porec-Fnot%j.o
#SBATCH -t 12:00:00
#SBATCH -p LM
#SBATCH --mem=200GB
#SBATCH --ntasks-per-node 14
#SBATCH --cpus-per-task 1

source ~/.bashrc
source activate ONT
cd /pylon5/mc5phkp/ljcohen/kfish_ONT

porechop -i F_nottii.fq.gz -o F_nottii_trimmed.fq.gz
