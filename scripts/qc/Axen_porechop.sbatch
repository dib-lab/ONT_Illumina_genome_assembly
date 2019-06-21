#!/bin/bash -l
#SBATCH -D /pylon5/mc5phkp/ljcohen/kfish_ONT
#SBATCH -J Axen-porechop
#SBATCH -o /pylon5/mc5phkp/ljcohen/kfish_ONT/porec-Axen%j.o
#SBATCH -e /pylon5/mc5phkp/ljcohen/kfish_ONT/porec-Axen%j.o
#SBATCH -t 12:00:00
#SBATCH -p LM
#SBATCH --mem=500GB
#SBATCH --ntasks-per-node 14
#SBATCH --cpus-per-task 1

source ~/.bashrc
source activate ONT
cd /pylon5/mc5phkp/ljcohen/kfish_ONT

porechop -i A_xenica_woLambda.fq.gz -o A_xenica_trimmed.fq.gz
