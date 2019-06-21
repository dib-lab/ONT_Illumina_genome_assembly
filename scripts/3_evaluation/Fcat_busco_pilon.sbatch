#!/bin/bash -l
#SBATCH -D /pylon5/mc5phkp/ljcohen/kfish_pilon/Fcat/
#SBATCH -J Fcat_busco
#SBATCH -o /pylon5/mc5phkp/ljcohen/kfish_pilon/Fcat/busco-Fcat%j.o
#SBATCH -e /pylon5/mc5phkp/ljcohen/kfish_pilon/Fcat/busco-Fcat%j.o
#SBATCH -t 6:00:00
#SBATCH -p LM
#SBATCH --mem=200GB
#SBATCH -c 32
source activate dammit_conda
export AUGUSTUS_CONFIG_PATH=/pylon5/mc5phkp/ljcohen/config
run_BUSCO.py --in Fcat_pilon.fasta -m genome --out Fcat_pilon --lineage_path /pylon5/bi5fpmp/ljcohen/dammit/busco2db/eukaryota_odb9
