Blog post (draft): Hybrid genome assembly of the killifish, *Fundulus olivaceus* with ONT and Illumina NovaSeq
=========================

# Thank you!

* C. Titus Brown and Moore Foundation for encouragement and funding to go to Porecamp!
* [Andrew Whitehead](https://whiteheadresearch.wordpress.com/) for motivating questions, and Tony Gil for assistance with samples and hmw DNA extraction: [Agilent protocol (Catalog #200600), Protocol II: Extraction from Whole Tissue](http://www.agilent.com/cs/library/usermanuals/public/200600.pdf)
* [TAMU seq facility](http://www.txgen.tamu.edu/), Charlie Johnson, Richard Metz, Joshua Hill. Hosted Porecamp, then they offered us free NovaSeq!
* [Porecamp](http://www.txgen.tamu.edu/porecamp_usa/): Nick Loman, Joshua Quick, Mick Watson, Matt Loose, John Tyson
* [David Duvernell, Southern Illinois State University](http://www.siue.edu/~dduvern/) for generously providing  *F. olivaceus* fish

# Purpose

To generate a genome assembly for each of several freshwater *Fundulus* killifish. In particular, *Fundulus olivaceus* data were collected for this analysis. Subsequent data will be collected for other freshwater *Fundulus* species. The final assemblies for each species will be used to compare genomic regions to brackish/marine congeners.

Want to be able to:

* Scan for regulatory elements, such as osmolality/salinity-responsive enhance (OSREs) ([Wang and Kültz 2017](http://www.pnas.org/content/114/13/E2729)).
* Detect gene duplications, deletions of regions
* Potentially help to resolve RNAseq assemblies, detect isoform switching, e.g. NKCC1 alpha1a and alpha1b subunit (assemblies with only short reads cannot resolve whether transcripts span multiple genomic regions):
https://link.springer.com/article/10.1007%2Fs00360-012-0719-y
http://ajpcell.physiology.org/content/ajpcell/287/2/C300.full.pdf
http://jeb.biologists.org/content/212/24/3994

Long reads can potentially help to solve problems with short-read assemblies.

# Questions
1. Do we have enough data from *F. olivaceus* to meet our goals to detect regulatory elements and isoforms? Do we need more ONT reads? More Illumina reads? Could we get away with fewer Illumina reads?
2. What is the most cost-efficient way to generate a reference genome when there is a closely-related sister species? ONT and Illumina? ONT only?

# Data Generation

* NovaSeq: 600 million reads of 2x150 (300b) = 180Gb
* 3 ONT flowcells, albacore base calling

```
stats for ../fastq2/porecamp_killifish2.fastq
sum = 4962626713, n = 740248, ave = 6704.01, largest = 973552
N50 = 12726, n = 117202
N60 = 10357, n = 160433
N70 = 8098, n = 214460
N80 = 5724, n = 286845
N90 = 3229, n = 400661
N100 = 5, n = 740248
N_count = 0
Gaps = 0
```

*Some* files are backed up on [OSF](https://osf.io/gek4p/). See [Instructions for installing osf client](http://osfclient.readthedocs.io/en/stable/).

```
export OSF_PASSWORD=password
export OSF_USERNAME=ljcohen@ucdavis.edu
# list files in OSF
osf -p zjv86 ls
# copy files from OSF to local (or hpc - wherever you're working)
osf -p zjv86 clone Folivaceus_hybrid_genome_assembly
osf -p zjv86 upload /mnt/scratch/ljcohen/onp_Porecamp_killifish/fastq2/porecamp_killifish_ONT_reads.fastq porecamp_killifish_ONT_reads.fastq 

```

12/18/2017: Need to backup ONT and NovaSeq raw files to a location for safekeeping. NCBI, S3?

Subset data for analysis runthrough.

```
cd 
head -2000000 F_olivaceus_ONT.fasta > subset_F_olivaceus_ONT_1m.fasta
zcat F_olivaceus_R1.fastq.gz | head -2000000 > subset_F_olivaceus_R1.fastq
zcat F_olivaceus_R2.fastq.gz | head -2000000 > subset_F_olivaceus_F2.fastq
```

## De novo genome assembly of *Fundulus olivaceus*
Generated ~5x coverage. With a 1.1Gb genome size for *F. olivaceus* (based on 1.1 Gb size of [*F. heteroclitus*](https://academic.oup.com/gbe/article/9/3/659/2992614/The-Landscape-of-Extreme-Genomic-Variation-in-the?guestAccessKey=8e254afa-b25d-4436-a4f8-8e54951d78f9)). 

Which methods to use? Could assembly long reads, then align short reads. Or, assemble short reads, then align contigs to long reads. Lots of tools to do this.

There are examples where long ONT reads improve Illumina assemblies. The [Murray cod genome](https://academic.oup.com/view-large/figure/95750508/gix063fig1.jpg) improved with low ONT coverage, 804 Mb combined with 70.6 Gb Illumina (HiSeq and MiSeq) ([Austin et al. 2017](https://academic.oup.com/gigascience/article/6/8/1/3979350/De-novo-genome-assembly-and-annotation-of)).

# Available Tools:

With varying degrees of difficulty installing and using:

### Assemblers:

* [Canu v1.6](https://github.com/marbl/canu/releases): [Quick start](http://canu.readthedocs.io/en/latest/quick-start.html), [tutorial](https://github.com/marbl/canu/blob/master/documentation/source/tutorial.rst) and [paper](http://genome.cshlp.org/content/27/5/722.full)
* Masurca: [new version 3.2.2: (ftp://ftp.genome.umd.edu/pub/MaSuRCA/latest/) (http://www.genome.umd.edu/masurca.html) ([supports ONT reads](http://masurca.blogspot.com.au/2017/05/please-upgrade-to-masurca-release-322.html) but...no instructions for ONT yet) here is a [presentation](https://www.youtube.com/watch?v=BEsDkjYKHCE&feature=youtu.be) and [tweet](https://twitter.com/conchoecia/status/920771453728854016) indicating possible to use ONT reads...), [quick start](http://bit.ly/2zl7G5Q)
* [SPAdes hybrid assembly](http://cab.spbu.ru/files/release3.11.1/manual.html) (more for bacteria, according to Lex Nederbragt)
* [Alpaca, hybrid assembly](https://github.com/VicugnaPacos/ALPACA/) and [paper](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3927-8)

### Evaluation:

* [Nanopolish](https://github.com/jts/nanopolish)
* [sourmash](http://sourmash.readthedocs.io/en/latest/tutorials.html) to classify and compare contig sequence similarity between assemblies, compare reads to assemblies. IDEA FROM TITUS (12/5/2017): use sourmash to cluster contigs by similarity, then do separate dotplots:
https://github.com/ctb/2017-sourmash-cluster/blob/master/cluster-big-hashsets.ipynb
* [BUSCO](chrome-extension://oemmndcbldboiebfnladdacbdfmadadm/https://gitlab.com/ezlab/busco/raw/master/BUSCO_v3_userguide.pdf) content evaluation
* [Porechop](https://github.com/rrwick/Porechop): removing adapters
* [Nanocomp](https://github.com/wdecoster/nanocomp), compare multiple ONT runs
* [Nanocorr](https://github.com/jgurtowski/nanocorr): error correction (See [presentation](http://schatzlab.cshl.edu/presentations/2015/2015.02.28.AGBT%20Nanocorr%20Assembly.pdf))

### Polishing, scaffolding and correcting:
* [Unicycler polish](https://github.com/rrwick/Unicycler/blob/master/docs/unicycler-polish.md) ([paper])(http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005595)
Polishes completed assembly with Illumina data:
https://github.com/nanoporetech/ont-assembly-polish
Nanopolish:
https://github.com/jts/nanopolish
* SSPACE-LongRead scaffolding: https://www.baseclear.com/genomics/bioinformatics/basetools/SSPACE-longread
(no instructions, but a [paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-211)!)

### Other

* [npScarf](https://www.nature.com/articles/ncomms14515)
# Results:

Method | Tool | bases | n_contigs | average | largest | N50
--- | --- | --- | --- | --- | --- | ---
ONT, Long reads | Canu | 9804264 | 540 | 18156.04 | 365191 | 40681, n = 43 
ONT, Long reads | Masurca | xxx | xxx
ONT, Long reads | Miniasm | 4917546 | 153 | 32140.82 | 233136 | 50056, n = 25
Illumina, Short reads | Megahit | 1183861293 | 1038799 | 1139.64 | 88218 | 3846, n = 77800 
Illumina, Short reads | ABySS | xxx | xxx
Hybrid | SPAdes | xxx | xxx
Hybrid | Alpaca | xxx | xxx

# Pipeline Workflow


## 1. Assess reads


### ONT adapters

[Porechop](https://github.com/rrwick/Porechop#known-adapters) has a list of known adapters found in data.

### Which long reads do not align to the short reads?


#### bwa alignments 

fasta of onp reads and fasta of Illumina NovaSeq reads:

```
module load bwa
cd /mnt/scratch/ljcohen/onp_Porecamp_killifish/fastq2/
# the following command took >5655.04 seconds (1.5 hrs)
bwa index porecamp_killifish.fasta
bwa mem -t4 -x ont2d porecamp_killifish  | samtools view -bS - | samtools sort - mapped_reads.sorted
samtools index mapped_reads.sorted

```

Map Illumina reads to assemblies, Miniasm and Canu assemblies:

```
bwa index Miniasm
bwa index Canu
bwa mem 

```

Map Illumina reads to Fhet ncbi version

```
bwa mem /mnt/home/ljcohen/reference/Fhet_ncbi

```

Map ONT assembly to F. heteroclitus reference:

```
bwa mem /mnt/home/ljcohen/reference/Fhet_ncbi

```

### [Nanocorr](https://github.com/jgurtowski/nanocorr): error correction pipeline
1. BLAST Illumina reads to all raw Oxford Nanopore reads 
2. Select non-repetitive alignments 
    * First pass scans to remove “contained” alignments 
    * Second pass uses Dynamic Programming (LIS) to select set of high-identity alignments with minimal overlaps 
3. Compute consensus of each Oxford Nanopore read (uses Pacbio’s pbdagcon)

To install Nanocorr, created and activated a py2.7 conda env:

```
cd ~/bin
git clone https://github.com/jgurtowski/nanocorr
cd nanocorr
# tried py3 env but did not work
# pbcore requires python 2.7

conda create -n nanocorr-py27env python=2.7
source activate nanocorr-py27env
pip install git+https://github.com/cython/cython
pip install numpy
pip install h5py
pip install git+https://github.com/jgurtowski/pbcore_python
pip install git+https://github.com/jgurtowski/pbdagcon_python
pip install git+https://github.com/jgurtowski/jbio
pip install git+https://github.com/jgurtowski/jptools
python setup.py install
echo export PATH=$PATH:$(pwd) >> ~/.bashrc
source ~/.bashrc
```

## 2. Prepare ONT reads:

Installed [screed](http://screed.readthedocs.io/en/latest/screed.html) and converted `.fastq` to `.fasta`:

```
pip install git+https://github.com/dib-lab/screed.git
screed db porecamp_killifish2.fastq 
python -m screed dump_fasta porecamp_killifish2.fastq_screed porecamp_killifish.fasta
```


According to [Nanocorr instructions](https://github.com/jgurtowski/nanocorr): BLAST in $PATH

```
ssh dev-intel14
source activate nanocorr-py27env
module load GNU/4.8.3
module unload python
module load BLAST+/2.3.0
qsub -cwd -V -t 1:500 -j y -o nanocorr_out nanocorr.py F_olivaceus_Illumina.fasta /mnt/scratch/ljcohen/onp_Porecamp_killifish/fastq2/F_olivaceus_ONT.fasta
```
#### Prepare Illumina reads:

Combine NovaSeq read files into one fasta file for alignment to ONT reads fasta file. We received wwo lanes of NovaSeq data from TAMU sequencing facility. These are the 4 files, about 20-25G each:
```
Killi_Fish_S19_L001_R1_001.fastq.gz
Killi_Fish_S19_L001_R2_001.fastq.gz
Killi_Fish_S19_L002_R1_001.fastq.gz
Killi_Fish_S19_L002_R2_001.fastq.gz
```
Combine reads from 2 lanes, 40G and 47G each:

```
cd /mnt/scratch/ljcohen/NovaSeq_killifish
cat Killi_Fish_S19_L001_R1_001.fastq.gz Killi_Fish_S19_L002_R1_001.fastq.gz > F_olivaceus_R1.fastq.gz
cat Killi_Fish_S19_L001_R2_001.fastq.gz 
Killi_Fish_S19_L002_R2_001.fastq.gz > F_olivaceus_R2.fastq.gz
```
Then combine R1 and R2 into a separate file for alignment, 400 G file:
```
zcat F_olivaceus_R1.fastq.gz F_olivaceus_R2.fastq.gz > F_olivaceus.fastq.gz
```
Convert combined Illumina reads into .fa file with [screed](http://screed.readthedocs.io/en/latest/screed.html):
```
screed db F_olivaceus.fastq
python -m screed dump_fasta F_olivaceus.fastq.fastq_screed F_olivaceus_Illumina.fasta
```

## 3. Assemble ONT long reads

## Miniasm

Download and install [miniasm](https://github.com/lh3/miniasm) and [minimap](https://github.com/lh3/minimap):
```
cd ~/bin
git clone https://github.com/lh3/minimap
cd minimap
make
echo export PATH=$PATH:$(pwd) >> ~/.bashrc
cd ..
git clone https://github.com/lh3/miniasm
cd miniasm
make
echo export PATH=$PATH:$(pwd) >> ~/.bashrc
source ~/.bashrc
```
# Overlap

```
minimap -Sw5 -L100 -m0 -t8 /mnt/scratch/ljcohen/onp_Porecamp_killifish/fastq2/F_olivaceus_ONT.fasta /mnt/scratch/ljcohen/onp_Porecamp_killifish/fastq2/F_olivaceus_ONT.fasta | gzip -1 > Fo.paf.gz
```
# Layout
```
miniasm -f /mnt/scratch/ljcohen/onp_Porecamp_killifish/fastq2/F_olivaceus_ONT.fasta Fo.paf.gz > Fo.gfa
```
```
awk '/^S/{print ">"$2"\n"$3}' Fo.gfa | fold > Fo.fa
assembly-stats Fo.fa
```

Minimap

https://github.com/lh3/miniasm

## Canu

For some reason I had a very difficult time using the integration with batch management scheduler. Filed a [Github Issue](https://github.com/marbl/canu/issues/566) with canu.

Command

```
canu \
-p killifish -d killifish_assembly \
genomeSize=1.1g useGrid=false \
-nanopore-raw porecamp_killifish2.fastq
```


#### On MSU hpcc:

Install canu v1.6:

```
wget https://github.com/marbl/canu/archive/v1.6.tar.gz
gunzip -dc v1.6.tar.gz | tar -xf -
cd canu-1.6/src
make -j 8
```
Add to path:
```
cd /mnt/home/ljcohen/bin/canu-1.6/Linux-amd64/bin
echo export PATH=$PATH:$(pwd) >> ~/.bashrc
source ~/.bashrc
```
(had previously installed `canu` with `brew install canu` and had version 1.5, so uninstalled with `brew uninstall canu` then `source ~/.bashrc` to put latest version in `$PATH`)


#### On UCD farm:

Ran test command (this worked):
```
srun -A millermrgrp -p bigmemh -t 24:00:00 --mem=100000 --pty bash
source ~/.bashrc
cd /home/ljcohen/ONP_killifish/test
#curl -L -o mix.tar.gz http://gembox.cbcb.umd.edu/mhap/raw/ecoliP6Oxford.tar.gz
#tar xvzf mix.tar.gz
module load java/1.8
canu  -p ecoli -d ecoli-mix useGrid=false genomeSize=4.8m  -pacbio-raw pacbio.part?.fastq.gz  -nanopore-raw oxford.fasta.gz
```


With killifish data, ran command:
```
# start interactive job
srun -A millermrgrp -p bigmemh -t 24:00:00 --mem=100000 --pty bash
# put canu in $PATH
source ~/.bashrc
# canu requires specific version of java
module load java/1.8
cd /home/ljcohen/ONP_killifish
canu -p F_olivaceus -d F_olivaceus_assembly useGrid=false ovsMethod=sequential genomeSize=1.1g -nanopore-raw porecamp_killifish2.fastq
```

Got a message from farm administrator, taking up too much I/O.

```
#!/bin/bash -l
#SBATCH -D /home/ljcohen/ONP_killifish/slurm_out
#SBATCH -J canu_F_olivaceus
#SBATCH -A millermrgrp
#SBATCH -p bigmeml
#SBATCH -t 96:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mem=100000

source ~/.bashrc
module load java/1.8

PROJECTDIR=/home/ljcohen/ONP_killifish

mkdir /scratch/$SLURM_JOBID
cd /scratch/$SLURM_JOBID
cp $PROJECTDIR/porecamp_killifish2.fastq .

canu -p F_olivaceus -d F_olivaceus_assembly useGrid=false ovsMethod=sequential maxMemory=140 genomeSize=1.1g -nanopore-raw porecamp_killifish2.fastq
cp -r /scratch/$SLURM_JOBID/F_olivaceus_assembly/ $PROJECTDIR/F_olivaceus_assembly
rm -rf /scratch/$SLURM_JOBID*
```

Job killed, not enough memory. Increased to --mem=150000

Same thing.

Increased to --mem=300000

Same thing. Gave up.

NOTE to try: IF that doesn't work...next, try:
minReadLength=100 minOverlapLength=50

Tried 500GB on MSU. Kept getting this error message:

```
ERROR: Input read file 'porecamp_killifish2.fastq' not found.
ERROR:  Invalid command line option 'porecamp_killifish2.fastq'.  Did you forget quotes around options with spaces?
```

Ended up putting full absolute path, which worked:

```
#!/bin/bash
#PBS -l walltime=48:00:00,nodes=1:ppn=32
#PBS -l mem=500GB
#PBS -j oe
#PBS -A ged
#PBS -M ljcohen@msu.edu
#PBS -m ae
set -e
set -x

source ~/.bashrc
canu --version
canu -p F_olivaceus -d F_olivaceus_assembly \
useGrid=false ovsMethod=sequential genomeSize=1.1g \
-nanopore-raw /mnt/scratch/ljcohen/onp_Porecamp_killifish/fastq2/porecamp_killifish2.fastq
```

Ran for 42 of allocated 48 hrs and realized this wasn't going to work. Increased to 750 GB and 168 hrs (1 week):

```
#!/bin/bash
#PBS -l walltime=168:00:00,nodes=1:ppn=32
#PBS -l mem=750GB
#PBS -j oe
#PBS -A ged
#PBS -M ljcohen@msu.edu
#PBS -m ae
set -e
set -x

source ~/.bashrc
canu --version
canu -p F_olivaceus -d /tmp/F_olivaceus_assembly \
useGrid=false ovsMethod=sequential genomeSize=1.1g \
-nanopore-raw /mnt/scratch/ljcohen/onp_Porecamp_killifish/fastq2/porecamp_killifish2.fastq

cp -r /tmp/F_olivaceus_assembly/ /mnt/scratch/ljcohen/onp_Porecamp_killifish/
rm -rf /tmp/F_olivaceus_assembly/
```
This worked!!!

```
[ljcohen@dev-intel14 F_olivaceus_assembly]$ assembly-stats F_olivaceus.contigs.fasta
stats for F_olivaceus.contigs.fasta
sum = 9804264, n = 540, ave = 18156.04, largest = 365191
N50 = 40681, n = 43
N60 = 28110, n = 72
N70 = 20781, n = 114
N80 = 14262, n = 169
N90 = 7269, n = 266
N100 = 1018, n = 540
N_count = 0
Gaps = 0
```

Resources used on MSU hpcc:
```

PBS Job Id: 49124933.mgr-04.i
Job Name:   canu_F_olivaceus.qsub
Exec host:  ifi-002/0-31
Execution terminated
Exit_status=0
resources_used.cput=746:12:08
resources_used.vmem=139095308kb
resources_used.walltime=41:52:06
resources_used.mem=134300784kb
resources_used.energy_used=0
```

## 4. Assembly with Illumina short reads

### ABySS

```
module load GNU/4.4.5 OpenMPI/1.4.3 Boost/1.47.0 ABySS/1.9.0
abyss-pe k=96 name=F_olivaceus_Illumina in="F_olivaceus_R1.fastq.gz F_olivaceus_R2.fastq.gz" contigs
```

### Megahit

Installed [megahit](https://github.com/voutcn/megahit) on MSU, from dibsi [tutorial](http://angus.readthedocs.io/en/2017/genome-assembly.html):

```
cd ~/bin
git clone https://github.com/voutcn/megahit.git
cd megahit
make
echo export PATH=$PATH:$(pwd) >> ~/.bashrc
source ~/.bashrc
```

Ran this without diginorm or trimming, just to see. Took ~21 hrs. 

## 5. Hybrid assembly with both Illumina and ONT

### SSPACE-LongRead

Installation:

Go to [BaseClear website](https://www.baseclear.com/services/bioinformatics/basetools/sspace-longread/) and signup for an academic license. This will take you to a [page](https://www.baseclear.com/services/bioinformatics/basetools/sspace-standard/thank-you/?business=academic-581227) for downloading the files.

```
wget https://github.com/ljcohen/hybrid_genome_assembly/blob/master/SSPACE/SSPACE-longread-v.-1-1.tar.gz
tar -xvzf SSPACE-longread-v.-1-1.tar.gz 
cd SSPACE-LongRead_v1-1
perl ./SSPACE-LongRead.pl #displays usage information
```

Command:
```
perl SSPACE-LongRead.pl -c <contig-sequences> -p <pacbio-reads>
```
Not sure what format `-p` reads have to be in. The "help" specified: "File containing PacBio CLR sequences to be used scaffolding". Tried a fastq:
```
ssh dev-intel14
cd /mnt/scratch/ljcohen/NovaSeq_killifish/F_olivaceus/
perl ~/bin/SSPACE-LongRead_v1-1/SSPACE-LongRead.pl -c final.contigs.fa -p ../../onp_Porecamp_killifish/fastq2/porecamp_killifish_ONT_reads.fastq
```
Did not throw an error!

## MaSurCa

Install:
```
# MaSuRCA 3.2.6 installation on MSU hpcc
ssh dev-intel14
cd bin
wget https://github.com/alekseyzimin/masurca/releases/download/3.2.6/MaSuRCA-3.2.6.tar.gz
tar -xvzf MaSuRCA-3.2.6.tar.gz
cd MaSuRCA-3.2.6
module load GNU/4.8.3
BOOST_ROOT=install ./install.sh
```
Config file:
```
# example configuration file 

# DATA is specified as type {PE,JUMP,OTHER,PACBIO} and 5 fields:
# 1)two_letter_prefix 2)mean 3)stdev 4)fastq(.gz)_fwd_reads
# 5)fastq(.gz)_rev_reads. The PE reads are always assumed to be
# innies, i.e. --->.<---, and JUMP are assumed to be outties
# <---.--->. If there are any jump libraries that are innies, such as
# longjump, specify them as JUMP and specify NEGATIVE mean. Reverse reads
# are optional for PE libraries and mandatory for JUMP libraries. Any
# OTHER sequence data (454, Sanger, Ion torrent, etc) must be first
# converted into Celera Assembler compatible .frg files (see
# http://wgs-assembler.sourceforge.com)
DATA
PE= pe 150 20  /mnt/research/ged/lisa/killifish_onp/Folivaceus/subset_Folivaceus/subset_F_olivaceus_R1.fastq /mnt/research/ged/lisa/killifish_onp/Folivaceus/subset_Folivaceus/subset_F_olivaceus_R2.fastq
#JUMP= sh 3600 200  /FULL_PATH/short_1.fastq  /FULL_PATH/short_2.fastq
#pacbio reads must be in a single fasta file! make sure you provide absolute path
NANOPORE=/mnt/research/ged/lisa/killifish_onp/Folivaceus/subset_Folivaceus/subset_F_olivaceus_ONT_1m.fasta
#OTHER=/FULL_PATH/file.frg
END

PARAMETERS
#set this to 1 if your Illumina jumping library reads are shorter than 100bp
EXTEND_JUMP_READS=0
#this is k-mer size for deBruijn graph values between 25 and 127 are supported, auto will compute the optimal size based on the read data and GC content
GRAPH_KMER_SIZE = auto
#set this to 1 for all Illumina-only assemblies
#set this to 1 if you have less than 20x long reads (454, Sanger, Pacbio) and less than 50x CLONE coverage by Illumina, Sanger or 454 mate pairs
#otherwise keep at 0
USE_LINKING_MATES = 1
#specifies whether to run mega-reads correction on the grid
USE_GRID=0
#specifies queue to use when running on the grid MANDATORY
GRID_QUEUE=all.q
#batch size in the amount of long read sequence for each batch on the grid
GRID_BATCH_SIZE=300000000
#coverage by the longest Long reads to use
LHE_COVERAGE=30
#this parameter is useful if you have too many Illumina jumping library mates. Typically set it to 60 for bacteria and 300 for the other organisms 
LIMIT_JUMP_COVERAGE = 300
#these are the additional parameters to Celera Assembler.  do not worry about performance, number or processors or batch sizes -- these are computed automatically. 
#set cgwErrorRate=0.25 for bacteria and 0.1<=cgwErrorRate<=0.15 for other organisms.
CA_PARAMETERS =  cgwErrorRate=0.15
#minimum count k-mers used in error correction 1 means all k-mers are used.  one can increase to 2 if Illumina coverage >100
KMER_COUNT_THRESHOLD = 1
#whether to attempt to close gaps in scaffolds with Illumina data
CLOSE_GAPS=1
#auto-detected number of cpus to use
NUM_THREADS = 16
#this is mandatory jellyfish hash size -- a safe value is estimated_genome_size*estimated_coverage
JF_SIZE = 200000000
#set this to 1 to use SOAPdenovo contigging/scaffolding module.  Assembly will be worse but will run faster. Useful for very large (>5Gbp) genomes from Illumina-only data
SOAP_ASSEMBLY=0
END
```
Run (on msu hpcc):
```
ssh dev-intel14
module load GNU/4.8.3
mkdir subset_run/
cd subset_run/
cp ../config.txt .
../bin/masurca config.txt
```
This will give you output like this:
```
[ljcohen@dev-intel14 subset_run]$ ../bin/masurca config.txt
Verifying PATHS...
jellyfish OK
runCA OK
createSuperReadsForDirectory.perl OK
creating script file for the actions...done.
execute assemble.sh to run assembly
```
Then:
```
./assemble.sh
```


Output:
```
[ljcohen@dev-intel14 subset_run]$ ./assemble.sh
[Fri Jul 20 20:01:14 EDT 2018] Processing pe library reads
[Fri Jul 20 20:01:24 EDT 2018] Average PE read length 150
[Fri Jul 20 20:01:24 EDT 2018] Using kmer size of 99 for the graph
[Fri Jul 20 20:01:25 EDT 2018] MIN_Q_CHAR: 33
[Fri Jul 20 20:01:25 EDT 2018] Creating mer database for Quorum
[Fri Jul 20 20:01:49 EDT 2018] Error correct PE.
[Fri Jul 20 20:02:10 EDT 2018] Estimating genome size.
[Fri Jul 20 20:02:13 EDT 2018] Estimated genome size: 5165648
[Fri Jul 20 20:02:13 EDT 2018] Creating k-unitigs with k=99
[Fri Jul 20 20:02:17 EDT 2018] Computing super reads from PE 
[Fri Jul 20 20:02:22 EDT 2018] Using linking mates
Using CABOG from is /mnt/home/ljcohen/bin/MaSuRCA-3.2.6/bin/../CA8/Linux-amd64/bin
Running mega-reads correction/assembly
Using mer size 15 for mapping, B=15, d=0.02
Estimated Genome Size 5165648
Estimated Ploidy 1
Using 16 threads
Output prefix mr.41.15.15.0.02
Using 30x of the longest ONT reads
Reducing super-read k-mer size
Mega-reads pass 1
Running locally in 1 batch
compute_psa 23286 4730687
Mega-reads pass 2
Running locally in 1 batch
compute_psa 8897 4870528
Refining alignments
Joining
Generating assembly input files
Coverage of the mega-reads less than 5 -- using the super reads as well
Coverage threshold for splitting unitigs is 15 minimum ovl 98
Running assembly
Recomputing A-stat
recomputing A-stat for super-reads
Mega-reads initial assembly complete.
[Fri Jul 20 20:05:22 EDT 2018] Gap closing.
[Fri Jul 20 20:05:42 EDT 2018] Gap close success.
[Fri Jul 20 20:05:43 EDT 2018] Assembly complete, final scaffold sequences are in CA.mr.41.15.15.0.02/final.genome.scf.fasta
[Fri Jul 20 20:05:43 EDT 2018] All done
```

stats:
```
source activate dammit_new
(dammit_new) [ljcohen@dev-intel14 subset_run]$ assembly-stats guillaumeKUnitigsAtLeast32bases_all.fasta
stats for guillaumeKUnitigsAtLeast32bases_all.fasta
sum = 9701456, n = 73629, ave = 131.76, largest = 2968
N50 = 134, n = 30256
N60 = 124, n = 37803
N70 = 115, n = 45954
N80 = 108, n = 54688
N90 = 102, n = 63946
N100 = 99, n = 73629
N_count = 0
Gaps = 0
```
## 6. Evaluation of the Assemblies

### Quast

Evaluated with Quast (from [dibsi tutorial](http://angus.readthedocs.io/en/2017/genome-assembly.html)):

```
cd ~/bin
git clone https://github.com/ablab/quast.git -b release_4.5
/mnt/scratch/ljcohen/onp_Porecamp_killifish/F_olivaceus_assembly/Foli_canu_quast_report/
export PYTHONPATH=$(pwd)/quast/libs/
quast/quast.py /mnt/scratch/ljcohen/onp_Porecamp_killifish/F_olivaceus_assembly/F_olivaceus.contigs.fasta -o /mnt/scratch/ljcohen/onp_Porecamp_killifish/F_olivaceus_assembly/Foli_canu_quast_report
quast/quast.py /mnt/scratch/ljcohen/onp_Porecamp_killifish/fastq2/miniasm/Fo.fa -o /mnt/scratch/ljcohen/onp_Porecamp_killifish/fastq2/miniasm/Foli_miniasm_quast_report
```

Then downloaded locally to visualize quast report output files:
```
scp -r ljcohen@rsync.hpcc.msu.edu:/mnt/scratch/ljcohen/onp_Porecamp_killifish/fastq2/miniasm/Foli_miniasm_quast_report/ .
scp -r ljcohen@rsync.hpcc.msu.edu:/mnt/scratch/ljcohen/onp_Porecamp_killifish/F_olivaceus_assembly/Foli_canu_quast_report/ .
scp -r ljcohen@rsync.hpcc.msu.edu:/mnt/scratch/ljcohen/NovaSeq_killifish/F_olivaceus/Foli_Illumina_megahit_report .
```

## Compare reads to assemblies with [sourmash](http://sourmash.readthedocs.io/en/latest/tutorials.html)

[Installation](https://gist.github.com/ljcohen/aa2c327b7f57f60edd2f9b7824268fbb) on msu hpcc, setup [conda env](https://conda.io/docs/user-guide/tasks/manage-environments.html).
```
module load GNU/4.8.3  anaconda parallel
module unload python
source activate sourmash_nov2017
```

Compute signatures of raw ONT reads:
```
sourmash compute --scaled 1000 ../porecamp_killifish2.fastq -o Foli_ONT.sig -k 31
```
Compute signatures of raw Illumina NovaSeq reads:
```
sourmash compute --scaled 10000 ../F_olivaceus_R1.fastq.gz -o Foli_Illumina_R1.sig -k 31
```

## Congener genome alignment

We have a closely-related reference genome from [*Fundulus heteroclitus*](https://www.ncbi.nlm.nih.gov/genome/?term=txid8078[Organism:exp]). Which parts of the genome do not align? Are these unique features or sequencing error?

1. [Downloaded](https://www.ncbi.nlm.nih.gov/assembly?LinkName=genome_assembly&from_uid=743) genome from NCBI (293M compressed, 989M uncompressed).

```
cd ~/reference
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/826/765/GCF_000826765.1_Fundulus_heteroclitus-3.0.2/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fna.gz
gunzip GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fna.gz
```
Visualize, sanity check:
```
head /mnt/home/ljcohen/reference/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fna
```

2. [Installed](https://github.com/mummer4/mummer/blob/master/INSTALL.md) mummer [v4.0.0.beta2](https://github.com/mummer4/mummer/releases/download/v4.0.0beta2/mummer-4.0.0beta2.tar.gz), for nucmer alignment of two genomes.

```
wget https://github.com/mummer4/mummer/releases/download/v4.0.0beta2/mummer-4.0.0beta2.tar.gz
tar -xzvf mummer-4.0.0beta2.tar.gz
cd mummer-4.0.0beta2
./configure --prefix=/mnt/home/ljcohen/bin/mumm-4.0.0beta2
make
make install
echo export PATH=$PATH:$(pwd) >> ~/.bashrc
source ~/.bashrc
```
3. All by all comparison with ONT reads:

Align:

```
mkdir /mnt/home/ljcohen/onp_killifish/mummer
cd /mnt/home/ljcohen/onp_killifish/mummer
nucmer --maxmatch --threads=7 \
-p F_oli_ONT_Fhet \
/mnt/scratch/ljcohen/onp_Porecamp_killifish/fastq2/F_olivaceus_ONT.fasta \
/mnt/home/ljcohen/reference/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fna
```

Reverse:
```
cd /mnt/home/ljcohen/onp_killifish/mummer
nucmer --maxmatch \
-p F_oli_Fhet_ONT /mnt/home/ljcohen/reference/GCF_000826765.1_Fundulus_heteroclitus-3.0.2_genomic.fna /mnt/scratch/ljcohen/onp_Porecamp_killifish/fastq2/F_olivaceus_ONT.fasta
```
Plot:

```
mummerplot --fat --filter --png -l --large -p F_olivaceus F_oli_ONT_Fhet.delta
```

Before running mummerplot, edit the ~/bin/mummer-3.9.4alpha/mummerplot script to comment out the three lines that have the word "mouse", because [reasons](https://sourceforge.net/p/mummer/mailman/message/34939032/). (Search for `clipboardformat`.)


Had to install different version of gnuplot. Contributed to this [issue](https://github.com/mummer4/mummer/issues/36), where it seems this is a common occurance.
```
cd ~/bin
wget wget https://downloads.sourceforge.net/project/gnuplot/gnuplot/5.2.0/gnuplot-5.2.0.tar.gz
tar -xvzf gnuplot-5.2.0.tar.gz
cd gnuplot-5.2.0
./configure --prefix=/mnt/home/ljcohen/bin/
make install
cd src
./gnuplot --version # confirms the version installed
echo export PATH=$(pwd):$PATH >> ~/.bashrc
source ~/.bashrc
```

Plots of F. heteroclitus vs. canu:

Based on this, can see portions of F. heteroclitus genome not covered by ONT canu assembly.


Plots of Miniasm and Canu:

```
module load GNU/4.8.3
source ~/.bashrc
nucmer --maxmatch -p F_oli_Miniasm_Canu /mnt/scratch/ljcohen/onp_Porecamp_killifish/fastq2/miniasm/Fo.fa /mnt/scratch/ljcohen/onp_Porecamp_killifish/F_olivaceus_assembly/F_olivaceus.contigs.fasta

mummerplot --fat --filter --png -l --large -p F_olivaceus_Miniasm_Canu F_oli_Miniasm_Canu.delta
```
# Remaining Questions

1. Are there better ways to easily assess reads, align, visualize, make decisions about data?

2. 
# References

* [Tutorial](http://angus.readthedocs.io/en/2017/analyzing_nanopore_data.html) on assessing ONT reads, assembling, and evaluation from [DIBSI](http://ivory.idyll.org/dibsi/). 
* Jessica Mizzi, [Thule elk draft genome](https://github.com/jessicamizzi/tule-elk) and [blogpost](http://ivory.idyll.org/blog/2016-tule-elk-draft.html)
* Taylor Reiter, [olive genome assessment](https://github.com/taylorreiter/olive_genome)
* Connor Tiffany, [E. coli genome assembly](https://github.com/recursive-deletion/EHEC_genome)
* [PacBio large genome assembly references](https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Large-Genome-Assembly-with-PacBio-Long-Reads)
* Motivation paper: [Hybrid assembly with long and short reads improves discovery of gene family expansions](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3927-8)
* [Wiki: Hybrid genome assembly](https://en.wikipedia.org/wiki/Hybrid_genome_assembly)
* Methods comparison paper: [De novo yeast genome assemblies from MinION, PacBio and MiSeq platforms](https://www.nature.com/articles/s41598-017-03996-z)
* [A "forward genomics" approach links genotype to phenotype using independent phenotypic losses among related species.](https://www.ncbi.nlm.nih.gov/pubmed/23022484#)
* [MAKER](http://www.yandell-lab.org/software/maker.html) and [MAKER-P](http://www.yandell-lab.org/software/maker-p.html), [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4286374/) and [Annotation Edit Distance (AED) evaluation metric](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-67)
* Alignment of whole genomes https://mummer4.github.io/publications/MUMmer.pdf
* MUMMer tutorial: https://mummer4.github.io/tutorial/tutorial.html
* What is finished and why does it matter? http://genome.cshlp.org/content/12/5/669.full.html
* De novo genome assembly: what every biologist should know http://www.nature.com/nmeth/journal/v9/n4/full/nmeth.1935.html
* A field guide to whole-genome sequencing, assembly and annotation http://onlinelibrary.wiley.com/doi/10.1111/eva.12178/full
* Titus Brown, Family Relations: http://family.caltech.edu/tutorial/
* Titus Brown, Family Jewels: http://family.caltech.edu/
