THIS REPOSITORY IS A WORK_IN_PROGRESS.

Hybrid genome assembly of the killifish, *Fundulus olivaceus* with ONT and Illumina NovaSeq
===

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

[Methods: Workflow, commands used and notes](https://github.com/ljcohen/hybrid_genome_assembly/blob/master/workflow.md)

Johnson, L. (2018, July 20). Hybrid genome assembly of the freshwater killifish, Fundulus olivaceus with ONT and Illumina NovaSeq. https://doi.org/10.17605/osf.io/zjv86
