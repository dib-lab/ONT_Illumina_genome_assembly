## SSPACE install files 

Install files for SSPACE, including SSPACE-LongRead for scaffolding genomes with long read sequence data.
https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-211

Go to the BaseClear website and enter your information for an academic license:
https://www.baseclear.com/services/bioinformatics/basetools/sspace-longread/

License information:
https://www.baseclear.com/wp-content/uploads/BaseTools-License-agreement-v.-2014-03-31.pdf

Installation:
```
wget https://github.com/ljcohen/hybrid_genome_assembly/blob/master/SSPACE/SSPACE-longread-v.-1-1.tar.gz
tar -xvzf SSPACE-longread-v.-1-1.tar.gz 
cd SSPACE-LongRead_v1-1
perl ./SSPACE-LongRead.pl #displays usage information
```
Requires perl.

Command:
```
perl SSPACE-LongRead.pl -c <contig-sequences> -p <pacbio-reads>
```

SSPACE (standard):

Boetzer M, Henkel CV, Jansen HJ, Butler D, Pirovano W (2011), Scaffolding pre-assembled contigs
using SSPACE, Bioinformatics, 27(4):578-9.

SSPACE-LongRead:
Boetzer M, Pirovano W (2014), SSPACE-LongRead: Finishing bacterial draft genomes using long read
sequence information, submitted for publication.

GapFiller:
Boetzer M, Pirovano W (2012), Toward almost closed genomes with GapFiller, Genome Biology, 13(6).

Usage information:
```
Usage SSPACE-LongRead scaffolder version 1-1

perl SSPACE-LongRead.pl -c <contig-sequences> -p <pacbio-reads>

General options:
-c  Fasta file containing contig sequences used for scaffolding (REQUIRED)
-p  File containing PacBio CLR sequences to be used scaffolding (REQUIRED)
-b  Output folder name where the results are stored (optional, default -b 'PacBio_scaffolder_results')

Alignment options:
-a  Minimum alignment length to allow a contig to be included for scaffolding (default -a 0, optional)
-i  Minimum identity of the alignment of the PacBio reads to the contig sequences. Alignment below this value will be filtered out (default -i 70, optional)
-t  The number of threads to run BLASR with
-g  Minimmum gap between two contigs

Scaffolding options:
-l  Minimum number of links (PacBio reads) to allow contig-pairs for scaffolding (default -k 3, optional)
-r  Maximum link ratio between two best contig pairs *higher values lead to least accurate scaffolding* (default -r 0.3, optional)
-o  Minimum overlap length to merge two contigs (default -o 10, optional)

Other options:
-k  Store inner-scaffold sequences in a file. These are the long-read sequences spanning over a contig-link (default no output, set '-k 1' to store inner-scaffold sequences. If set, a folder is generated named 'inner-scaffold-sequences'
-s  Skip the alignment step and use a previous alignment file. Note that the results of a previous run will be overwritten. Set '-s 1' to skip the alignment.
-h  Prints this help message

ERROR: Please insert a file with contig sequences. You've inserted '' which either does not exist or is not filled in
```
