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
