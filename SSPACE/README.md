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
