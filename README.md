# README for TISCA

2021-3-26 Kazuya Ichihara

#### Obtain genome annotation files

* Download human genome and transcriptome annotation files from GENCODE (Harrow et al., 2012).

```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.transcripts.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/GRCh38.primary_assembly.genome.fa.gz
```

* Get gtfToGenePred command line from UCSC Genome Browser, and use the tool to convert GTF file to genePred format. 

```bash
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
gtfToGenePred gencode.v31.annotation.gtf gencode.v31.annotation.genePred.txt
```

#### Execute RibORF

* Download RibORF package from https://github.com/zhejilab/RibORF/. 

* Generate candidate ORFs in genePred format.

```bash
mkdir RibORF_output
perl ./RibORF/RibORF.1.0/ORFannotate.pl \
-g GRCh38.primary_assembly.genome.fa -t gencode.v31.annotation.genePred.txt \
-o ./RibORF_output -s ATG/CTG/TTG/GTG/AGG/ACG/AAG/ATC/ATA/ATT
```

* Examine the read distribution

```bash
samtools view -@ 8 -h ./data/Ribo.bam > Ribo.sam
perl ./RibORF/RibORF.1.0/readDist.pl -f ./Ribo.sam -g gencode.v31.annotation.genePred.txt -o RibORF_output
```

* Offset correction

```bash
perl ./RibORF/RibORF.1.0/offsetCorrect.pl -r Ribo.sam -p ./data/offset.correction.parameters.txt -o Ribo_corrected.sam
perl ./RibORF/RibORF.1.0/readDist.pl -f Ribo_corrected.sam -g gencode.v31.annotation.genePred.txt -o RibORF_output
```

* Calculate pred.pvalue

```bash
perl ./RibORF/RibORF.1.0/ribORF.pl -f Ribo_corrected.sam -c ./RibORF_output/candidateORF.genepred.txt -o RibORF_output
```

#### Prepare annotation files for TISCA

```bash
python CreateAnnotation.py
```

#### Create read aggregation plots

```
python readdist.py ./data/Ribo.bam ./annotation/Ribo
python readdist.py ./data/GTI.bam ./annotation/GTI
python readdist.py ./data/Sel.bam ./annotation/Sel
python readend.py ./data/Sel.bam ./annotation/Sel
```

#### Execute TISCA

```
python TISCA.py
```
