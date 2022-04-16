# Marker Tech

[![platform](https://img.shields.io/badge/platform-Website-blue)](http://zeasystemsbio.hzau.edu.cn:8080/) [![copyright](https://img.shields.io/badge/copyright-LiLab-%23298850)](https://faculty.hzau.edu.cn/lilin12/zh_CN/index.htm) [![doc](https://img.shields.io/badge/doc-Primer-red)](https://github.com/9-34platform/MarkerTech/blob/master/PrimerDesigner/引物设计开发文档.md) [![doc](https://img.shields.io/static/v1?label=doc&message=iBP-Seq&color=red)](https://github.com/9-34platform/MarkerTech/blob/master/iBP-Seq/iBP-Seq%20开发文档.md) 

A cost-saving genotyping pattern based on the next-generation sequencing (NGS) technology. We also build an online analysis system which provides several services for using the technology in an easier way. However, in order to meet diverse and highly customized needs, the main source codes are provided for users to modify as needed. The services we provide include primer design service and iBP-Seq (improved Bulked-PCR Seq) analysis service which can do both genotyping and methylation typing. The primer design service generates special primer sequences flanking the maker for each template sequence. The NGS data for the sequence analysis service should be generate by using the primer design service. 

## Preparation

### Primer design

1. **Environment**
   - python>=3.6.0
     - primer3-py>=0.6.1
3. **Data**
   - a template sequence file in fasta format containing template sequences (header format: \>Name|301-304; Sequence length>=601)

### BP-Seq & iBP-Seq analysis

Due to the limitation of some modules, such as pysam, we recommend running these main codes on a Linux operating system when you call them from the command line. 

Both the genotyping service for iBP-Seq and methylation typing service require these configurations: 

1. **Environment**
   - python>=3.6.0
     - openpyxl>=3.0.7
     - XlsxWriter>=1.4.0
     - seaborn
     - numpy
     - pandas
     - matplotlib==3.4.3
     - sklearn
     - pysam==0.16.0.1
   - fastp>=0.21.0
   - bwa>=0.7.17
   - samtools>=1.12
3. **Data**
   - a pair of paired-end sequencing data in fastq format
   - a barcode file which order is consistent with the actual order of use (separator='\n')
   - a mutation table with header include four columns (CHROM, POS, REF, ALT)
   - **optional**: a reference sequence file in fasta format 

### mBP-Seq analysis

1. **Environment**
   - python>=3.6.0
     - openpyxl>=3.0.7
     - XlsxWriter>=1.4.0
     - seaborn
     - numpy
     - pandas
     - matplotlib==3.4.3
     - sklearn
     - pysam==0.16.0.1
   - fastp>=0.21.0
   - bsmark>=0.23.1
   - samtools>=1.12

### Ready to work

Before starting the analysis, if the sequencing file is in .gz compressed format, it also needs to be decompressed. 

~~~shell
gzip -dk <filename>.fq.gz
~~~

And make sure you have the following empty directories in the working directory. If not, create them first.

**BP** : 

~~~shell
mkdir qc fq bam csv xlsx png
~~~

**iBP** : 

~~~shell
mkdir qc fq umi bam csv xlsx png
~~~

**mBP**: 

~~~shell
mkdir qc fq umi sam bam csv xlsx png
~~~

## Usage

### Primer design

The input parameters of PrimerDesigner.py are: strategy number (BP-Seq: 2; iBP-Seq: 3), input template file, file name (including .xlsx suffix), shortest primer length, optimal primer length, longest primer length, shortest amplicon length, optimal amplicon length, longest amplicon length, minimum Tm of primer, optimal Tm of primer, maximum Tm of primer. 

An example of the command line is as follows: 

~~~shell
python3 PrimerDesigner.py 3 testSeq.fa ./filename.xlsx 18 20 25 280 286 320 50 60 65
~~~

### BP-Seq analysis

> A pipeline for the base genotyping strategy. 

#### Quality control

Filter out low quality (lower than 20) and sequence length less than 154bp (20bp EndR + 2\*8bp barcode + 60bp YFP + 18bp bridge + 2\*20bp primer = 154bp) reads. However, due to the read length of next-generation sequencing, we removed the 40bp length of the primers to relax the constraints. 

~~~shell
fastp -i <filename_1>.fq -I <filename_2>.fq -o qc/clean_1.fq -O qc/clean_2.fq -q 20 --length_required 114 -j qc/<filename>.json -h qc/<filename>.html
~~~

#### Data division

Divide the mixed sequencing data to several files by using the barcode file and build-in bridge sequence (selected by the first parameter). The number of barcodes in the file decides the number of divided samples and the order of barcodes decides the order of each sample. Pair-end sequencing genome sequences are stored as pairs of fastq files in the fq directory. 

~~~shell
python3 SeqPurifier.py BP barcode/<barcode>.txt qc/clean_1.fq qc/clean_2.fq fq
~~~

#### Reference alignment

Align sequence to reference sequence or genome. 

~~~shell
for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++)); do bwa mem -t 4 -M <Reference>.fa fq/${f}_1.fq fq/${f}_2.fq | samtools view -S -b - | samtools sort -o bam/${f}.sorted.bam -; done
~~~

Build index for sorted bam file. 

~~~shell
for b in $(ls bam/*); do samtools index ${b}; done
~~~

#### Marker frequency

Calculate each mutated marker frequency for each sample. 

~~~shell
python3.7 AltFraction.py BP barcode/<barcode>.txt mutation/<mut>.csv bam csv
~~~

#### KDE Cluster

Identify different genotypes or methylation types by the minimum point of kernel density estimation (KDE) from the sample population. 

~~~shell
for t in $(ls csv/*); do python3.7 Minima.py ${t} xlsx png; done
~~~

### iBP-Seq analysis

> A pipeline for an improved solution based on BP-Seq.

#### Quality control

Filter out low quality (lower than 20) and sequence length less than 72bp (2*20bp primer + 6bp UMI + 18Bbp bridge + 8bp barcode = 72bp) reads. 

~~~shell
fastp -i <filename_1>.fq -I <filename_2>.fq -o qc/clean_1.fq -O qc/clean_2.fq -q 20 --length_required 72 -j qc/<filename>.json -h qc/<filename>.html
~~~

#### Data division

Divide the mixed sequencing data to several files by using the barcode file and build-in bridge sequence (selected by the first parameter). The number of barcodes in the file decides the number of divided samples and the order of barcodes decides the order of each sample. Pair-end sequencing genome sequences are stored as pairs of fastq files in the fq directory and UMI sequences are stored as single fastq files in the umi directory. 

~~~shell
python3 SeqPurifier.py iBP barcode/<barcode>.txt qc/clean_1.fq qc/clean_2.fq fq umi
~~~

#### Reference alignment

Align sequence to reference sequence or genome. 

~~~shell
for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++)); do bwa mem -t 4 -M <Reference>.fa fq/${f}_1.fq fq/${f}_2.fq | samtools view -S -b - | samtools sort -o bam/${f}.sorted.bam -; done
~~~

Build index for sorted bam file. 

~~~shell
for b in $(ls bam/*); do samtools index ${b}; done
~~~

#### Marker frequency

Calculate each mutated marker frequency for each sample. 

~~~shell
python3.7 AltFraction.py iBP barcode/<barcode>.txt mutation/<mut>.csv bam umi csv
~~~

#### KDE Cluster

Identify different genotypes or methylation types by the minimum point of kernel density estimation (KDE) from the sample population. 

~~~shell
for t in $(ls csv/*); do python3.7 Minima.py ${t} xlsx png; done
~~~

### mBP-Seq analysis

> This pipeline is an application of iBP-Seq method in processing methylation data. 

#### Quality control

~~~shell
fastp -i <filename_1>.fq -I <filename_2>.fq -o qc/clean_1.fq -O qc/clean_2.fq -q 20 --length_required 72 -j qc/<filename>.json -h qc/<filename>.html
~~~

#### Data division

~~~shell
python3 SeqPurifier.py mBP barcode/<barcode>.txt qc/clean_1.fq qc/clean_2.fq fq umi
~~~

#### Reference alignment

##### Reference split

~~~shell
python3 MethReference.py barcode/<barcode>.txt ref/<ref>.fa
# Build index
for i in $(ls -d ref/*);do bismark_genome_preparation ${i}/ ;done
~~~

##### Alignment

~~~shell
for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++)); do lineage=$(sed -n $[${f}+1]p barcode/<barcode>.txt | sed "s/.*,//" | sed "s/\r//" | sed "s/\//_/g");bismark -N 0 -L 20 --bowtie2 --non_bs_mm /<abs_path>/ref/${lineage}/ --sam -o sam/ --q -1 fq/${f}_1.fq -2 fq/${f}_2.fq; done
~~~

#### SAM processing

~~~shell
# Text processing
for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++)); 
do
    grep "@" sam/${f}_1_bismark_bt2_pe.sam > sam/${f}-header;
    grep "XA:Z:0" sam/${f}_1_bismark_bt2_pe.sam > sam/${f}-temp1;
    grep "XB:Z:0" sam/${f}_1_bismark_bt2_pe.sam > sam/${f}-temp2;
    awk 'NR==FNR {f1[$1]=$0;next} {idx=$1;if(idx in f1)print$0}' sam/${f}-temp1 sam/${f}-temp2 > sam/${f}-temp4;
    awk 'NR==FNR {f1[$1]=$0;next} {idx=$1;if(idx in f1)print$0}' sam/${f}-temp2 sam/${f}-temp1 > sam/${f}-temp3;
    cat sam/${f}-header sam/${f}-temp3 sam/${f}-temp4 > sam/${f}_1_bismark_bt2_pe_nonmismatch.sam; 
done

# SAM to BAM
for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++));
do
    samtools sort -n -O bam -o bam/${f}_1_bismark_bt2_pe_nonmismatch.bam sam/${f}_1_bismark_bt2_pe_nonmismatch.sam;
    lineage=$(sed -n $[${f}+1]p barcode/barcode.txt | sed "s/.*,//" | sed "s/\r//" | sed "s/\//_/g");
    bismark_methylation_extractor --no_overlap --comprehensive --bedGraph --cytosine_report --CX_context --genome_folder /<abs_path>/ref/${lineage}/ -p bam/${f}_1_bismark_bt2_pe_nonmismatch.bam -o bam/;
    gunzip bam/${f}_1_bismark_bt2_pe_nonmismatch.bismark.cov.gz;
    python3 magic.py bam/${f}_1_bismark_bt2_pe_nonmismatch.bismark.cov bam/${f}_1_bismark_bt2_pe_nonmismatch.CX_report.txt Me-information/${f}_result;done
    samtools sort bam/${f}_1_bismark_bt2_pe_nonmismatch.bam -o bam/${f}_1_bismark_bt2_pe_nonmismatch.sorted.bam;
    samtools index bam/${f}_1_bismark_bt2_pe_nonmismatch.sorted.bam;
done
~~~

#### Methylation activity

~~~shell
for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++)); do python3 MethFraction.py umi/${f}.fq bam/${f}_1_bismark_bt2_pe_nonmismatch.sorted.bam Me-information/${f}_result csv/${f}.csv;done
~~~

#### Histogram plot

~~~shell
for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++)); do python3 MethHist.py csv/${f}.csv png;done
~~~

## Copyright

copyright (c) Li Lab of HZAU, All rights reserved. 
