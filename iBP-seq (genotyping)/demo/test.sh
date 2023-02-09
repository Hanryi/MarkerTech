#!/bin/bash
gzip -dk 96raw_0.1G_R*.fq.gz
mkdir qc fq umi bam csv xlsx png

# Quality control
fastp -i 96raw_0.1G_R1.fq -I 96raw_0.1G_R2.fq -o qc/clean_1.fq -O qc/clean_2.fq -q 20 --length_required 72 -j qc/96raw_0.1G.json -h qc/96raw_0.1G.html

# Data division
python3 SeqPurifier.py iBP barcode/barcode.txt qc/clean_1.fq qc/clean_2.fq fq umi

# Reference alignment
for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++));
do
    bwa mem -t 4 -M Zea_mays.AGPv4.dna.toplevel.fa fq/${f}_1.fq fq/${f}_2.fq | samtools view -S -b - | samtools sort -o bam/${f}.sorted.bam -;
done
for b in $(ls bam/*);
do
    samtools index ${b};
done

# Marker frequency
python3 AltFraction.py iBP barcode/barcode.txt mutation/getSNP.csv bam umi csv

# KDE Cluster
for t in $(ls csv/*);
do
    python3 Minima.py ${t} xlsx png;
done