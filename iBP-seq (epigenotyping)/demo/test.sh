#!/bin/bash
gzip -dk outclean-R*.fq.gz
mkdir qc fq umi sam bam csv xlsx png

# Quality control
fastp -i outclean-R1.fq -I outclean-R2.fq -o qc/clean_1.fq -O qc/clean_2.fq -q 20 --length_required 72 -j qc/outclean.json -h qc/outclean.html

# Data division
python3 SeqPurifier.py mBP barcode/barcode.txt qc/clean_1.fq qc/clean_2.fq fq umi

# Reference alignment
python3 MethReference.py barcode/barcode.txt ref/BSseq-ref.fa
for i in $(ls -d ref/*);
do
    bismark_genome_preparation ${i}/ ;
done
for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++));
do
    ref_name=$(sed -n $[${f}+1]p barcode/barcode.txt | sed "s/.*,//" | sed "s/\r//");
    lineage=$(python3 -c "print('/'.join(sorted('${ref_name}'.split('/'))))" | sed "s/\//_/g");
    bismark -N 0 -L 20 --bowtie2 --non_bs_mm ref/"${lineage}"/ --sam -o sam/ --q -1 fq/${f}_1.fq -2 fq/${f}_2.fq;
done

# SAM processing
for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++)); 
do
    grep "@" sam/${f}_1_bismark_bt2_pe.sam > sam/${f}-header;
    grep "XA:Z:0" sam/${f}_1_bismark_bt2_pe.sam > sam/${f}-temp1;
    grep "XB:Z:0" sam/${f}_1_bismark_bt2_pe.sam > sam/${f}-temp2;
    awk 'NR==FNR {f1[$1]=$0;next} {idx=$1;if(idx in f1)print$0}' sam/${f}-temp1 sam/${f}-temp2 > sam/${f}-temp4;
    awk 'NR==FNR {f1[$1]=$0;next} {idx=$1;if(idx in f1)print$0}' sam/${f}-temp2 sam/${f}-temp1 > sam/${f}-temp3;
    cat sam/${f}-header sam/${f}-temp3 sam/${f}-temp4 > sam/${f}_1_bismark_bt2_pe_nonmismatch.sam; 
done
for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++));
do
    samtools sort -n -O bam -o bam/${f}_1_bismark_bt2_pe_nonmismatch.bam sam/${f}_1_bismark_bt2_pe_nonmismatch.sam;
    lineage=$(sed -n $[${f}+1]p barcode/barcode.txt | sed "s/.*,//" | sed "s/\r//" | sed "s/\//_/g");
    bismark_methylation_extractor --no_overlap --comprehensive --bedGraph --cytosine_report --CX_context --genome_folder ref/${lineage}/ -p bam/${f}_1_bismark_bt2_pe_nonmismatch.bam -o bam/;
    gunzip bam/${f}_1_bismark_bt2_pe_nonmismatch.bismark.cov.gz;
    python3 magic.py bam/${f}_1_bismark_bt2_pe_nonmismatch.bismark.cov bam/${f}_1_bismark_bt2_pe_nonmismatch.CX_report.txt Me-information/${f}_result;done
    samtools sort bam/${f}_1_bismark_bt2_pe_nonmismatch.bam -o bam/${f}_1_bismark_bt2_pe_nonmismatch.sorted.bam;
    samtools index bam/${f}_1_bismark_bt2_pe_nonmismatch.sorted.bam;
done

# Methylation activity
for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++));
do
    python3 MethFraction.py umi/${f}.fq bam/${f}_1_bismark_bt2_pe_nonmismatch.sorted.bam Me-information/${f}_result csv/${f}.csv;
done

# Histogram plot
for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++));
do
    python3 MethHist.py csv/${f}.csv png;
done