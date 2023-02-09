~~~shell
fastp -i BSPCR_R1.fq.gz -I BSPCR_R2.fq.gz -o qc/clean_1.fq -O qc/clean_2.fq -q 20 --length_required 72


bismark -N 0 -L 20 --bowtie2 --non_bs_mm ref/ --sam -o sam/ --q -1 qc/clean_1.fq -2 qc/clean_2.fq

grep "@" sam/clean_1_bismark_bt2_pe.sam > sam/clean-header
grep "XA:Z:0" sam/clean_1_bismark_bt2_pe.sam > sam/clean-temp1
grep "XB:Z:0" sam/clean_1_bismark_bt2_pe.sam > sam/clean-temp2
awk 'NR==FNR {f1[$1]=$0;next} {idx=$1;if(idx in f1)print$0}' sam/clean-temp1 sam/clean-temp2 > sam/clean-temp4
awk 'NR==FNR {f1[$1]=$0;next} {idx=$1;if(idx in f1)print$0}' sam/clean-temp2 sam/clean-temp1 > sam/clean-temp3
cat sam/clean-header sam/clean-temp3 sam/clean-temp4 > sam/clean_1_bismark_bt2_pe_nonmismatch.sam


samtools sort -n -O bam -o bam/clean_1_bismark_bt2_pe_nonmismatch.bam sam/clean_1_bismark_bt2_pe_nonmismatch.sam
bismark_methylation_extractor --no_overlap --comprehensive --bedGraph --cytosine_report --CX_context --genome_folder ref/ -p bam/clean_1_bismark_bt2_pe_nonmismatch.bam -o bam/
gunzip bam/clean_1_bismark_bt2_pe_nonmismatch.bismark.cov.gz
python3 magic.py bam/clean_1_bismark_bt2_pe_nonmismatch.bismark.cov bam/clean_1_bismark_bt2_pe_nonmismatch.CX_report.txt Me-information/clean_result;done
samtools sort bam/clean_1_bismark_bt2_pe_nonmismatch.bam -o bam/clean_1_bismark_bt2_pe_nonmismatch.sorted.bam
samtools index bam/clean_1_bismark_bt2_pe_nonmismatch.sorted.bam

awk '{print $4 "\t" $5 "\t" $6}' bam/clean_1_bismark_bt2_pe_nonmismatch.bismark.cov > 1.txt
cat Me-information/clean_result bam/1.txt > 2.txt
~~~

即为最终得到的甲基化水平

