# Workflow (2022.04.16)

## BP

**输入**：一对 fq 文件 (fastq.fq) 、barcode 文本文件 (barcode.txt) 、选择或上传参考序列 (Reference) 、突变位点文件 (mutation.csv)

~~~shell
# 准备阶段：测序文件解压缩，创建目录
> gzip -dk <filename>.fq.gz
> mkdir qc fq bam csv xlsx png

# 数据质控，注意 --length_required 的参数为 114
> fastp -i <filename_1>.fq -I <filename_2>.fq -o qc/clean_1.fq -O qc/clean_2.fq -q 20 --length_required 114 -j qc/<filename>.json -h qc/<filename>.html
# 测序数据拆分
> python3 SeqPurifier.py BP barcode/<barcode>.txt qc/clean_1.fq qc/clean_2.fq fq

# 比对参考序列或参考基因组，并建立 bam 文件索引
> for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++)); do bwa mem -t 4 -M <Reference>.fa fq/${f}_1.fq fq/${f}_2.fq | samtools view -S -b - | samtools sort -o bam/${f}.sorted.bam -; done
> for b in $(ls bam/*); do samtools index ${b}; done

# 等位基因频率计算
> python3.7 AltFraction.py BP barcode/<barcode>.txt mutation/<mut>.csv bam csv
# 基因分型
> for t in $(ls csv/*); do python3.7 Minima.py ${t} xlsx png; done
~~~

**输出**：分型图 (markers.png) 、分型结果 (markers.xlsx)

~~~mermaid
flowchart TB
    
    %% Definition of style class
    classDef endPoint fill:#FFFCED,stroke:#B3AA82,stroke-width:2px
    classDef startPoint stroke-width:2px
    
    %% Subgraph of sequence analysis
    subgraph BP-Seq
    Fastp(Fastp) --> clean.fq --> SPpy([SeqPurifier.py])
    %% Diverge of clean.fq
    SPpy --> fqSample[[fq/samples.fq]]
    fqSample --> BWA(BWA) --> SAMtools(SAMtools) --> sampleBam[[samples.sorted.bam]]
    %% Converge of molecular marker and umi marker
    sampleBam --> AFpy([AltFraction.py])
    AFpy --> Mcsv[[markers.csv]] --> Mpy([Minima.py])
    end
    
    %% Input start
    SAservice((Input)):::startPoint --> fastq.fq
    SAservice:::startPoint --> barcode.txt
    SAservice:::startPoint --> mutation.csv
    SAservice:::startPoint -.Optional.-> RefSeq[(Reference)]
    
    %% Flow of user inputs
    fastq.fq --> Fastp
    barcode.txt --> SPpy
    barcode.txt --> AFpy
    mutation.csv --> AFpy
    RefSeq --> BWA
    
    %% Output
    Mpy --> Mpng[[markers.png]] -.genotyping.-> Mxlsx[[markers.xlsx]]
    Mpy --> Mxlsx --> result((Genotype)):::endPoint
~~~

## iBP

**输入**：一对 fq 文件 (fastq.fq) 、barcode 文本文件 (barcode.txt) 、选择或上传参考序列 (Reference) 、突变位点文件 (mutation.csv)

~~~shell
# 准备阶段：测序文件解压缩，创建目录
> gzip -dk <filename>.fq.gz
> mkdir qc fq umi bam csv xlsx png

# 数据质控
> fastp -i <filename_1>.fq -I <filename_2>.fq -o qc/clean_1.fq -O qc/clean_2.fq -q 20 --length_required 72 -j qc/<filename>.json -h qc/<filename>.html
# 测序数据拆分
> python3 SeqPurifier.py iBP barcode/<barcode>.txt qc/clean_1.fq qc/clean_2.fq fq umi

# 比对参考序列或参考基因组，并建立 bam 文件索引
> for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++)); do bwa mem -t 4 -M <Reference>.fa fq/${f}_1.fq fq/${f}_2.fq | samtools view -S -b - | samtools sort -o bam/${f}.sorted.bam -; done
> for b in $(ls bam/*); do samtools index ${b}; done

# 等位基因频率计算
> python3.7 AltFraction.py iBP barcode/<barcode>.txt mutation/<mut>.csv bam umi csv
# 基因分型
> for t in $(ls csv/*); do python3.7 Minima.py ${t} xlsx png; done
~~~

**输出**：分型图 (markers.png) 、分型结果 (markers.xlsx)

~~~mermaid
flowchart TB
    
    %% Definition of style class
    classDef endPoint fill:#FFFCED,stroke:#B3AA82,stroke-width:2px
    classDef startPoint stroke-width:2px
    
    %% Subgraph of sequence analysis
    subgraph iBP-Seq
    Fastp(Fastp) --> clean.fq --> SPpy([SeqPurifier.py])
    %% Diverge of clean.fq
    SPpy --> fqSample[[fq/samples.fq]]
    SPpy --> umiSample[[umi/samples.fq]]
    fqSample --> BWA(BWA) --> SAMtools(SAMtools) --> sampleBam[[samples.sorted.bam]]
    %% Converge of molecular marker and umi marker
    umiSample --> AFpy([AltFraction.py])
    sampleBam --> AFpy
    AFpy --> Mcsv[[markers.csv]] --> Mpy([Minima.py])
    end
    
    %% Input start
    SAservice((Input)):::startPoint --> fastq.fq
    SAservice:::startPoint --> barcode.txt
    SAservice:::startPoint --> mutation.csv
    SAservice:::startPoint -.Optional.-> RefSeq[(Reference)]
    
    %% Flow of user inputs
    fastq.fq --> Fastp
    barcode.txt --> SPpy
    barcode.txt --> AFpy
    mutation.csv --> AFpy
    RefSeq --> BWA
    
    %% Output
    Mpy --> Mpng[[markers.png]] -.genotyping.-> Mxlsx[[markers.xlsx]]
    Mpy --> Mxlsx --> result((Genotype)):::endPoint
~~~

## mBP

**输入**：一对 fq 文件 (fastq.fq) 、barcode 文本文件 (barcode.txt) 、选择或上传参考序列 (Reference) 

~~~shell
# 准备阶段：测序文件解压缩，创建目录
> gzip -dk <filename>.fq.gz
> mkdir qc fq umi sam bam csv xlsx png

# 数据质控
> fastp -i <filename_1>.fq -I <filename_2>.fq -o qc/clean_1.fq -O qc/clean_2.fq -q 20 --length_required 72 -j qc/<filename>.json -h qc/<filename>.html
# 测序数据拆分
> python3 SeqPurifier.py mBP barcode/<barcode>.txt qc/clean_1.fq qc/clean_2.fq fq umi

# 参考序列按自交系拆分，并对每种类型的参考序列建立文件索引
> python3 MethReference.py barcode/<barcode>.txt ref/<ref>.fa
> for i in $(ls -d ref/*);do bismark_genome_preparation ${i}/ ;done
~~~

~~~shell
# 比对参考基因组，注意这里参考基因组作为参数输入应为绝对路径
> for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++)); do lineage=$(sed -n $[${f}+1]p barcode/<barcode>.txt | sed "s/.*,//" | sed "s/\r//" | sed "s/\//_/g");bismark -N 0 -L 20 --bowtie2 --non_bs_mm /<abs_path>/ref/${lineage}/ --sam -o sam/ --q -1 fq/${f}_1.fq -2 fq/${f}_2.fq; done
~~~

~~~
# SAM 文件处理（外层循环相同，可考虑在一个循环内用分号换行执行下列命令）
> for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++)); do grep "@" sam/${f}_1_bismark_bt2_pe.sam > sam/${f}-header; done
> for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++)); do grep "XA:Z:0" sam/${f}_1_bismark_bt2_pe.sam > sam/${f}-temp1; done
> for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++)); do grep "XB:Z:0" sam/${f}_1_bismark_bt2_pe.sam > sam/${f}-temp2; done
> for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++)); do awk 'NR==FNR {f1[$1]=$0;next} {idx=$1;if(idx in f1)print$0}' sam/${f}-temp1 sam/${f}-temp2 > sam/${f}-temp4; done
> for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++)); do awk 'NR==FNR {f1[$1]=$0;next} {idx=$1;if(idx in f1)print$0}' sam/${f}-temp2 sam/${f}-temp1 > sam/${f}-temp3; done
> for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++)); do cat sam/${f}-header sam/${f}-temp3 sam/${f}-temp4 > sam/${f}_1_bismark_bt2_pe_nonmismatch.sam; done

# SAM 转 BAM，并提取甲基化信息
> for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++)); do samtools sort -n -O bam -o bam/${f}_1_bismark_bt2_pe_nonmismatch.bam sam/${f}_1_bismark_bt2_pe_nonmismatch.sam; done
> for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++)); do lineage=$(sed -n $[${f}+1]p barcode/barcode.txt | sed "s/.*,//" | sed "s/\r//" | sed "s/\//_/g"); bismark_methylation_extractor --no_overlap --comprehensive --bedGraph --cytosine_report --CX_context --genome_folder /<a>/ref/${lineage}/ -p bam/${f}_1_bismark_bt2_pe_nonmismatch.bam -o bam/;done
> for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++)); do gunzip bam/${f}_1_bismark_bt2_pe_nonmismatch.bismark.cov.gz;done
> for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++)); do python3 magic.py bam/${f}_1_bismark_bt2_pe_nonmismatch.bismark.cov bam/${f}_1_bismark_bt2_pe_nonmismatch.CX_report.txt Me-information/${f}_result;done
> for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++)); do samtools sort bam/${f}_1_bismark_bt2_pe_nonmismatch.bam -o bam/${f}_1_bismark_bt2_pe_nonmismatch.sorted.bam;done
> for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++)); do samtools index bam/${f}_1_bismark_bt2_pe_nonmismatch.sorted.bam;done
~~~

~~~shell
# 计算甲基化位点活性，得到甲基化分型结果
> for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++)); do python3 MethFraction.py umi/${f}.fq bam/${f}_1_bismark_bt2_pe_nonmismatch.sorted.bam Me-information/${f}_result csv/${f}.csv;done
~~~

~~~shell
# 甲基化位点活性图
> for ((f=0; f<$(ls -l fq | grep "_1.fq$" | wc -l); f++)); do python3 MethHist.py csv/${f}.csv png;done
~~~

**输出**：甲基化位点活性图 (active.png) 、甲基化分型结果 (markers.csv)

~~~mermaid
flowchart TB
    
    %% Definition of style class
    classDef endPoint fill:#FFFCED,stroke:#B3AA82,stroke-width:2px
    classDef startPoint stroke-width:2px
    
    %% Subgraph of sequence analysis
    subgraph epi-Genotyping
    Fastp(Fastp) --> clean.fq --> SPpy([SeqPurifier.py])
    %% Diverge of clean.fq
    SPpy --> fqSample[[fq/samples.fq]]
    SPpy --> umiSample[[umi/samples.fq]]
    fqSample --> Bismark(Bismark) --> Sam[[sam]] --> Samtools(Samtools) --> Bam[[Bam]] --> extractor([bismark_methylation_extractor]) --> sampleTxt[[text]] --> magic([magic.py]) --> methTxt[[text]]
    %% Converge of molecular marker and umi marker
    umiSample --> MFpy([MethFraction.py])
    methTxt --> MFpy
    %% Plot unit
    MethHist([MethHist.py])
    end
    
    %% Subgraph of sequence analysis
    subgraph sub-Referecne
    MethReference([MethReference.py]) --> subRefSeq[(Sub-Reference)] --> Bismark
    end
    
    %% Input start
    SAservice((Input)):::startPoint --> fastq.fq
    SAservice:::startPoint --> barcode.txt
    SAservice:::startPoint -.Optional.-> RefSeq[(Reference)]
    
    %% Flow of user inputs
    fastq.fq --> Fastp
    barcode.txt --> SPpy
    barcode.txt --> MFpy
    RefSeq --> MethReference
    
    MFpy --> Mcsv[[markers.csv]]
    Mcsv[[markers.csv]] --> MethHist --> Mpng[[active.png]]
    Mpng[[active.png]] --> result((epi-Genotype)):::endPoint
    Mcsv[[markers.csv]] --> result((epi-Genotype)):::endPoint
    result((epi-Genotype)):::endPoint
~~~

