# 需求分析

- `质控` 对于测序数据，首先进行质控，过滤掉低质量的测序 read ，过滤后的测序数据由于后续分析，后续分析仅对于一对无重复 barcode 序列的双端测序文件而言，具有重复 barcode 序列的多对测序文件不在本次需求的讨论范围。

- `编写程序` 使之能够将混合测序数据拆分为若干不同样本的测序数据，同时对该程序有以下几点要求：

  1. 根据 barcode 序列和内置的 bridge 序列完成混合数据分流时，在样本识别时要求 barcode 序列精确匹配，而在 bridge 序列模糊匹配，即要求最多存在一个碱基的错配（允许错配碱基为 N 或其他任意字符）
  2. 每条完成样本识别的 read 要先将非基因组序列区段 (barcode + bridge + UMI) 和对应的质量值字符串切除，再输出到 fq 目录下样本编号从 0 开始的双端 fastq 文件中，作为编号对应样本的测序数据
  3. 每条完成样本识别的 read 要同时将测序 ID 、切除的 UMI 序列、+ 分隔符和切除 UMI 序列对应的质量值字符串输出到 umi 目录下样本编号从 0 开始的双端 fastq 文件中，没有丢失 UMI 相关的任何信息

  因此，该程序不仅需要完成不同样本数据的分流，对于同一个样本也要将基因组序列与非基因组序列分离开，需要具有两层拆分的功能。

- `比对` 用 BWA 将 fq 目录下每个样品的双端测序文件比对到参考序列上，避免了非基因组序列的干扰，提高了比对准确度；用 SAMtools 将比对后的 SAM 文件进行排序，生成对应编号的 sorted.bam 文件。该步骤可以用管道符，以避免产生中间文件。

- `编写程序` 使之能够根据该样本的 reads 信息，计算每个样本在各变异位点上的等位基因频率，同时做如下说明：

  1. 根据  barcode 文件中的 barcode 数量作为样本数，无需关注任何具体的 barcode 序列
  2. 根据每条 read ID 可以从对应编号的 BAM 文件中确认对应 ID 的 read 各碱基的位置信息，并判断变异的序列类型；也可以用 ID 从 umi 目录下对应编号的 fastq 中获取 UMI 序列信息
  3. 该程序需要允许存储有变异位信息的文件至少以 csv, txt, xlsx 和 vcf 格式存在
  4. 对于每个变异位点，在计算每个样品的等位基因频率的同时统计不重复的 UMI 序列数量，并将该位点的所有样本对应的信息输出为 csv 格式的表格文件。该文件具有表头 (Num, Frequency, UMI) 。

- `编写程序` 使之能够根据每个变异位点输出的 Frequency 信息进行核密度估计，以 KDE (kernel density estimation) 曲线的若干极值点作为分型阈值，根据 Frequency 的数值，将实验群体划分为不同的基因型或甲基化类型。

## Flowchart

由上述分析表述构建如下分析流程，其中 PrimerDesign Service 和 SeqAnalysis Service 是两个开始节点，genotype 是结束节点：

~~~mermaid
flowchart TB
    
    %% Definition of style class
    classDef endPoint fill:#FFFCED,stroke:#B3AA82,stroke-width:2px
    classDef startPoint stroke-width:2px
    
    %% Subgraph of primer design
    subgraph PrimerDesign
    PDpy([PrimerDesigner.py])
    end
    
    %% Subgraph of sequence analysis
    subgraph SeqAnalysis
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
    
    %% Primer design
    PDservice((PrimerDesign<br>Service)):::startPoint --> template.fa
    template.fa --> PDpy --> primers.xlsx
    %% Wet experiment
    primers.xlsx -.Amplicon<br>Sequencing.-> fastq.fq --> Fastp
    %% Sequence analysis
    SAservice((SeqAnalysis<br>Service)):::startPoint --> fastq.fq
    SAservice:::startPoint --> barcode.txt
    SAservice:::startPoint --> mutation.csv
    SAservice:::startPoint -.Optional.-> RefSeq[(Reference)]
    
    %% Flow of user inputs
    barcode.txt --> SPpy
    barcode.txt --> AFpy
    mutation.csv --> AFpy
    RefSeq --> BWA
    
    Mpy --> Mpng[[markers.png]] -.genotyping.-> Mxlsx[[markers.xlsx]]
    Mpy --> Mxlsx --> result((Genotype)):::endPoint
~~~

