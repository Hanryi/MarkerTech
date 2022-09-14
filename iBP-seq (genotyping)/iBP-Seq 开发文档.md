# 需求分析

- `质控` 对于测序数据，首先进行质控，过滤掉低质量的测序 reads ，过滤后的测序数据用于后续分析，后续分析仅对于一对无重复 barcode 序列的双端测序文件而言，而具有重复 barcode 序列的多对测序文件则不在本次需求的讨论范围。
- `编写程序` 使之能够将混合测序数据拆分为若干不同样本的测序数据，同时对该程序有以下几点要求：
  1. 根据 barcode 序列和内置的 bridge 序列完成混合数据分流，在样本识别时要求 barcode 序列精确匹配，而在 bridge 序列模糊匹配，即要求最多存在 0 碱基的错配（允许错配碱基为 N 或其他任意字符）；
  2. 每条完成样本识别的 read 要先将非基因组序列区段 (barcode + bridge + UMI) 和对应的质量值字符串切除，再输出到 fq 目录下样本编号从 0 开始的一对 fastq 文件中，作为编号对应样本的测序数据；
  3. 每条完成样本识别的 read 要同时将测序 ID 、切除的 UMI 序列、+ 分隔符和切除 UMI 序列对应的质量值字符串输出到 umi 目录下样本编号从 0 开始的一个 fastq 文件中，没有丢失 UMI 相关的任何信息；
  3. 对混合数据分流的结果进行百分比统计，体现分流的质量
- `比对` 用 BWA 将 fq 目录下每个样品的双端测序文件比对到参考序列上，避免了非基因组序列的干扰，提高了比对准确度；用 SAMtools 将比对后的 SAM 文件进行排序，生成对应编号的 sorted.bam 文件。该步骤可以用管道符，以避免产生中间文件。
- `编写程序` 使之能够根据该样本的 reads 信息，计算每个样本在各变异位点上的等位基因频率，同时做如下说明：

  1. 根据  barcode 文件中的 barcode 数量作为样本数，无需关注任何具体的 barcode 序列；
  2. 根据每条 read ID 可以从对应编号的 BAM 文件中确认对应 ID 的 read 各碱基的位置信息，并判断变异的序列类型；也可以用 ID 从 umi 目录下对应编号的 fastq 中获取 UMI 序列信息；
  3. 该程序需要允许存储有变异位信息的文件至少以 csv, txt, xlsx 和 vcf 格式存在；
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

# 混合数据分流

根据[需求分析](# 需求分析)中第二项的要求，该程序不仅需要在样本层面，完成不同测序数据的分流；对于同一个样本也要序列层面，将基因组序列与非基因组序列分离开，SeqPurifier 具有两层拆分的功能。

在测序数据分流的过程中，对于测序 read 的判读采用每四行作为一个 read 单位的读取方式，对同一 read 中的各行按需进行判断和处理，将符合条件的 read 输出到文件中。根据要求，在判读过程中对 read 提出的条件应有如下项目：

## 特异序列识别

SeqPurifier 的输入有 barcode 文件和一对双端测序文件，barcode 文件里的 barcode 序列是样本识别的关键特征；而程序内置的 bridge 序列信息则是辅助识别 barcode 序列的标志序列，以提高 barcode 识别的精确度，减少在 read 两端都能匹配到 barcode 的情况，但过度严格的匹配又会使得 read 两端都匹配不到的结果增加，对 bridge 的识别设置模糊匹配机制，默认允许 bridge 序列区段存在最多一个碱基。

定义 sample_num(read_seq, tolerance=1) 方法，必要参数 read_seq 对于输入的测序序列进行判断，是否测序序列是否符合前 8bp 的 barcode 序列精确匹配同时 9-18bp 的 bridge 序列模糊匹配（最多错配一个碱基）的匹配模式。若匹配，则根据匹配到的 barcode 序列返回样本编号；否则不匹配，则返回 None 。默认参数 tolerance 设置默认允许的最多错配碱基数量为 0 。

使用 sample_num 方法对一对双端测序 read 进行判断，只有一个返回样本编号，另一个返回 None 值的情况才是正确的情况。并对于 _1 端和 _2 端返回样本编号的情况分别进行标记，用 cache 变量存储两种情况，便于后续为拆分非基因组序列指明需要处理的那一端 read 。

## 基因组序列拆分

对于符合要求的 read ，将 _1 端和 _2 端序列分别输出，并将其中包含附加序列一端的基因组和非基因组序列分别输出，这一端的部分基因组序列和另一端全部的基因组序列输出到 fq 目录下，这一端的部分非基因组序列中取出 UMI 序列输出到 umi 目录下，所有的输出都遵照 fastq 文件格式，包含 ID、输出的序列、分隔符和对应的质量值。

定义 fq_assign(num, tag, f_read, r_read) 方法进行 fastq 文件的输出，该方法需要三个必要输入参数。num 参数指定当前输出的序列到指定样本编号的文件；tag 参数用于传递 cache 变量的数值，指定需要拆分哪一端的 read 序列，另一端的序列无需拆分，按原样输出；r_read 和 r_read 参数分别为 forward read 和 reverse read 的列表。该方法仅进行了 fastq 文件输出的操作，返回值为 None 。

在接收到非 0 的 cache 数值后，调用 fa_assign 方法进行文件输出。在输出前先对测序序列和相应的质量值序列进行长度的比较，确保长度相等，否则不进行后续的输出。在保证两行序列相等的情况下，根据 tag 值动态选择 f_read 或者 r_read 作为序列拆分的处理对象。该方法将一对测序序列分为三份文件输出，其中两份文件是从一端测序文件中拆分出来的，另一份文件与另一端测序文件一致，三份文件都按 fastq 文件格式输出，文件换行符为 '\n'。此外，UMI 序列的输出也遵照 fastq 文件格式，提供的 ID 信息为后续查找 UMI 序列提供了便利，在没有丢失测序信息的情况下提高了该问题下的处理效率。
