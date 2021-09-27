

# 1. Pro-seq/Gro-seq Analyis



### PEPPRO 

<font color = blue><font size= 8> 2021.06.18 </font></font>

* Install PEPPRO

```shell
# git clone
git clone https://github.com/databio/peppro.git
cd peppro
# create env
conda create -n peppro
conda install --file requirements.txt
#### requirements ###
#cutadapt>=2.9
#numpy
#pandas>=0.20.2
#piper>=0.12.1
#refgenconf
#refgenie>=0.9.1
pip install pararead
pip install looper

# install R library
Rscript -e 'install.packages("devtools")'
Rscript -e 'devtools::install_github("pepkit/pepr")'
Rscript -e 'install.packages("BiocManager")'
Rscript -e 'BiocManager::install("GenomicRanges")'
Rscript -e 'devtools::install_github("databio/GenomicDistributions")'
Rscript -e 'BiocManager::install(c("BSgenome", "GenomicFeatures", "ensembldb"))'
Rscript -e 'install.packages("http://big.databio.org/GenomicDistributionsData/GenomicDistributionsData_0.0.1.tar.gz", repos=NULL)'
```

* Data Description

  | Dataset   | CellLine | Treatment | Experimrnt | Replicates |
  | --------- | -------- | --------- | ---------- | ---------- |
  | GSE128081 | A375     | DMSO      | Pro-seq    | 3          |
  | GSE57430  | A375     | DMSO      | Gro-seq    | 2          |

* Run Peppro

```shell
# Pro-seq
python ./peppro/pipelines/peppro.py \
  --sample-name A375_Proseq_r1 \
  --genome hg38 \
  --prealignments human_rDNA \
  --input /h/menghan/eRNA/Pro-seq/GSE128081/SRR8706163_GSM3664682_dmso_r1_Homo_sapiens_OTHER_1.fastq.gz /h/menghan/eRNA/Pro-seq/GSE128081/SRR8706163_GSM3664682_dmso_r1_Homo_sapiens_OTHER_2.fastq.gz\
  --single-or-paired paired-end \
  -P 40 \
  -O /h/menghan/eRNA/Pro-seq/res  # replicate 1
  
python ./peppro/pipelines/peppro.py \
  --sample-name A375_Proseq_r2 \
  --genome hg38 \
  --prealignments human_rDNA \
  --input /h/menghan/eRNA/Pro-seq/GSE128081/SRR8706164_GSM3664683_dmso_r2_Homo_sapiens_OTHER_1.fastq.gz /h/menghan/eRNA/Pro-seq/GSE128081/SRR8706164_GSM3664683_dmso_r2_Homo_sapiens_OTHER_2.fastq.gz\
  --single-or-paired paired-end \
  -P 40 \
  -O /h/menghan/eRNA/Pro-seq/res  # replicate 2
  
  python ./peppro/pipelines/peppro.py \
  --sample-name A375_Proseq_r3 \
  --genome hg38 \
  --prealignments human_rDNA \
  --input /h/menghan/eRNA/Pro-seq/GSE128081/SRR8706165_GSM3664684_dmso_r3_Homo_sapiens_OTHER_1.fastq.gz /h/menghan/eRNA/Pro-seq/GSE128081/SRR8706165_GSM3664684_dmso_r3_Homo_sapiens_OTHER_2.fastq.gz \
  --single-or-paired paired-end \
  -P 40 \
  -O /h/menghan/eRNA/Pro-seq/res  # replicate 3
  
  
  # Gro-seq
  
  python ./peppro/pipelines/peppro.py \
  --sample-name A375_Groseq_r1 \
  --genome hg38 \
  --prealignments human_rDNA \
  --input /h/menghan/eRNA/Gro-seq/GSE57430/SRR1275489_GSM1382433_48h_DMSO_treated_A375_cells_replicate_a_Homo_sapiens_RNA-Seq.fastq.gz \
  --single-or-paired Single \
  -P 40 \
  --protocol gro \
  -O /h/menghan/eRNA/Gro-seq/res  # replicate a
  
python ./peppro/pipelines/peppro.py \
  --sample-name A375_Groseq_r2 \
  --genome hg38 \
  --prealignments human_rDNA \
  --input /h/menghan/eRNA/Gro-seq/GSE57430/SRR1275490_GSM1382434_48h_DMSO_treated_A375_cells_replicate_b_Homo_sapiens_RNA-Seq.fastq.gz \
  --single-or-paired Single \
  -P 40 \
  --protocol gro \
  -O /h/menghan/eRNA/Gro-seq/res  # replicate b
```

* QC report

1.  Pro-seq

![image-20210622220027868](/Users/menghan/Library/Application Support/typora-user-images/image-20210622220027868.png)

![image-20210622221155532](/Users/menghan/Library/Application Support/typora-user-images/image-20210622221155532.png)

![image-20210622221834643](/Users/menghan/Library/Application Support/typora-user-images/image-20210622221834643.png)

2. Gro-seq

![image-20210624165455237](/Users/menghan/Library/Application Support/typora-user-images/image-20210624165455237.png)

![image-20210624170813412](/Users/menghan/Library/Application Support/typora-user-images/image-20210624170813412.png)

![image-20210728115747090](/Users/menghan/Library/Application Support/typora-user-images/image-20210728115747090.png)

<font color = blue><font size= 8> 2021.06.24 </font></font>

### dREG

![image-20210624171214805](/Users/menghan/Library/Application Support/typora-user-images/image-20210624171214805.png)

![image-20210624171113519](/Users/menghan/Library/Application Support/typora-user-images/image-20210624171113519.png)

# 2. Pro-seq/Gro-seq Signal in Hanleng reference

<font color = blue><font size= 8> 2021.06.27 </font></font>

![image-20210624171332793](/Users/menghan/Library/Application Support/typora-user-images/image-20210624171332793.png)

![image-20210624171353783](/Users/menghan/Library/Application Support/typora-user-images/image-20210624171353783.png)





## 

![image-20210728144425996](/Users/menghan/Library/Application Support/typora-user-images/image-20210728144425996.png)

<font color = blue><font size= 8> 2021.07.01 </font></font>

# 3. A375 eRNA Analysis



## Define A375 eRNA

![image-20210728144531201](/Users/menghan/Library/Application Support/typora-user-images/image-20210728144531201.png)

![image-20210728144552054](/Users/menghan/Library/Application Support/typora-user-images/image-20210728144552054.png)



<font color = blue><font size= 8> 2021.07.05 </font></font>

## Nascent Signal in A375 eRNA candidate

### 1). Pro-seq Signal

![image-20210728144906899](/Users/menghan/Library/Application Support/typora-user-images/image-20210728144906899.png)

![image-20210728144919054](/Users/menghan/Library/Application Support/typora-user-images/image-20210728144919054.png)

### 2). Gro-seq Signal

![image-20210728145053237](/Users/menghan/Library/Application Support/typora-user-images/image-20210728145053237.png)

![image-20210728145100263](/Users/menghan/Library/Application Support/typora-user-images/image-20210728145100263.png)

​																										*Gro-seq数据质量较差，后续不进行分析*



## Functional Annotation

<font color = blue><font size= 8> 2021.07.07 </font></font>

### 1). Assocatied with SNV hotspots

![image-20210728151441703](/Users/menghan/Library/Application Support/typora-user-images/image-20210728151441703.png)

![image-20210728151457099](/Users/menghan/Library/Application Support/typora-user-images/image-20210728151457099.png)

### 2). Associated with A375 super enhancer

![image-20210728151538678](/Users/menghan/Library/Application Support/typora-user-images/image-20210728151538678.png)

![image-20210728151546170](/Users/menghan/Library/Application Support/typora-user-images/image-20210728151546170.png)

![image-20210728163057368](/Users/menghan/Library/Application Support/typora-user-images/image-20210728163057368.png)

<font color = blue><font size= 8> 2021.07.09 </font></font>

​	

## Epigenetic Marker

![image-20210728162827664](/Users/menghan/Library/Application Support/typora-user-images/image-20210728162827664.png)

### 1). H3K4me1 Marker

![image-20210728162916630](/Users/menghan/Library/Application Support/typora-user-images/image-20210728162916630.png)

### 2). H3K4me3 Marker

![image-20210728162942374](/Users/menghan/Library/Application Support/typora-user-images/image-20210728162942374.png)



## Target Gene Annotation

![image-20210728163037737](/Users/menghan/Library/Application Support/typora-user-images/image-20210728163037737.png)





<font color = blue><font size= 8> 2021.07.15 </font></font>

## eRNA in CRIRSPi Screen 

### 1). Cell growth context

![image-20210728164455329](/Users/menghan/Library/Application Support/typora-user-images/image-20210728164455329.png)

### 2). Drug resistance context

![image-20210728164530217](/Users/menghan/Library/Application Support/typora-user-images/image-20210728164530217.png)



## Result

![image-20210728164736897](/Users/menghan/Library/Application Support/typora-user-images/image-20210728164736897.png)



# 4. sgRNA Library Design

<font color = blue><font size= 8> 2021.07.21 </font></font>

* tool: Cas13design.v0.2 (https://gitlab.com/sanjanalab/cas13/-/blob/master/Cas13designGuidePredictor/Readme.md)

  ![image-20210728171143971](/Users/menghan/Library/Application Support/typora-user-images/image-20210728171143971.png)

## A375 genome-wide library

```sh
#!/bin/bash


for file in /h/menghan/eRNA/Cas13design/fasta/*.fasta
#for file in *
do
   Rscript ../scripts/RfxCas13d_GuideScoring.R $file ../data/Cas13designGuidePredictorInput.csv true
done
```

### 1). Output Example

![image-20210728171744855](/Users/menghan/Library/Application Support/typora-user-images/image-20210728171744855.png)

![image-20210728171759664](/Users/menghan/Library/Application Support/typora-user-images/image-20210728171759664.png)

### 2). Combine and Filter

<font color = blue><font size= 8> 2021.07.23 </font></font>

*长转录本设计4个sgRNA，短转录本设计3个sgRNA*

![image-20210728171944130](/Users/menghan/Library/Application Support/typora-user-images/image-20210728171944130.png)

![image-20210728172126262](/Users/menghan/Library/Application Support/typora-user-images/image-20210728172126262.png)

![image-20210728180018211](/Users/menghan/Library/Application Support/typora-user-images/image-20210728180018211.png)

<font color = blue><font size= 8> 2021.07.24 </font></font>

## Specific eRNA Library

### 1) .Tiling Library

![image-20210728222241639](/Users/menghan/Library/Application Support/typora-user-images/image-20210728222241639.png)

### 2). Mismatch Library

![image-20210728222332469](/Users/menghan/Library/Application Support/typora-user-images/image-20210728222332469.png)

### 3). Lengthen&Shorten Library

![image-20210728222419615](/Users/menghan/Library/Application Support/typora-user-images/image-20210728222419615.png)

