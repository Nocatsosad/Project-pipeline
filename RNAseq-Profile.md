# RNA-seq Profile of Melanoma BRAFi-resistance Patients/Cells 

**Gene of Interest:**

***MIR17HG***
***TP53***
***MYC***
***E2F3***
***STAT3***
***IRF1***

***WWC1***

***TGFB2***

***TEAD3***

***TGFBR1***

<font color = blue><font size= 8> 2021.05.05 </font></font>

## GSE50509

* Data Description

  **`Patients Pre untreated` VS  `Patients Prog Treated with vemurafenib`**

* Download Data 

```R
data <- getGEO('GSE50509',GSEMatrix=TRUE,AnnotGPL=TRUE,destdir="./RNA-seq/Public/GSE50509")
# 下载表达信息
exprSet = exprs(data[[1]])
write.table(exprSet,'./RNA-seq/Public/GSE50509/expression.txt',sep='\t')

# 下载探针注释信息
gpl<-getGEO('GPL10558',destdir ="./RNA-seq/Public/GSE50509")
genename <- Table(gpl)[,c(1,10,13)] 
write.table(genename,'./RNA-seq/Public/GSE50509/genename.txt',sep='\t',row.names = FALSE)
```

* Gene Expression (interest)

![image-20210505162133396](/Users/menghan/Library/Application Support/typora-user-images/image-20210505162133396.png)

![image-20210511215958798](/Users/menghan/Library/Application Support/typora-user-images/image-20210511215958798.png)

## GSE68840

* Data Description

  **`LM16`  VS  `LM16VR`**

* Download Data

```R
data <- getGEO('GSE68840',GSEMatrix=TRUE,AnnotGPL=TRUE,destdir="./RNA-seq/Public/GSE68840")
# 下载表达信息
exprSet = exprs(data[[1]])
write.table(exprSet,'./RNA-seq/Public/GSE68840/expression.txt',sep='\t')

# 下载探针注释信息
gpl<-getGEO('GPL13938',destdir ="./RNA-seq/Public/GSE68840")
genename <- Table(gpl)[,c(1,9,12)]
write.table(genename,'./RNA-seq/Public/GSE68840/genename.txt',sep='\t',row.names = FALSE)
```

* Differential Expression(interest)

| Symbol | logFC      | AveExpr    | t          | P.Value    | adj.P.Val  | B          |
| ------ | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- |
| IRF1   | 1.64862547 | 8.12442048 | 13.6515047 | 3.31E-07   | 5.01E-06   | 6.71723279 |
| TGFBR1 | -0.5576797 | 8.43233612 | -9.8780654 | 4.81E-06   | 3.83E-05   | 3.77831795 |
| STAT3  | 0.60710969 | 11.2199148 | 9.62207309 | 5.95E-06   | 4.49E-05   | 3.54495727 |
| MYC    | -0.4001254 | 12.7422201 | -5.5202888 | 0.0004064  | 0.00131145 | -1.0709995 |
| WWC1   | 0.71297508 | 6.73712102 | 5.10369655 | 0.00069677 | 0.00205404 | -1.6541302 |
| TGFB2  | 0.2434645  | 11.1488783 | 4.21660766 | 0.00238666 | 0.00582255 | -2.9757977 |
| TEAD3  | 0.59563383 | 8.33477887 | 4.14377078 | 0.00265398 | 0.00636459 | -3.0889454 |
| E2F3   | -0.1335483 | 10.6528223 | -2.0935726 | 0.06661504 | 0.10367317 | -6.4024698 |
| TP53   | -0.1284795 | 9.29450766 | -1.105471  | 0.29836895 | 0.37414806 | -7.7441967 |

<font color = blue><font size= 8> 2021.05.06 </font></font>

## GSE99923

* Data Description

  **`A375`  VS  `A375VR`**

* Download Data

![image-20210506114704544](/Users/menghan/Library/Application Support/typora-user-images/image-20210506114704544.png)

* Differential Expression

| Symbol  | logFC      | logCPM     | PValue     |
| ------- | ---------- | ---------- | ---------- |
| E2F3    | 0.18929224 | 6.08268454 | 0.00013803 |
| HRAS    | -0.5141126 | 4.88737211 | 1.66E-10   |
| IRF1    | 0.06972539 | 4.27390138 | 0.36657996 |
| KRAS    | -0.0984657 | 6.37300321 | 0.05836715 |
| MIR17HG | -0.387355  | 2.09201982 | 0.00947483 |
| MYC     | -0.5061856 | 3.82646107 | 2.85E-08   |
| NRAS    | 0.23537769 | 8.01930555 | 8.24E-06   |
| STAT3   | 0.61536538 | 7.02728993 | 6.67E-43   |
| TEAD3   | -0.063607  | 3.5452602  | 0.52617761 |
| TGFB2   | 3.84396152 | 5.86523389 | 0          |
| TGFBR1  | 0.83375744 | 6.7094004  | 4.97E-70   |
| TP53    | 0.6827542  | 6.23630455 | 3.57E-39   |
| WWC1    | 0.54419528 | 5.63982615 | 5.83E-23   |

![image-20210506211615544](/Users/menghan/Library/Application Support/typora-user-images/image-20210506211615544.png![image-20210511221824096](/Users/menghan/Library/Application Support/typora-user-images/image-20210511221824096.png)



<font color = blue><font size= 8> 2021.05.07 </font></font>

## GSE50535

* Data Description

  ​																		*Differntial gene expression between pre- and post-treatment biopsies.* 

![image-20210507151138300](/Users/menghan/Library/Application Support/typora-user-images/image-20210507151138300.png)

**Error:** <u>SRR961664.fastq 数据损坏</u>

* Alignment and Quantification

```python
## RNA-Seq Pipeline 

import os
from subprocess import call

datadir = '/h/menghan/drug-resistance/RNA-seq/Public/GSE50535/raw' # 原始数据
data_list = os.listdir(datadir)
star_output = '/h/menghan/drug-resistance/RNA-seq/Public/GSE50535/star' 
rsem_output = '/h/menghan/drug-resistance/RNA-seq/Public/GSE50535/rsem'


data_list.sort()

DATA= {'patient1': data_list[0:2],
	'patient2':data_list[2:4],
	'patient3':data_list[4:6]}


for names , samples in DATA.items():
	## STAR
	samples.sort()
	fastq_path = [os.path.join(datadir,fastq) for fastq in samples]
	print('Starting STAR processing %s,%s'%(fastq_path[0],fastq_path[1]))
  ## Processing pre sample
	call('STAR --runThreadN 40 \
		--sjdbOverhang 50 \
		--outSAMtype BAM Unsorted \
		--quantMode TranscriptomeSAM GeneCounts \
		--readFilesCommand gunzip -c \
		--outFileNamePrefix %s_pre \
		--genomeDir /h/menghan/drug-resistance/RNA-seq/Public/GSE50535/index \
		--readFilesIn %s'%(os.path.join(star_output,names),fastq_path[0]),
		shell=True)
  ## Processing post sample
	call('STAR --runThreadN 40 \
		--sjdbOverhang 50 \
		--outSAMtype BAM Unsorted \
		--quantMode TranscriptomeSAM GeneCounts \
		--readFilesCommand gunzip -c \
		--outFileNamePrefix %s_post \
		--genomeDir /h/menghan/drug-resistance/RNA-seq/Public/GSE50535/index \
		--readFilesIn %s'%(os.path.join(star_output,names),fastq_path[1]),
		shell=True)

	## RSEM
	pre_output_path = os.path.join(star_output,'%s_preAligned.toTranscriptome.out.bam'%names)
	post_output_path = os.path.join(star_output,'%s_postAligned.toTranscriptome.out.bam'%names)
	call('rsem-calculate-expression -p 40 \
		--bam \
		--alignments \
		--append-names %s /h/menghan/ref/hg19/index/rsem/rsem %s_pre'%(pre_output_path,names),
		shell=True)

	call('rsem-calculate-expression -p 40 \
		--bam \
		--alignments \
		--append-names %s /h/menghan/ref/hg19/index/rsem/rsem %s_post'%(post_output_path,names),
		shell=True)
```

<font color = blue><font size= 8> 2021.05.09 </font></font>

*  Expression Matrix

```R
library(DESeq2)
library(readr)
library(tximport)
library(dplyr)
library(tidyr)
library(tibble)
library(apeglm)

rsemdir <- '/h/menghan/drug-resistance/RNA-seq/Public/GSE50535/rsem/'

specify_group <- '/h/menghan/drug-resistance/RNA-seq/Public/GSE50535/design.txt'

samples <- read_tsv(specify_group,col_names = T)

files <- file.path(rsemdir,paste0(samples$run,".genes.results"))

names(files) <- samples$sample

txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

write.table(txi$counts,'/h/menghan/drug-resistance/RNA-seq/Public/GSE50535/expression.txt',sep='\t')
```

| id                        | p2_pre | p2_post | p3_pre | p3_post |
| ------------------------- | ------ | ------- | ------ | ------- |
| ENSG00000215417.6_MIR17HG | 6      | 47      | 16.07  | 18      |
| ENSG00000141510.11_TP53   | 211.03 | 787.92  | 856.11 | 883.45  |
| ENSG00000112242.10_E2F3   | 113.77 | 536.68  | 422.99 | 376.18  |
| ENSG00000136997.10_MYC    | 232    | 202     | 453.98 | 184     |
| ENSG00000168610.10_STAT3  | 917    | 2904.2  | 1468   | 3715.2  |
| ENSG00000125347.9_IRF1    | 209.51 | 130.01  | 146.26 | 410     |

<font color = blue><font size= 8> 2021.05.11 </font></font>

## GSE110948

* Data Description

![image-20210511112804160](/Users/menghan/Library/Application Support/typora-user-images/image-20210511112804160.png)

* Alignment and Quantification

```python
## RNA-Seq Pipeline 

import os
from subprocess import call

datadir = '/h/menghan/drug-resistance/RNA-seq/Public/GSE110948/raw'
data_list = os.listdir(datadir)
star_output = '/h/menghan/drug-resistance/RNA-seq/Public/GSE110948/star'
rsem_output = '/h/menghan/drug-resistance/RNA-seq/Public/GSE110948/rsem'


data_list.sort()

DATA= {
	'A375': data_list[0],
	'A375R': data_list[1],
	'A375DR' : data_list[2],
	'patientA_pre': data_list[3],
	'patientA_post1': data_list[4],
	'patientA_post2': data_list[5],
	'patientB_pre': data_list[6],
	'patientB_post': data_list[7],
	'patientC_pre': data_list[8],
	'patientC_post': data_list[9],
	}


for names , samples in DATA.items():
	## STAR
	fastq_path = os.path.join(datadir,samples)
	print('Starting STAR processing %s : %s'%(names,fastq_path))
	call('STAR --runThreadN 40 \
		--sjdbOverhang 50 \
		--outSAMtype BAM Unsorted \
		--quantMode TranscriptomeSAM GeneCounts \
		--readFilesCommand gunzip -c \
		--outFileNamePrefix %s \
		--genomeDir /h/menghan/drug-resistance/RNA-seq/index \
		--readFilesIn %s'%(os.path.join(star_output,names),fastq_path),
		shell=True)


	# RSEM
	output_path = os.path.join(star_output,'%sAligned.toTranscriptome.out.bam'%names)
	call('rsem-calculate-expression -p 40 \
		--bam \
		--alignments \
		--append-names %s /h/menghan/ref/hg19/index/rsem/rsem %s'%(output_path,names),
		shell=True)

```

* Expression Matrix

  ```R
  rsemdir <- '/h/menghan/drug-resistance/RNA-seq/Public/GSE110948/rsem/'
  
  specify_group <- '/h/menghan/drug-resistance/RNA-seq/Public/GSE110948/design.txt'
  
  samples <- read_tsv(specify_group,col_names = T)
  
  files <- file.path(rsemdir,paste0(samples$run,".genes.results"))
  
  names(files) <- samples$sample
  
  txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
  
  write.table(txi$counts,'/h/menghan/drug-resistance/RNA-seq/Public/GSE110948/expression.txt',sep='\t')
  ```

  

| id                        | A375    | A375R   | A375DR  | patientA_pre | patientA_post1 | patientA_post2 | patientB_pre | patientB_post | patientC_pre | patientC_post |
| ------------------------- | ------- | ------- | ------- | ------------ | -------------- | -------------- | ------------ | ------------- | ------------ | ------------- |
| ENSG00000215417.6_MIR17HG | 90      | 118     | 87      | 48           | 20             | 29             | 110          | 29            | 46           | 20            |
| ENSG00000141510.11_TP53   | 741.42  | 1716.96 | 1088.89 | 573          | 718.96         | 1366.93        | 3835.77      | 1178.99       | 782.97       | 609           |
| ENSG00000136997.10_MYC    | 1369.18 | 1068    | 2650    | 634          | 1562           | 2200           | 2480         | 382           | 3681         | 2316          |
| ENSG00000112242.10_E2F3   | 903.44  | 2009    | 1694    | 379          | 992            | 699            | 1594         | 535           | 777          | 587           |
| ENSG00000168610.10_STAT3  | 1821.01 | 1857    | 3909    | 1610         | 3279           | 3046           | 8482.8       | 5057          | 5610         | 3390          |
| ENSG00000125347.9_IRF1    | 146     | 280     | 182     | 1771         | 1279           | 568            | 593          | 940           | 1441         | 1473          |

![image-20210511154139483](/Users/menghan/Library/Application Support/typora-user-images/image-20210511154139483.png)

<font color = blue><font size= 8> 2021.05.12 </font></font>

![image-20210512172746159](/Users/menghan/Library/Application Support/typora-user-images/image-20210512172746159.png)

![image-20210512172800437](/Users/menghan/Library/Application Support/typora-user-images/image-20210512172800437.png)



<font color = blue><font size= 8> 2021.05.20 </font></font>

## GSE65185

* Data Description

| sample          | Run        | mapki_treatment |
| --------------- | ---------- | --------------- |
| Pt1-baseline_1  | SRR1768821 | none            |
| Pt1-baseline_2  | SRR1768822 | none            |
| Pt1-DP1_1       | SRR1768823 | BRAFi           |
| Pt1-DP1_2       | SRR1768824 | BRAFi           |
| Pt2-baseline_1  | SRR1768825 | none            |
| Pt2-baseline_2  | SRR1768826 | none            |
| Pt2-DP1_1       | SRR1768827 | BRAFi           |
| Pt2-DP1_2       | SRR1768828 | BRAFi           |
| Pt2-DP2_1       | SRR1768829 | BRAFi           |
| Pt2-DP2_2       | SRR1768830 | BRAFi           |
| Pt3-baseline_1  | SRR1768831 | none            |
| Pt3-baseline_2  | SRR1768832 | none            |
| Pt3-DP1_1       | SRR1768833 | BRAFi           |
| Pt3-DP1_2       | SRR1768834 | BRAFi           |
| Pt3-DP1_3       | SRR1768835 | BRAFi           |
| Pt3-DP2_1       | SRR1768836 | BRAFi           |
| Pt3-DP2_2       | SRR1768837 | BRAFi           |
| Pt3-DP3_1       | SRR1768838 | BRAFi           |
| Pt3-DP3_2       | SRR1768839 | BRAFi           |
| Pt4-baseline_1  | SRR1768840 | none            |
| Pt4-baseline_2  | SRR1768841 | none            |
| Pt4-DP1_1       | SRR1768842 | BRAFi           |
| Pt4-DP1_2       | SRR1768843 | BRAFi           |
| Pt5-baseline_1  | SRR1768844 | none            |
| Pt5-baseline_2  | SRR1768845 | none            |
| Pt5-DP1_1       | SRR1768846 | BRAFi           |
| Pt5-DP1_2       | SRR1768847 | BRAFi           |
| Pt5-DP2_1       | SRR1768848 | BRAFi           |
| Pt5-DP2_2       | SRR1768849 | BRAFi           |
| Pt5-DP3_1       | SRR1768850 | BRAFi           |
| Pt5-DP3_2       | SRR1768851 | BRAFi           |
| Pt6-baseline_1  | SRR1768852 | none            |
| Pt6-baseline_2  | SRR1768853 | none            |
| Pt6-baseline_3  | SRR1768854 | none            |
| Pt6-DP1_1       | SRR1768855 | BRAFi           |
| Pt6-DP1_2       | SRR1768856 | BRAFi           |
| Pt6-DP1_3       | SRR1768857 | BRAFi           |
| Pt6-DP2_1       | SRR1768858 | BRAFi           |
| Pt6-DP2_2       | SRR1768859 | BRAFi           |
| Pt6-DP2_3       | SRR1768860 | BRAFi           |
| Pt8-baseline_1  | SRR1768861 | none            |
| Pt8-baseline_2  | SRR1768862 | none            |
| Pt8-DP1_1       | SRR1768863 | BRAFi           |
| Pt8-DP1_2       | SRR1768864 | BRAFi           |
| Pt8-DP2_1       | SRR1768865 | BRAFi           |
| Pt8-DP2_2       | SRR1768866 | BRAFi           |
| Pt8-DP3_1       | SRR1768867 | BRAFi           |
| Pt8-DP3_2       | SRR1768868 | BRAFi           |
| Pt9-baseline_1  | SRR1768869 | none            |
| Pt9-baseline_2  | SRR1768870 | none            |
| Pt9-DP1_1       | SRR1768871 | BRAFi           |
| Pt9-DP1_2       | SRR1768872 | BRAFi           |
| Pt9-DP2_1       | SRR1768873 | BRAFi           |
| Pt9-DP2_2       | SRR1768874 | BRAFi           |
| Pt10-baseline_1 | SRR1768875 | none            |
| Pt10-baseline_2 | SRR1768876 | none            |
| Pt10-DP1_1      | SRR1768877 | BRAFi           |
| Pt10-DP1_2      | SRR1768878 | BRAFi           |
| Pt10-DP2_1      | SRR1768879 | BRAFi           |
| Pt10-DP2_2      | SRR1768880 | BRAFi           |
| Pt10-DP3_1      | SRR1768881 | BRAFi           |
| Pt10-DP3_2      | SRR1768882 | BRAFi           |
| Pt10-DP4_1      | SRR1768883 | BRAFi           |
| Pt10-DP4_2      | SRR1768884 | BRAFi           |
| Pt10-DP5_1      | SRR1768885 | BRAFi           |
| Pt10-DP5_2      | SRR1768886 | BRAFi           |
| Pt10-DP6_1      | SRR1768887 | BRAFi           |
| Pt10-DP6_2      | SRR1768888 | BRAFi           |
| Pt10-DP7_1      | SRR1768889 | BRAFi           |
| Pt10-DP7_2      | SRR1768890 | BRAFi           |
| Pt10-DP8_1      | SRR1768891 | BRAFi           |
| Pt10-DP8_2      | SRR1768892 | BRAFi           |
| Pt10-DP9_1      | SRR1768893 | BRAFi           |
| Pt10-DP9_2      | SRR1768894 | BRAFi           |
| Pt15-baseline_1 | SRR1768895 | none            |
| Pt15-baseline_2 | SRR1768896 | none            |
| Pt15-DP1_1      | SRR1768897 | BRAFi           |
| Pt15-DP1_2      | SRR1768898 | BRAFi           |
| Pt15-DP2_1      | SRR1768899 | BRAFi           |
| Pt15-DP2_2      | SRR1768900 | BRAFi           |
| Pt16-baseline_1 | SRR1768901 | none            |
| Pt16-baseline_2 | SRR1768902 | none            |
| Pt16-baseline   | SRR1768903 | none            |
| Pt16-DP1_1      | SRR1768904 | BRAFi           |
| Pt16-DP1_2      | SRR1768905 | BRAFi           |
| Pt16-DP1_3      | SRR1768906 | BRAFi           |
| Pt17-baseline_1 | SRR1768907 | none            |
| Pt17-baseline_2 | SRR1768908 | none            |
| Pt17-baseline_3 | SRR1768909 | none            |
| Pt17-DP1_1      | SRR1768910 | BRAFi           |
| Pt17-DP1_2      | SRR1768911 | BRAFi           |
| Pt17-DP1_3      | SRR1768912 | BRAFi           |
| Pt17-DP2_1      | SRR1768913 | BRAFi           |
| Pt17-DP2_2      | SRR1768914 | BRAFi           |
| Pt17-DP2_3      | SRR1768915 | BRAFi           |
| Pt20-baseline_1 | SRR1768924 | none            |
| Pt20-baseline_2 | SRR1768925 | none            |
| Pt20-DP1_1      | SRR1768928 | BRAFi           |
| Pt20-DP1_2      | SRR1768929 | BRAFi           |

* Alignment and Quantification

```python
## RNA-Seq Pipeline 

import os
from subprocess import call
import pandas as pd

datadir = '/h/menghan/drug-resistance/RNA-seq/Public/GSE65185/raw'
data_list = os.listdir(datadir)
star_output = '/h/menghan/drug-resistance/RNA-seq/Public/GSE65185/star'
rsem_output = '/h/menghan/drug-resistance/RNA-seq/Public/GSE65185/rsem'
meta = pd.read_csv('/h/menghan/drug-resistance/RNA-seq/Public/GSE65185/meta2.csv')

samples = meta['Run'].tolist()

for sample  in samples:
	## STAR
	fastq1= os.path.join(datadir,'%s_1.fastq.gz'%sample)
	fastq2= os.path.join(datadir,'%s_2.fastq.gz'%sample)
	print('Starting STAR processing %s : %s & %s'%(sample,fastq1,fastq2))

	call('STAR --runThreadN 40 \
		--sjdbOverhang 99 \
		--outSAMtype BAM Unsorted \
		--quantMode TranscriptomeSAM GeneCounts \
		--readFilesCommand gunzip -c \
		--outFileNamePrefix %s \
		--genomeDir /h/menghan/drug-resistance/RNA-seq/index99 \
		--readFilesIn %s %s'%(os.path.join(star_output,sample),fastq1,fastq2),
		shell=True)

	# # RSEM
	output_path = os.path.join(star_output,'%sAligned.toTranscriptome.out.bam'%sample)
	call('rsem-calculate-expression -p 40 \
		--paired-end \
		--bam \
		--alignments \
		--append-names %s /h/menghan/ref/hg19/index/rsem/rsem %s'%(output_path,sample),
		shell=True)

```



* Expression Matrix













<font color = blue><font size= 8> 2021.06.05 </font></font>

## EGAD00001001306

* Alignment and Quantification

```python
## RNA-Seq Pipeline 

import os
from subprocess import call
import pandas as pd

datadir = '/h/menghan/drug-resistance/RNA-seq/Public/EGAD00001001306/raw'
star_output = '/h/menghan/drug-resistance/RNA-seq/Public/EGAD00001001306/star'
fastqdir = '/h/menghan/drug-resistance/RNA-seq/Public/EGAD00001001306/fq'
rsem_output = '/h/menghan/drug-resistance/RNA-seq/Public/EGAD00001001306/rsem'
meta = pd.read_csv('/h/menghan/drug-resistance/RNA-seq/Public/EGAD00001001306/meta.txt',sep='\t')



for idx, row  in meta.iterrows():

	
	FilePath = row['FilePath']
	FileName = row['FileName']
	File = os.path.join(datadir,FilePath,FileName)
	label = str(row['Patient']) + '_' + str(row['TimePoint'])
	print('Starting RSEM processing %s '%label)

	## Bam to Fastq
	call('samtools index -@ 40 %s '%(File),shell=True)
	call('bazam -bam %s -r1 %s_1.fq -r2 %s_2.fq '%(File,label,label),shell=True)


	## STAR
	fastq1 = os.path.join(fastqdir,'%s_1.fq'%label)
	fastq2 = os.path.join(fastqdir,'%s_2.fq'%label)
	call('STAR --runThreadN 40 \
		--sjdbOverhang 74 \
		--outSAMtype BAM Unsorted \
		--quantMode TranscriptomeSAM GeneCounts \
		--readFilesCommand gunzip -c \
		--outFileNamePrefix %s \
		--genomeDir /h/menghan/drug-resistance/RNA-seq/index74 \
		--readFilesIn %s %s'%(os.path.join(star_output,label),fastq1,fastq2),
		shell=True)


	## RSEM
	File = os.path.join(star_output,'%sAligned.toTranscriptome.out.bam'%label)

	call('rsem-calculate-expression -p 40 \
		--paired-end \
		--bam \
		--alignments \
		--append-names %s /h/menghan/ref/hg19/index/rsem/rsem %s'%(File,label),
		shell=True)
```







| Dataset   | DataType      | NumSample | Result             |
| --------- | ------------- | --------- | ------------------ |
| GSE50509  | Patients      | 5         | MIR17HG 下调       |
| GSE110948 | Patients      | 3         | MIR17HG 下调       |
| GSE65185  | Patients      | 10        | MIR17HG 变化不明显 |
| GSE50535  | Patients      | 2         | MIR17HG 上调       |
| GSE99923  | A375 cellline | -         | MIR17HG 下调       |



```sh
rsem-calculate-expression -p 40 \
		--bam \
		--alignments \
		--append-names sample_name \ 
		./rsem_index_dir/ \
		./star/*Aligned.toTranscriptome.out.bam
```

