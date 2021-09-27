

<font color = blue><font size= 8> 2021.09.03 </font></font>

# 1. CellRanger process

* Lane 1 :

```shell
cellranger count \
--id=A375VR_lane1 \
--fastqs=/f/menghan/Single-cell/data/210624_TJJ_A375KO_Cropseq/ \
--sample=P_1-3,P_1-4,P_1-5,P_1-6 \
--transcriptome=/f/menghan/ref/cellranger/refdata-gex-GRCh38-2020-A \
--nosecondary --localcores=30 --localmem=30 \
--lanes=1 --chemistry=SC3Pv3
```

* Lane 4:

```shell
cellranger count \
--id=A375VR_lane4 \
--fastqs=/f/menghan/Single-cell/data/210624_TJJ_A375KO_Cropseq/ \
--sample=P_1-1,P_1-2,P_1-7,P_1-8 \
--transcriptome=/f/menghan/ref/cellranger/refdata-gex-GRCh38-2020-A \
--nosecondary --localcores=30 --localmem=30 \
--lanes=4 --chemistry=SC3Pv3
```

```
cellranger count \
--id=A375VR_all \
--fastqs=/f/menghan/Single-cell/data/210624_TJJ_A375KO_Cropseq/ \
--sample=P_1-1,P_1-2,P_1-3,P_1-4,P_1-5,P_1-6,P_1-7,P_1-8 \
--transcriptome=/f/menghan/ref/cellranger/refdata-gex-GRCh38-2020-A \
--nosecondary --localcores=30 --localmem=30 \
--lanes=1,4 --chemistry=SC3Pv3
```



* cDNA:

```
cellranger count \
--id=A375VR_cDNA \
--fastqs=/f/menghan/Single-cell/data/210818_TJJ_A375VR_Cropseq_cDNA/ \
--sample=CROP_seq_A375-1 \
--transcriptome=/f/menghan/ref/cellranger/refdata-gex-GRCh38-2020-A \
--nosecondary --localcores=30 --localmem=30 \
--lanes=2 --chemistry=SC3Pv3
```



* Output example:

![image-20210831162648739](/Users/menghan/Library/Application Support/typora-user-images/image-20210831162648739.png)



# 2. Quality control

<font color = blue><font size= 8> 2021.09.04 </font></font>

## Get barcode from cDNA enrichment

```shell
python ./single-cell-ko-screens/get_barcodes.py \
--input_bams ./batch2/A375VR_cDNA/outs/possorted_genome_bam.bam \
-o ko_barcodes.txt \
--whitelist whitelist.txt \
--search_seq GTGGAAAGGACGAAACACCG \
--force_correction 2
```

```sh
python ./single-cell-ko-screens/get_barcodes.py \
--input_bams ./batch2/A375VR_cDNA/outs/possorted_genome_bam.bam \
-o ko_barcodes.no_wl.txt \
--search_seq GTGGAAAGGACGAAACACCG \
--force_correction 3
```



```shell
Rscript ./single-cell-ko-screens/preprocess_cfg.R sample_metadata.txt output_cds.rds output_metadata.txt --guide_metadata gene_barcode_associations.txt --barcode_enrichment_qc_plot qc_plot.png
```

![image-20210917143944412](/Users/menghan/Library/Application Support/typora-user-images/image-20210917143944412.png)

<font color =red> 只进入1条sgRNA序列 </font>

