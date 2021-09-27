# HiChIP Analysis

<font color = blue><font size= 8> 2021.05.13 </font></font>

## 1. Extract gene-enhancer Pairs (P-value < 0.05)

```python
import pandas as pd
import numpy as np
from subprocess import call
from functools import reduce
from matplotlib_venn import venn2
from matplotlib_venn import venn3

wt = pd.read_csv('./Wt_L4.interactions_FitHiC.bed',sep='\t')
wt = wt[wt['P-Value_Bias']<0.05]  

wt = wt.iloc[:,0:6]

wt['chr1'] = wt['chr1'].map(lambda x : x[3:])
wt['chr2'] = wt['chr2'].map(lambda x : x[3:])

wt.iloc[:,0:3].drop_duplicates().to_csv('./loops.binA.txt',index=False,sep='\t',header=None)
wt.iloc[:,3:6].drop_duplicates().to_csv('./loops.binB.txt',index=False,sep='\t',header=None)

call('bedtools intersect -a ./enhancer2.extend.nochr.bed -b ./loops.binA.txt -wb > loops.enhancer.binA.txt',shell=True)
call('bedtools intersect -a ./enhancer2.extend.nochr.bed -b ./loops.binB.txt -wb > loops.enhancer.binB.txt',shell=True)
call('bedtools intersect -a ./genes.info.txt -b ./loops.binA.txt -wb > loops.gene.binA.txt',shell=True)
call('bedtools intersect -a ./genes.info.txt -b ./loops.binB.txt -wb > loops.gene.binB.txt',shell=True)

enbinA = pd.read_csv('./loops.enhancer.binA.txt',sep='\t',header=None)
enbinB = pd.read_csv('./loops.enhancer.binB.txt',sep='\t',header=None)
enbinA.columns = ['enChr','enStart','enEnd','enhancer','hichipChr','hichipStart','hichipEnd']
enbinB.columns = ['enChr','enStart','enEnd','enhancer','hichipChr','hichipStart','hichipEnd']
enbinA = enbinA.drop_duplicates()
enbinB = enbinB.drop_duplicates()

genebinA = pd.read_csv('./loops.gene.binA.txt',sep='\t',header=None)
genebinB = pd.read_csv('./loops.gene.binB.txt',sep='\t',header=None)

genebinA = genebinA[[5,6,7,4]]
genebinB = genebinB[[5,6,7,4]]
genebinA.columns = ['chrA','startA','endA','gene']
genebinB.columns = ['chrB','startB','endB','gene']
genebinA = genebinA.drop_duplicates()
genebinB = genebinB.drop_duplicates()

wt.columns =['chrA','startA','endA','chrB','startB','endB']
tmp = pd.merge(wt , genebinA , on = ['chrA','startA','endA'] , how = 'left')
hichip_genes = pd.merge(tmp,genebinB , on = ['chrB','startB','endB'] , how = 'left')

tmp2 = pd.merge(hichip_genes , enbinA , left_on =['chrA','startA','endA'] ,right_on = ['hichipChr','hichipStart','hichipEnd'] ,how='left')
hichipEGpais = pd.merge(tmp2 , enbinB , left_on =['chrB','startB','endB'] ,right_on = ['hichipChr','hichipStart','hichipEnd'] ,how='left')

hichipEGpais = hichipEGpais[['chrA', 'startA', 'endA', 'chrB', 'startB', 'endB', 
                             'gene_x', 'gene_y',
                            'enhancer_x','enhancer_y']]

hichipEGpais = hichipEGpais.dropna(subset=['enhancer_x','enhancer_y'],how='all')

pairA = hichipEGpais[['gene_x','enhancer_y']]
pairB = hichipEGpais[['gene_y','enhancer_x']]
pairA = pairA.dropna(how='any')
pairB = pairB.dropna(how='any')

pairA.columns = ['gene','enhancer']
pairB.columns = ['gene','enhancer']
pairs = pd.concat([pairA,pairB])
pairs = pairs.drop_duplicates()
pairs = pairs[pairs['enhancer'] != '']

pairs.to_csv('./enhancer_gene.AllPairs.txt',sep='\t',index=False)

```

![image-20210519214634730](/Users/menghan/Library/Application Support/typora-user-images/image-20210519214634730.png)

## 2. Compare to Hi-C Results

* Associated with SV

![image-20210519214723276](/Users/menghan/Library/Application Support/typora-user-images/image-20210519214723276.png)

* Associated with SNV/Indel

![image-20210524103506413](/Users/menghan/Library/Application Support/typora-user-images/image-20210524103506413.png)

## 3. SV Outliers (*Fishhook* FDR < 0.5)

| Hotspot               | enhancer                                                     | gene                                                         | variant |
| --------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ | ------- |
| 10:3138000-3140000    | 10:3131938-3141173                                           | LOC101927824, PFKP                                           | TRA     |
| 15:63769000-63771000  | 15:63766838-63769439                                         | USP3, USP3-AS1, FBXL22                                       | TRA     |
| 19:53178000-53390000  | 19:53192766-53194794                                         | ZNF83, ZNF611, ZNF160                                        | DUP     |
| 1:150593000-150595000 | 1:150574207-150603862                                        | H2BC18, H3C13, H3P4, H3C15, H2AC18, H3C14, H2AC19, H2BC20P, H2AC20,  H2BC21, H2AC21, APH1A, C1orf54, RPRD2, MIR4257, ADAMTSL4, ADAMTSL4-AS1, MCL1,  CTSS, ARNT, SETDB1, ENSA, CTSK, CTXND2, MINDY1, GOLPH3L, HORMAD1, CERS2,  ANXA9, GABPB2, MLLT11 | TRA     |
| 6:56172000-56208000   | 6:56193775-56194743, 6:56198723-56202615                     | COL21A1, DST                                                 | DUP     |
| 7:131319000-131320000 | 7:131302439-131324264                                        | PODXL                                                        | TRA     |
| 8:19524000-19525000   | 8:19519734-19524167                                          | CSGALNACT1                                                   | TRA     |
| 9:21736000-22030000   | 9:21764067-21771931, 9:21827010-21839028, 9:21862923-21864570,  9:21917545-21920934 | MIR31HG, MTAP, IFNA8, CDKN2B-AS1                             | DEL     |



<font color = blue><font size= 8> 2021.05.24 </font></font>

## 4. Overlap with Our lib-gene

![image-20210524103804918](/Users/menghan/Library/Application Support/typora-user-images/image-20210524103804918.png)

![image-20210524114355875](/Users/menghan/Library/Application Support/typora-user-images/image-20210524114355875.png)

