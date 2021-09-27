<font color = Blue><font size = 12> 2021.04.04 </font></font>

# Input Data 

* DESeq2 LogFC (miRNA over-expressed)

* Target Scan Information

# Analysis

```R
library(TargetScore)

lfc <-mydata$log2FoldChange
cs <- mydata$Cumulative.weighted.context...score
pct <- mydata$Aggregate.PCT
ss <- cbind(cs,pct)
p.ts <- targetScore(lfc,ss,tol=1e-3)
mydata$targetScore <- c(p.ts)
```

<font color = Blue><font size = 12> 2021.04.05 </font></font>

# Hippo Target

Download: https://www.genome.jp/dbget-bin/www_bget?hsa04390

```
CRB2	 crumbs cell polarity complex component 2 [KO:K16681]
CRB1	 crumbs cell polarity complex component 1 [KO:K16681]
PARD3	 par-3 family cell polarity regulator [KO:K04237]
PARD6A	 par-6 family cell polarity regulator alpha [KO:K06093]
PARD6G	 par-6 family cell polarity regulator gamma [KO:K06093]
PARD6B	 par-6 family cell polarity regulator beta [KO:K06093]
PRKCZ	 protein kinase C zeta [KO:K18952] [EC:2.7.11.13]
PRKCI	 protein kinase C iota [KO:K06069] [EC:2.7.11.13]
PATJ	 PATJ crumbs cell polarity complex component [KO:K06092]
MPP5	 membrane palmitoylated protein 5 [KO:K06091]
AMOT	 angiomotin [KO:K16819]
YAP1	 Yes1 associated transcriptional regulator [KO:K16687]
WWTR1	 WW domain containing transcription regulator 1 [KO:K16820]
CDH1	 cadherin 1 [KO:K05689]
LIMD1	 LIM domain containing 1 [KO:K16682]
AJUBA	 ajuba LIM protein [KO:K16682]
WTIP	 WT1 interacting protein [KO:K16682]
NF2	 neurofibromin 2 [KO:K16684]
WWC1	 WW and C2 domain containing 1 [KO:K16685]
FRMD1	 FERM domain containing 1 [KO:K16821]
FRMD6	 FERM domain containing 6 [KO:K16822]
SAV1	 salvador family WW domain containing protein 1 [KO:K16686]
STK3	 serine/threonine kinase 3 [KO:K04412]
RASSF6	 Ras association domain family member 6 [KO:K09854]
RASSF1	 Ras association domain family member 1 [KO:K09850]
PPP2CA	 protein phosphatase 2 catalytic subunit alpha [KO:K04382] [EC:3.1.3.16]
PPP2CB	 protein phosphatase 2 catalytic subunit beta [KO:K04382] [EC:3.1.3.16]
PPP2R1B	 protein phosphatase 2 scaffold subunit Abeta [KO:K03456]
PPP2R1A	 protein phosphatase 2 scaffold subunit Aalpha [KO:K03456]
PPP2R2A	 protein phosphatase 2 regulatory subunit Balpha [KO:K04354]
PPP2R2B	 protein phosphatase 2 regulatory subunit Bbeta [KO:K04354]
PPP2R2C	 protein phosphatase 2 regulatory subunit Bgamma [KO:K04354]
PPP2R2D	 protein phosphatase 2 regulatory subunit Bdelta [KO:K04354]
LATS2	 large tumor suppressor kinase 2 [KO:K08791] [EC:2.7.11.1]
LATS1	 large tumor suppressor kinase 1 [KO:K08791] [EC:2.7.11.1]
MOB1A	 MOB kinase activator 1A [KO:K06685]
MOB1B	 MOB kinase activator 1B [KO:K06685]
PPP1CA	 protein phosphatase 1 catalytic subunit alpha [KO:K06269] [EC:3.1.3.16]
PPP1CB	 protein phosphatase 1 catalytic subunit beta [KO:K06269] [EC:3.1.3.16]
PPP1CC	 protein phosphatase 1 catalytic subunit gamma [KO:K06269] [EC:3.1.3.16]
TP53BP2	 tumor protein p53 binding protein 2 [KO:K16823]
LLGL2	 LLGL scribble cell polarity complex component 2 [KO:K06094]
LLGL1	 LLGL scribble cell polarity complex component 1 [KO:K06094]
SCRIB	 scribble planar cell polarity protein [KO:K16175]
DLG1	 discs large MAGUK scaffold protein 1 [KO:K12076]
DLG2	 discs large MAGUK scaffold protein 2 [KO:K12075]
DLG3	 discs large MAGUK scaffold protein 3 [KO:K21098]
DLG4	 discs large MAGUK scaffold protein 4 [KO:K11828]
DLG5	 discs large MAGUK scaffold protein 5 [KO:K24050]
CSNK1D	 casein kinase 1 delta [KO:K08959] [EC:2.7.11.1]
CSNK1E	 casein kinase 1 epsilon [KO:K08960] [EC:2.7.11.1]
TPTEP2-CSNK1E	 TPTEP2-CSNK1E readthrough [KO:K08960] [EC:2.7.11.1]
BTRC	 beta-transducin repeat containing E3 ubiquitin protein ligase [KO:K03362]
FBXW11	 F-box and WD repeat domain containing 11 [KO:K03362]
TP73	 tumor protein p73 [KO:K10148]
BBC3	 BCL2 binding component 3 [KO:K10132]
TEAD1	 TEA domain transcription factor 1 [KO:K09448]
TEAD4	 TEA domain transcription factor 4 [KO:K09448]
TEAD3	 TEA domain transcription factor 3 [KO:K09448]
TEAD2	 TEA domain transcription factor 2 [KO:K09448]
CCN2	 cellular communication network factor 2 [KO:K06827]
GLI2	 GLI family zinc finger 2 [KO:K16798]
AREG	 amphiregulin [KO:K09782]
BIRC5	 baculoviral IAP repeat containing 5 [KO:K08731]
AFP	 alpha fetoprotein [KO:K16144]
ITGB2	 integrin subunit beta 2 [KO:K06464]
FGF1	 fibroblast growth factor 1 [KO:K18496]
TGFB1	 transforming growth factor beta 1 [KO:K13375]
TGFB2	 transforming growth factor beta 2 [KO:K13376]
TGFB3	 transforming growth factor beta 3 [KO:K13377]
TGFBR1	 transforming growth factor beta receptor 1 [KO:K04674] [EC:2.7.11.30]
TGFBR2	 transforming growth factor beta receptor 2 [KO:K04388] [EC:2.7.11.30]
SMAD7	 SMAD family member 7 [KO:K19631]
SMAD2	 SMAD family member 2 [KO:K04500]
SMAD3	 SMAD family member 3 [KO:K23605]
SMAD4	 SMAD family member 4 [KO:K04501]
SERPINE1	 serpin family E member 1 [KO:K03982]
BMP2	 bone morphogenetic protein 2 [KO:K21283]
BMP4	 bone morphogenetic protein 4 [KO:K04662]
BMP5	 bone morphogenetic protein 5 [KO:K04663]
BMP6	 bone morphogenetic protein 6 [KO:K16620]
BMP7	 bone morphogenetic protein 7 [KO:K16621]
BMP8B	 bone morphogenetic protein 8b [KO:K16622]
BMP8A	 bone morphogenetic protein 8a [KO:K16622]
GDF5	 growth differentiation factor 5 [KO:K04664]
GDF6	 growth differentiation factor 6 [KO:K20012]
GDF7	 growth differentiation factor 7 [KO:K20013]
AMH	 anti-Mullerian hormone [KO:K04665]
BMPR1A	 bone morphogenetic protein receptor type 1A [KO:K04673] [EC:2.7.11.30]
BMPR1B	 bone morphogenetic protein receptor type 1B [KO:K13578] [EC:2.7.11.30]
BMPR2	 bone morphogenetic protein receptor type 2 [KO:K04671] [EC:2.7.11.30]
SMAD1	 SMAD family member 1 [KO:K04676]
ID1	 inhibitor of DNA binding 1, HLH protein [KO:K04680]
ID2	 inhibitor of DNA binding 2 [KO:K17693]
WNT1	 Wnt family member 1 [KO:K03209]
WNT2	 Wnt family member 2 [KO:K00182]
WNT2B	 Wnt family member 2B [KO:K00182]
WNT3	 Wnt family member 3 [KO:K00312]
WNT3A	 Wnt family member 3A [KO:K00312]
WNT4	 Wnt family member 4 [KO:K00408]
WNT5A	 Wnt family member 5A [KO:K00444]
WNT5B	 Wnt family member 5B [KO:K00444]
WNT6	 Wnt family member 6 [KO:K00445]
WNT7A	 Wnt family member 7A [KO:K00572]
WNT7B	 Wnt family member 7B [KO:K00572]
WNT8A	 Wnt family member 8A [KO:K00714]
WNT8B	 Wnt family member 8B [KO:K00714]
WNT9A	 Wnt family member 9A [KO:K01064]
WNT9B	 Wnt family member 9B [KO:K01064]
WNT10B	 Wnt family member 10B [KO:K01357]
WNT10A	 Wnt family member 10A [KO:K01357]
WNT11	 Wnt family member 11 [KO:K01384]
WNT16	 Wnt family member 16 [KO:K01558]
FZD1	 frizzled class receptor 1 [KO:K02432]
FZD7	 frizzled class receptor 7 [KO:K02432]
FZD2	 frizzled class receptor 2 [KO:K02235]
FZD3	 frizzled class receptor 3 [KO:K02329]
FZD4	 frizzled class receptor 4 [KO:K02354]
FZD5	 frizzled class receptor 5 [KO:K02375]
FZD8	 frizzled class receptor 8 [KO:K02375]
FZD6	 frizzled class receptor 6 [KO:K02376]
FZD10	 frizzled class receptor 10 [KO:K02842]
FZD9	 frizzled class receptor 9 [KO:K02842]
DVL3	 dishevelled segment polarity protein 3 [KO:K02353]
DVL2	 dishevelled segment polarity protein 2 [KO:K02353]
DVL1	 dishevelled segment polarity protein 1 [KO:K02353]
YWHAZ	 tyrosine 3-monooxygenase/tryptophan 5-monooxygenase activation protein zeta [KO:K16197]
YWHAB	 tyrosine 3-monooxygenase/tryptophan 5-monooxygenase activation protein beta [KO:K16197]
YWHAQ	 tyrosine 3-monooxygenase/tryptophan 5-monooxygenase activation protein theta [KO:K16197]
YWHAE	 tyrosine 3-monooxygenase/tryptophan 5-monooxygenase activation protein epsilon [KO:K06630]
YWHAH	 tyrosine 3-monooxygenase/tryptophan 5-monooxygenase activation protein eta [KO:K16198]
YWHAG	 tyrosine 3-monooxygenase/tryptophan 5-monooxygenase activation protein gamma [KO:K16198]
GSK3B	 glycogen synthase kinase 3 beta [KO:K03083] [EC:2.7.11.26]
CTNNB1	 catenin beta 1 [KO:K02105]
APC	 APC regulator of WNT signaling pathway [KO:K02085]
APC2	 APC regulator of WNT signaling pathway 2 [KO:K02085]
AXIN1	 axin 1 [KO:K02157]
AXIN2	 axin 2 [KO:K04385]
NKD1	 NKD inhibitor of WNT signaling pathway 1 [KO:K03213]
NKD2	 NKD inhibitor of WNT signaling pathway 2 [KO:K03213]
TCF7	 transcription factor 7 [KO:K02620]
TCF7L1	 transcription factor 7 like 1 [KO:K04490]
TCF7L2	 transcription factor 7 like 2 [KO:K04491]
LEF1	 lymphoid enhancer binding factor 1 [KO:K04492]
MYC	 MYC proto-oncogene, bHLH transcription factor [KO:K04377]
CCND1	 cyclin D1 [KO:K04503]
CCND2	 cyclin D2 [KO:K10151]
CCND3	 cyclin D3 [KO:K10152]
SOX2	 SRY-box transcription factor 2 [KO:K16796]
SNAI2	 snail family transcriptional repressor 2 [KO:K05706]
BIRC2	 baculoviral IAP repeat containing 2 [KO:K16060]
BIRC3	 baculoviral IAP repeat containing 3 [KO:K16060]
ACTG1	 actin gamma 1 [KO:K05692]
ACTB	 actin beta [KO:K05692]
CTNNA3	 catenin alpha 3 [KO:K05691]
CTNNA1	 catenin alpha 1 [KO:K05691]
CTNNA2	 catenin alpha 2 [KO:K05691]
```



## **Result**

* miR-18a-5p target in Hippo

| id                       | log2FoldChange_group1 | targetScore_group1 | log2FoldChange_group3 | targetScore_group3 | details                                                      |
| ------------------------ | --------------------- | ------------------ | --------------------- | ------------------ | ------------------------------------------------------------ |
| ENSG00000092969.7_TGFB2  | -5.1187494            | 0.32767766         | -1.0185056            | 0.2458069          | transforming growth factor beta 2  [KO:K13376]               |
| ENSG00000113645.9_WWC1   | -1.0227723            | 0.2469461          | -0.8960158            | 0.24079815         | WW and C2 domain containing 1  [KO:K16685]                   |
| ENSG00000106799.8_TGFBR1 | -0.4074718            | 0.20774407         | -0.8361008            | 0.23885901         | transforming growth factor beta  receptor 1 [KO:K04674] [EC:2.7.11.30] |
| ENSG00000111186.8_WNT5B  | -0.0426645            | 0.17019448         | -1.4864186            | 0.27079335         | Wnt family member 5B [KO:K00444]                             |
| ENSG00000073350.9_LLGL2  | -0.3317444            | 0.19916192         | -0.7630824            | 0.23651433         | LLGL scribble cell polarity  complex component 2 [KO:K06094] |
| ENSG00000105327.11_BBC3  | -0.8566662            | 0.24252322         | -0.0731898            | 0.17466424         | BCL2 binding component 3  [KO:K10132]                        |
| ENSG00000018408.10_WWTR1 | -0.1123128            | 0.17861599         | -0.7246143            | 0.23588735         | WW domain containing transcription  regulator 1 [KO:K16820]  |
| ENSG00000108379.5_WNT3   | -0.5537085            | 0.22670574         | -0.160383             | 0.18408353         | Wnt family member 3 [KO:K00312]                              |
| ENSG00000180340.5_FZD2   | -0.7429185            | 0.2389287          | -0.0420326            | 0.17165861         | frizzled class receptor 2  [KO:K02235]                       |
| ENSG00000138696.6_BMPR1B | 0.40724357            | 0.1318733          | -1.3809264            | 0.26536482         | bone morphogenetic protein  receptor type 1B [KO:K13578] [EC:2.7.11.30] |

* miR-18a-3p target in Hippo

| id                       | symbol | log2FoldChange_group1 | targetScore_group1 | log2FoldChange_group3 | targetScore_group3 | details                                                      |
| ------------------------ | ------ | --------------------- | ------------------ | --------------------- | ------------------ | ------------------------------------------------------------ |
| ENSG00000180340.5_FZD2   | FZD2   | -0.7429185            | 0.33885161         | -0.0420326            | 0.46744324         | frizzled class receptor 2  [KO:K02235]                       |
| ENSG00000007866.14_TEAD3 | TEAD3  | -0.23847              | 0.28789535         | -0.1241151            | 0.48608651         | TEA domain transcription factor 3  [KO:K09448]               |
| ENSG00000213923.6_CSNK1E | CSNK1E | 0.12028671            | 0.24074666         | -0.2538614            | 0.51004672         | casein kinase 1 epsilon  [KO:K08960] [EC:2.7.11.1]           |
| ENSG00000188763.3_FZD9   | FZD9   | -0.2038459            | 0.2845078          | 0.00573943            | 0.45565445         | frizzled class receptor 9  [KO:K02842]                       |
| ENSG00000107404.13_DVL1  | DVL1   | 0.03283499            | 0.25415675         | -0.1004991            | 0.48095177         | dishevelled segment polarity  protein 1 [KO:K02353]          |
| ENSG00000204217.8_BMPR2  | BMPR2  | -0.1086116            | 0.27393227         | -0.0115453            | 0.45999163         | bone morphogenetic protein  receptor type 2 [KO:K04671] [EC:2.7.11.30] |
| ENSG00000074047.16_GLI2  | GLI2   | 0.34146536            | 0.20865354         | -0.1209201            | 0.48540395         | GLI family zinc finger 2  [KO:K16798]                        |
| ENSG00000101665.4_SMAD7  | SMAD7  | 0.12861559            | 0.2393347          | 0.07494456            | 0.43753692         | SMAD family member 7 [KO:K19631]                             |
| ENSG00000163558.8_PRKCI  | PRKCI  | 0.50212597            | 0.18852949         | -0.1310934            | 0.48756178         | protein kinase C iota [KO:K06069]  [EC:2.7.11.13]            |
| ENSG00000108379.5_WNT3   | WNT3   | -0.5537085            | 0.17695485         | -0.160383             | 0.48851173         | Wnt family member 3 [KO:K00312]                              |

![image-20210430170148430](/Users/menghan/Library/Application Support/typora-user-images/image-20210430170148430.png)

# PI3K-AKT Target

<font color = blue><font size= 8> 2021.06.06 </font></font>

* PI3K-AKT gene list

![image-20210603175510080](/Users/menghan/Library/Application Support/typora-user-images/image-20210603175510080.png)

```python
pi3k = pd.read_csv('./TargetScore/PI3K.genelist.txt',sep='\t')

p = '3p' # 5p or 3p

g1 = pd.read_csv('./TargetScore/%s/group1.TSresults.csv'%p)

g1_pi3k = pd.merge(g1,pi3k,on='symbol',how='inner')

g1_pi3k.to_csv('./TargetScore/%s/group1.PI3K.txt'%p,sep='\t',index=False)



g3 = pd.read_csv('./TargetScore/%s/group3.TSresults.csv'%p)

g3_pi3k = pd.merge(g3,pi3k,on='symbol',how='inner')

g3_pi3k.to_csv('./TargetScore/%s/group3.PI3K.txt'%p,sep='\t',index=False)

overlap = pd.merge(g1_pi3k,g3_pi3k,on='symbol',how='inner')

overlap.to_csv('./TargetScore/%s/overlap.PI3K.txt'%p,sep='\t',index=False)
```

* Result

![image-20210603222309515](/Users/menghan/Library/Application Support/typora-user-images/image-20210603222309515.png)

