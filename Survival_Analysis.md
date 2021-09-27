

# Survival Analysis of melanoma patients

<font color = blue><font size= 8> 2021.05.30 </font></font>

## GSE65185

| Patient | baseline    | treat      | LFC        | Dose of BRAF/MEK inhibitor (mg) | PFS  |
| ------- | ----------- | ---------- | ---------- | ------------------------------- | ---- |
| Pt6     | 1.111541333 | 0.22814108 | -2.2845634 | Vemu 960 bid                    | 132  |
| Pt20    | 1.78863     | 0.4154375  | -2.1061516 | Vemu 960 bid                    | 153  |
| Pt15    | 5.44073     | 1.62232125 | -1.7457407 | Vemu 960 bid                    | 373  |
| Pt1     | 8.433645    | 3.18179    | -1.4063177 | Vemu 960 bid                    | 212  |
| Pt22    | 0.7214725   | 0.284947   | -1.3402508 | Dabra 150mg bid/trame 2mg qd    | 60   |
| Pt17    | 0.567559    | 0.27357117 | -1.0528542 | Vemu 960 bid                    | 327  |
| Pt5     | 0.187667    | 0.10629148 | -0.820149  | Vemu 960 bid                    | 258  |
| Pt3     | 4.430295    | 2.80487067 | -0.6594685 | Vemu 960 bid                    | 137  |
| Pt4     | 1.39435     | 0.912606   | -0.6115287 | Vemu 960 bid                    | 196  |
| Pt9     | 7.224455    | 7.8662     | 0.12277811 | Vemu 960 bid                    | 97   |
| Pt19    | 0.3214345   | 0.3553572  | 0.14474515 | Vemu 960mg/GDC0973 60mg         | 562  |
| Pt16    | 0.395088    | 0.49623133 | 0.32883881 | Vemu 960 bid                    | 108  |
| Pt8     | 3.420945    | 4.53423167 | 0.4064632  | Vemu 960 bid                    | 238  |
| Pt23    | 0.1139718   | 0.189304   | 0.73202799 | Dabra 150 mg bid/trame 2 mg qd  | 300  |
| Pt24    | 1.2986      | 3.31229    | 1.35087188 | Dabra 150 mg bid/trame 2mg qd   | 365  |
| Pt2     | 0.8550825   | 2.327995   | 1.44495243 | Vemu 960 bid                    | 161  |
| Pt18    | 0.726275    | 2.69811    | 1.89336134 | Dabra 150mg/Trame 1mg           | 307  |
| Pt10    | 0.1044889   | 0.85217217 | 3.02779525 | Dabra 150 bid                   | 383  |

```R
mydata <- read.table('./RNA-seq/Public/GSE65185/survival_brafi_data.txt',header = 1,sep = '\t',stringsAsFactors = F)
mydata <- data.frame(mydata)


group <- ifelse(mydata$LFC < 0,'MIR17HG-down','MIR17HG-up')

sfit <- survfit(Surv(PFS)~group, data=mydata)

png('./RNA-seq/Public/GSE65185/mir17hg_brafi_pfs.png',width=800,height=600)
ggsurvplot(sfit, conf.int=F, pval=TRUE,font.size=12)
dev.off()
```

![image-20210608102126122](/Users/menghan/Library/Application Support/typora-user-images/image-20210608102126122.png)

## EGAD00001001306

<font color = blue><font size= 8> 2021.06.08 </font></font>

* Divede on_treatemnt & progression patient , compute log2foldchange

```python
data = pd.read_csv('./RNA-seq/Public/EGAD00001001306/expression.txt',sep='\t')
data['symbol'] = data['id'].map(lambda x: ''.join(x.split('_')[1:]))
data = data[data['symbol']=='MIR17HG']
clin = pd.read_csv('./RNA-seq/Public/EGAD00001001306/clin.txt',sep='\t')

clin_on = clin[(clin['TimePoint']=='on') |(clin['TimePoint']=='both')]
clin_pro = clin[(clin['TimePoint']=='pro') |(clin['TimePoint']=='both')]


def pre_sample(pt,df):
    pre = [p for p in data.columns[1:-1] if p == '%s_pre'%pt]

    return df.reset_index(drop=True).loc[0,pre][0]


def on_sample(pt,df):
    on = [p for p in data.columns[1:-1] if p == '%s_on'%pt]
    if len(on) > 0:
        return df.reset_index(drop=True).loc[0,on][0]
    else: 
        return 'lack'

def pro_sample(pt,df):
    pro = [p for p in data.columns[1:-2] if p == '%s_pro'%pt]
    if len(pro) > 0:
        return df.reset_index(drop=True).loc[0,pro][0]
    else: 
        return 'lack'
      
pstat_on = {}
for pt in clin_on['PT'].to_list():
    pstat_on[pt] = {}
    pstat_on[pt]['pre'] = pre_sample(pt, data)
    pstat_on[pt]['on'] = on_sample(pt,data)
    pstat_on[pt]['LFC'] = np.log2(on_sample(pt,data)+1) - np.log2(pre_sample(pt, data)+1)
    
pstat_pro = {}
for pt in clin_pro['PT'].to_list():
    pstat_pro[pt] = {}
    pstat_pro[pt]['pre'] = pre_sample(pt, data)
    pstat_pro[pt]['pro'] = pro_sample(pt,data)
    pstat_pro[pt]['LFC'] = np.log2(pro_sample(pt,data)+1) - np.log2(pre_sample(pt, data)+1)
    
pstat_on = pd.DataFrame.from_dict(pstat_on,orient='index')
pstat_pro = pd.DataFrame.from_dict(pstat_pro,orient='index')

sur_on = pd.merge(pstat_on,clin_on,left_index=True,right_on='PT',how='left')
sur_pro = pd.merge(pstat_pro,clin_pro,left_index=True,right_on='PT',how='left')

sur_on.to_csv('./RNA-seq/Public/EGAD00001001306/survival_data_on.txt',sep='\t',index=False)
sur_pro.to_csv('./RNA-seq/Public/EGAD00001001306/survival_data_pro.txt',sep='\t',index=False)
```

![image-20210608103354451](/Users/menghan/Library/Application Support/typora-user-images/image-20210608103354451.png)

![image-20210608103402366](/Users/menghan/Library/Application Support/typora-user-images/image-20210608103402366.png)

* Combine on_treatment patient and progression patient, filter patients sample of MIR17HG does not expressed

```python
sur_on = sur_on[sur_on['LFC']!=0]
sur_pro = sur_pro[sur_pro['LFC']!=0]

sur_all = pd.concat([sur_on,sur_pro])
sur_all.to_csv('./RNA-seq/Public/EGAD00001001306/survival_data_all.txt',sep='\t',index=False)
```

![image-20210608104238227](/Users/menghan/Library/Application Support/typora-user-images/image-20210608104238227.png)

* Survival Analysis of 3 group

1. *On_treatment VS Pre_treatment*

```R
library(dplyr)
library(survival)
library(ggplot2)
library(ggpubr)
library(magrittr)
library(survminer)

mydata <- read.table('./RNA-seq/Public/EGAD00001001306/survival_data_on.txt',header = 1,sep = '\t',stringsAsFactors = F)
mydata <- data.frame(mydata)

group <- ifelse(mydata$LFC > 0,'MIR17HG up','MIR17HG not up')

sfit <- survfit(Surv(PFS)~group, data=mydata)

png('./RNA-seq/Public/EGAD00001001306/mir17hg_pfs_on.png',width=800,height=600)
ggsurvplot(sfit, conf.int=F, pval=TRUE,font.size=12)
dev.off()
```

![image-20210608104114644](/Users/menghan/Library/Application Support/typora-user-images/image-20210608104114644.png)

2. *Progression VS pre_treatment*

```R
mydata <- read.table('./RNA-seq/Public/EGAD00001001306/survival_data_pro.txt',header = 1,sep = '\t',stringsAsFactors = F)
mydata <- data.frame(mydata)

group <- ifelse(mydata$LFC < 0,'MIR17HG down','MIR17HG up')

sfit <- survfit(Surv(PFS)~group, data=mydata)

png('./RNA-seq/Public/EGAD00001001306/mir17hg_pfs_pro.png',width=800,height=600)
ggsurvplot(sfit, conf.int=F, pval=TRUE,font.size=12)
dev.off()
```

![image-20210608104215052](/Users/menghan/Library/Application Support/typora-user-images/image-20210608104215052.png)

3. On/pro VS pre_treatment

```R
mydata <- read.table('./RNA-seq/Public/EGAD00001001306/survival_data_all.txt',header = 1,sep = '\t',stringsAsFactors = F)
mydata <- data.frame(mydata)

group <- ifelse(mydata$LFC < 0,'MIR17HG down','MIR17HG up')

sfit <- survfit(Surv(PFS)~group, data=mydata)

png('./RNA-seq/Public/EGAD00001001306/mir17hg_pfs_all.png',width=800,height=600)
ggsurvplot(sfit, conf.int=F, pval=TRUE,font.size=12)
dev.off()
```

​                

![image-20210608142113764](/Users/menghan/Library/Application Support/typora-user-images/image-20210608142113764.png)

## GSE50509

```R
mydata <- read.table('./RNA-seq/Public/GSE50509/survival_data.txt',header = 1,sep = '\t',stringsAsFactors = F)
mydata <- data.frame(mydata)

group <- ifelse(mydata$diff < 0,'MIR17HG-down','others')

sfit <- survfit(Surv(PFS)~group, data=mydata)

CairoPDF('./RNA-seq/Public/GSE50509/mir17hg_pfs.pdf',width=8,height=6,family='Arial')
ggsurvplot(sfit, conf.int=F, pval=TRUE,font.size=12)
dev.off()
```

![image-20210802114953810](/Users/menghan/Library/Application Support/typora-user-images/image-20210802114953810.png)

## Combine All Patients

<font color = blue><font size= 8> 2021.07.01 </font></font>

### Pre VS On-treatment

```R
mydata <- read.table('./RNA-seq/Public/survival_Pre_On.txt',header = 1,sep = '\t',stringsAsFactors = F)
mydata <- data.frame(mydata)

group <- ifelse(mydata$LFC < 0 ,'MIR17HG-down','others')

sfit <- survfit(Surv(PFS/7)~group, data=mydata)
CairoPDF('./RNA-seq/Public/mir17hg_pfs.pre_on.pdf',width=6,height=5,family='Arial')
ggsurvplot(sfit, conf.int=F, pval=TRUE,font.size=12,xlab=c('PFS(week)'))
dev.off()
```

![image-20210802115135894](/Users/menghan/Library/Application Support/typora-user-images/image-20210802115135894.png)

<font color = blue><font size= 8> 2021.07.20 </font></font>

### Pre VS Progression

```R
mydata <- read.table('./RNA-seq/Public/survival_Pre_Pro2.txt',header = 1,sep = '\t',stringsAsFactors = F)
mydata <- data.frame(mydata)
#mydata <- mydata[which(mydata$pre+mydata$pro > 10),]
group <- ifelse(mydata$Direction=='down' ,'MIR17HG Down','others')

sfit <- survfit(Surv(PFS/7)~group, data=mydata)
CairoPDF('./RNA-seq/Public/mir17hg_pfs.pre_pro.pdf',width=6,height=5,family='Arial')
ggsurvplot(sfit, conf.int=F, pval=TRUE,font.size=12, xlab=c('PFS(week)'))
dev.off()
```

![image-20210802115228897](/Users/menghan/Library/Application Support/typora-user-images/image-20210802115228897.png)

![image-20210802115320418](/Users/menghan/Library/Application Support/typora-user-images/image-20210802115320418.png)