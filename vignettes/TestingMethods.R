## ------------------------------------------------------------------------
library(edgeR)
library(compcodeR)
library(DESeq2)

samples<-6
OD<-"OverDisperssion6"
fractionNonOD<-0
ODdif<-FALSE
DEperc<-0.3
genesnum<-15000
upreg<-0.5
librarysiz<-1e6

name_results<-paste(fractionNonOD,OD,ODdif,"_",as.character(samples),"_",as.character(DEperc),"_", as.character(upreg),"_",as.character(librarysiz), sep="")
main<-name_results

set.seed(123456987)
expr3<-generateSyntheticData(dataset = name_results, n.vars = genesnum,samples.per.cond = samples, n.diffexp = DEperc*genesnum,repl.id = 1, seqdepth = librarysiz,fraction.upregulated = upreg,between.group.diffdisp = ODdif,filter.threshold.total = 1,filter.threshold.mediancpm = 0,fraction.non.overdispersed = fractionNonOD)

#Extract information of interest from simulated RNA-seq data sets
expr<-expr3@count.matrix
conditions<-(expr3@sample.annotations$condition-1)


## ----dEMO test, message=FALSE, warning=FALSE-----------------------------
#TMM normalization
TMMfac<-calcNormFactors.default(expr,method = "TMM")
exprT<-t(t(expr)*TMMfac)

#dEMO test
library(dEMO)
testdEMOTMM<-dEMObuODlmTest(expr=exprT, condition = conditions,original = exprT)

## ----edgeR, message = FALSE, warning = FALSE-----------------------------
library(edgeR)
edgeR.dgelist = DGEList(counts = expr, group = factor(conditions))
edgeR.dgelist = calcNormFactors(edgeR.dgelist, method = "TMM")
edgeR.dgelist = estimateCommonDisp(edgeR.dgelist)
edgeR.dgelist = estimateTagwiseDisp(edgeR.dgelist, trend = "movingave")
edgeR.test = exactTest(edgeR.dgelist)
edgeR.pvalues = edgeR.test$table$PValue
edgeR.adjpvalues = p.adjust(edgeR.pvalues, method = "BH")
testedgeRTMM<-data.frame(FC=edgeR.test$table$logFC, PVALUE=edgeR.test$table$PValue, ADJ.PVAL=edgeR.adjpvalues)


## ----DESeq, message = FALSE, warning = FALSE-----------------------------
library(DESeq)
DESeq.cds = newCountDataSet(countData = expr,conditions = factor(conditions))
DESeq.cds = estimateSizeFactors(DESeq.cds)
DESeq.cds = estimateDispersions(DESeq.cds, sharingMode = "maximum",method = "pooled", fitType = "local")
DESeq.test = nbinomTest(DESeq.cds, "0", "1")
DESeq.pvalues = DESeq.test$pval
DESeq.adjpvalues = p.adjust(DESeq.pvalues, method = "BH")
testDESeq<-data.frame(FC=DESeq.test$log2FoldChange, PVALUE=DESeq.test$pval, ADJ.PVAL=DESeq.adjpvalues)


## ----NBPseq, message = FALSE, warning = FALSE----------------------------
library(edgeR)
library(NBPSeq)
NBPSeq.dgelist = DGEList(counts = expr, group = factor(conditions))
NBPSeq.dgelist = calcNormFactors(NBPSeq.dgelist, method = "TMM")
NBPSeq.norm.factors = as.vector(NBPSeq.dgelist$samples$norm.factors)
NBPSeq.test = nbp.test(counts = expr, grp.ids = conditions,grp1 = 0, grp2 = 1, norm.factors = NBPSeq.norm.factors)
NBPSeq.pvalues = NBPSeq.test$p.values
NBPSeq.adjpvalues = NBPSeq.test$q.values
testNBPseq<-data.frame(FC=NBPSeq.test$log.fc, PVALUE=NBPSeq.test$p.values, ADJ.PVAL=NBPSeq.test$q.values)

## ----pvalues stored, message = FALSE, warning = FALSE--------------------
dataRoc<-data.frame(EMO=testdEMOTMM$ADJ.PVAL,NBPseq=testNBPseq$ADJ.PVAL,edgeR=testedgeRTMM$ADJ.PVAL, DESeq=testDESeq$ADJ.PVAL, DE=expr3@variable.annotations$differential.expression)

dataRocp<-data.frame(EMO=testdEMOTMM$PVALUE,NBPseq=testNBPseq$PVALUE,edgeR=testedgeRTMM$PVALUE, DESeq=testDESeq$PVALUE, DE=expr3@variable.annotations$differential.expression)
dataRocp<-dataRocp[expr3@variable.annotations$differential.expression==0,]

## ----ROCs, message = FALSE, warning = FALSE------------------------------
ROCcurve(dataRoc = dataRoc, main="ROCs")

ROCAreas(dataRoc = dataRoc, main="ROCs Areas",resample = 10,DE = DEperc)
areasdata<-ROCAreasData(dataRoc = dataRoc, resample = 10,DE = DEperc)

ROCAreas005(dataRoc = dataRoc, main="ROCs Areas 0.05",resample = 10,DE = DEperc)

## ----FDR curves, message= FALSE, warning= FALSE--------------------------
FDRcurve(dataRoc = dataRoc, main="")

