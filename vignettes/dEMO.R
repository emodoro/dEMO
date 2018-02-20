## ----data sets read, message=FALSE---------------------------------------
library(airway)
data("airway")

## ----ExpressionSet, message=FALSE----------------------------------------
#Expression data set
expression<-assay(airway)
#phenotype datas
design<-colData(airway)
design<-design[,c("SampleName", "cell", "dex")]
design$GROUP<-ifelse(design$dex == "trt", 1, 0)

design$Sample<-c("A","B", "C", "D", "E", "F", "G", "H")
design$ID<-rownames(design)
rownames(design)<-design$Sample

phenodatas<-design

#change names of samples
colnames(expression)<-rownames(phenodatas)

#ExpressionSet
eset<-new("ExpressionSet", exprs=expression,annotation="hsa")
pData(eset)<-data.frame(phenodatas)


## ----Expression Zero, message=FALSE--------------------------------------
tabla0<-prop.table(table(rowSums(exprs(eset))==0))
tabla0

## ----filtered, message=FALSE, fig.align='center', cache=TRUE-------------
#At first, differences between library sizes are corrected, and then those genes whise expression profile is zero are removed. Two ExpressionSets are filtered, where the first one stores the raw data (eset) and the second one stores the corrected data (esetF) 
esetF <- eset
exprs(esetF) <- (t(t(exprs(esetF)) / (apply(exprs(esetF), 2, sum)))) * 10 ^ 6
esetF <- esetF[rowSums(exprs(eset)) > 0]
eset <- eset[rowSums(exprs(eset)) > 0]
#Quantiles are alculated and the minimum exression value A is fixed. In this case, the first quantile is selected as A value. 
A <- quantile(exprs(esetF), 0.5)
datos_condicion <- exprs(esetF) > A
#Filtred. The minimum number of K samples which have to show an expression greater or equal to A. K is the number of individuals of the smaller group tested.  
K=4
esetF<-esetF[apply(datos_condicion,1,sum)>=K]
eset<-eset[apply(datos_condicion,1,sum)>=K]
#histogram
par(mfrow=c(1,2))
x<-hist(assay(airway), breaks="Scott", main="a) Histogram \n Raw Datas", xlab = "counts")
hist(exprs(esetF), breaks="Scott", main="b) Histogram rlog \n dFiltered Datas",xlab = "counts)", ylim=c(0,max(x$mids)[1]))


## ----boxplot, message=FALSE, cache=TRUE, fig.align='center'--------------

boxplot(log(exprs(eset)+1,2), main="Boxplot log2", ylab="A", xlab="Sample",names=as.character(pData(eset)$Sample), col=as.numeric(pData(eset)$cell)) 


## ----Boxplot Norm, cache=TRUE, message=FALSE, fig.align='center'---------
#TMM 

library(edgeR)
par(mfrow=c(1,3))
nf<-calcNormFactors(exprs(eset))
esetTMM<-eset
exprs(esetTMM)<-t(t(exprs(esetTMM))*nf)

boxplot(log(exprs(eset)+1,2), main="Boxplot log2", ylab="A", xlab="Sample",names=as.character(pData(eset)$Sample), col=as.numeric(pData(eset)$cell)) 

boxplot(log(exprs(esetF)+0.1,2), main="a) \n Boxplot CMM", ylab="A", xlab="Sample",names=as.character(pData(eset)$Sample), col=as.numeric(pData(eset)$cell)) 

boxplot(log(exprs(esetTMM)+1,2), main="c) \n Boxplot TMM", ylab="A", xlab="Sample",names=as.character(pData(eset)$Sample), col=as.numeric(pData(eset)$cell)) 

## ----dEMO, message=FALSE-------------------------------------------------
library(dEMO)
library(edgeR)

expr<-exprs(eset)
conditions<-(pData(eset)[,"GROUP"])
TMMfac<-calcNormFactors.default(expr,method = "TMM")
exprT<-t(t(expr)*TMMfac)
exprorig<-t(t(assay(airway))*TMMfac)

testdEMOTMM<-dEMObuODlmTest(expr=exprT, condition = conditions,original = exprorig)

round(prop.table(table(testdEMOTMM$ADJ.PVAL<0.05))*100,2)

## ----edgeR, message=FALSE------------------------------------------------
#TMM
edgeR.dgelist <- DGEList(counts = expr, group = factor(conditions))
edgeR.dgelist <- calcNormFactors(edgeR.dgelist, method = "TMM")
edgeR.dgelist <- estimateCommonDisp(edgeR.dgelist)
edgeR.dgelist <- estimateTagwiseDisp(edgeR.dgelist, trend = "movingave")
edgeR.test <- exactTest(edgeR.dgelist)
edgeR.pvalues <- edgeR.test$table$PValue
edgeR.adjpvalues <- p.adjust(edgeR.pvalues, method = "BH")
testedgeRTMM <- data.frame(FC=edgeR.test$table$logFC, PVALUE=edgeR.test$table$PValue, ADJ.PVAL=edgeR.adjpvalues)

round(prop.table(table(testedgeRTMM$ADJ.PVAL<0.05))*100,2)

## ----Venn, message=FALSE-------------------------------------------------
#venn
library(gplots)
input <- list(edgeR = rownames(expr)[testedgeRTMM$ADJ.PVAL<=0.05],  dEMO=rownames(exprT)[testdEMOTMM$ADJ.PVAL<=0.05])

venn(input)

