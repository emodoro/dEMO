#' @title ROCs
#' @description ROCcurve is a functions which allow plotting ROCs for several hypothesis testing methods, shuch a several Differential Expression analysis.
#' @author Enrique Perez_Riesgo
#' @param dataRoc a data.frame where each column stores p values from each of the hypothesis test methods which are wanted to compare. The last column store the information related to real differential expression, where "0" means no differential expression and "1" means differential expression.
#' @param main the title of the plot
#' @return a plot with the ROC for every and each one of hypothesis test methods
#' @export ROCcurve
#' @examples
#' library(compcodeR)
#' library(edgeR)
#' set.seed(123456987)
#' datasCI<-generateSyntheticData(dataset = EMOresults, n.vars = 15000,samples.per.cond = 6, n.diffexp = 1500,repl.id = 1, seqdepth = 1e7,fraction.upregulated = 0.5,between.group.diffdisp = FALSE,filter.threshold.total = 1,filter.threshold.mediancpm = 0,fraction.non.overdispersed = 0)
#' expressiondata<-datasCI@count.matrix
#' TMMfac<-calcNormFactors.default(expressiondata,method = "TMM")
#' exprT<-t(t(expressiondata)*TMMfac)
#' conditions<-(datasCI@sample.annotations$condition-1)
#' testEMOTMM<-EMObuODlmTest(expr=exprT, condition = conditions,original = exprT,originalm=exprT)
#' edgeR.dgelist = DGEList(counts = expr, group = factor(conditions))
#' edgeR.dgelist = calcNormFactors(edgeR.dgelist, method = "TMM")
#' edgeR.dgelist = estimateCommonDisp(edgeR.dgelist)
#' edgeR.dgelist = estimateTagwiseDisp(edgeR.dgelist, trend = "movingave")
#' edgeR.test = exactTest(edgeR.dgelist)
#' edgeR.pvalues = edgeR.test$table$PValue
#' edgeR.adjpvalues = p.adjust(edgeR.pvalues, method = "BH")
#' testedgeRTMM<-data.frame(FC=edgeR.test$table$logFC, PVALUE=edgeR.test$table$PValue, ADJ.PVAL=edgeR.adjpvalues)
#' dataRoc<-data.frame(EMO=testEMOTMM$ADJ.PVAL,edgeR=testedgeRTMM$ADJ.PVAL, DE=expr3@variable.annotations$differential.expression)
#' ROCcurve(dataRoc = dataRoc, main="ROC")

ROCcurve<-function(dataRoc, main){
  secuencians<-sort(c(0.05,unique(dataRoc[,(dim(dataRoc)[2]-1)])),decreasing = FALSE)
  FPr<-data.frame(matrix(nrow = length(secuencians), ncol = (dim(dataRoc)[2]-1)))
  TPr<-data.frame(matrix(nrow = length(secuencians), ncol = (dim(dataRoc)[2]-1)))

  for(i in (1:length(secuencians))){
    FPr[i,]<-apply(dataRoc[,1:(dim(dataRoc)[2]-1)]<=(secuencians[i])&dataRoc[,"DE"]==0, 2,function(x){sum(x,na.rm = TRUE)})/sum(dataRoc[,"DE"]==0,na.rm = TRUE)
    TPr[i,]<-apply(dataRoc[,1:(dim(dataRoc)[2]-1)]<=(secuencians[i])&dataRoc[,"DE"]==1, 2,function(x){sum(x,na.rm = TRUE)})/sum(dataRoc[,"DE"]==1,na.rm = TRUE)
  }
  par(mar=c(8,7,5,0.5))
  plot(y = TPr[,1],x=FPr[,1], type="l", col=1, main=paste("ROCs\n",main),ylab="True Positive Rate", xlab="False Positive Rate", cex.axis=2, cex.lab=2)
  lines(y=c(0,1), x=c(0,1), type="p")
  legendID<-c(paste(colnames(dataRoc)[1]))
  legendcol<-c(1:(dim(dataRoc)[2]-1))
  for(i in 2:(dim(dataRoc)[2]-1)){
    lines(y = TPr[,i],x=FPr[,i], type="l", col=i)
    legendID<-c(legendID,paste(colnames(dataRoc)[i]))
  }
  legend("bottomright", col = legendcol, legend = legendID, lty=1)
}
