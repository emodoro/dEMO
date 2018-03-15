#' @title FDR curves
#' @description FDR curves is a functions which allow plotting the Real False Discoveries commited as number of differential expressed genes are identified by several hypothesis testing methods, shuch a several Differential Expression analysis.
#' @author Enrique Perez_Riesgo
#' @param dataRoc a data.frame where each column stores p values from each of the hypothesis test methods which are wanted to compare. The last column store the information related to real differential expression, where "0" means no differential expression and "1" means differential expression.
#' @param main the title of the plot
#' @return a plot with the FDRs curves for every and each one of hypothesis test methods
#' @export FDRcurve
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
#' FDRcurve(dataRoc = dataRoc, main="FDR curves")

FDRcurve<-function(dataRoc, main){
  secuencians<-sort(c(0.05,unique(dataRoc[,(dim(dataRoc)[2]-1)])),decreasing = FALSE)
  FP<-data.frame(matrix(nrow = length(secuencians), ncol = (dim(dataRoc)[2]-1)))
  DEGs<-data.frame(matrix(nrow = length(secuencians), ncol = (dim(dataRoc)[2]-1)))
  for(i in (1:length(secuencians))){
    FP[i,]<-apply(dataRoc[,1:(dim(dataRoc)[2]-1)]<=(secuencians[i])&dataRoc[,"DE"]==0, 2,function(x){sum(x,na.rm = TRUE)})
    DEGs[i,]<-apply(dataRoc[,1:(dim(dataRoc)[2]-1)]<=(secuencians[i]),2,function(x){sum(x,na.rm = TRUE)})
  }
  par(mar=c(8,7,5,0.5))
  plot(x = DEGs[,1],y=FP[,1], type="l", col=1, main=paste("FDR curve\n",main),ylab="False Discovery Rate\n", xlab="DEGs", xlim=c(1,1000), ylim=c(1,1000), log = "y", cex.axis=2, cex.lab=2)
  legendID<-c(paste(colnames(dataRoc)[1]))
  legendcol<-c(1:(dim(dataRoc)[2]-1))
  for(i in 2:(dim(dataRoc)[2]-1)){
    lines(x = DEGs[,i],y=FP[,i], type="l", col=i)
    legendID<-c(legendID,paste(colnames(dataRoc)[i]))
  }
  legend("bottomright", col = legendcol, legend = legendID, lty=1)
  #return(list(FP=FP,DEGs=DEGs))
}
