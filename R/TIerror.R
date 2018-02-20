#' @title Type I error
#' @description TIerror is a functions which allow plotting the type I error commited by several hypothesis testing methods, shuch a several Differential Expression analysis.
#' @author Enrique Perez_Riesgo
#' @param dataRoc a data.frame where each column stores p values from each of the hypothesis test methods which are wanted to compare. The last column store the information related to real differential expression, where "0" means no differential expression and "1" means differential expression.
#' @param resample The number of resamplings of genes from data set in order to construct a boxplot from several samples from just one data set. Default value is 100.
#' @param main the title of the plot
#' @return a box plot with the Type I error for every and each one of hypothesis test methods
#' @export TIerror
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
#' TIerror(dataRoc = dataRoc, resample=100, main="Type I Error")

TIerror<-function(dataRoc, resample=100,main){
  TIe<-data.frame(matrix(ncol=(dim(dataRoc)[2]-1), nrow=resample))
  colnames(TIe)<-colnames(dataRoc)[1:(dim(dataRoc)[2]-1)]
  for(i in 1:resample){
    genessample<-sample(rownames(dataRoc),1000)
    TIe[i,]<-data.frame(t(as.matrix((apply(dataRoc[genessample,1:(dim(dataRoc)[2]-1)]<=0.05, 2,function(x){sum(x,na.rm = TRUE)}))/1000)))
  }

  par(mar=c(8,7,5,1))
  boxplot(TIe, main=paste("Type I error Rate\n",main),ylab="Typer I Error Rate\n", ylim=c(0,0.2), cex.axis=2, las=2, cex.lab=2)

  lines(x =c(0:(dim(dataRoc)[2])),y=c(rep(0.05,(dim(dataRoc)[2]+1))), type="l", col=2,cex.lab=2)

}
