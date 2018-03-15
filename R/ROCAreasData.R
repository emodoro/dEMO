#' @title ROCs AREA data
#' @description ROCAreasData is a functions which provide the values of Areas under ROCs, for both full interval and [0,0.05] interval of p values, for several hypothesis testing methods, shuch a several Differential Expression analysis.
#' @author Enrique Perez_Riesgo
#' @param dataRoc a data.frame where each column stores p values from each of the hypothesis test methods which are wanted to compare. The last column store the information related to real differential expression, where "0" means no differential expression and "1" means differential expression.
#' @param resample The number of resamplings of genes from data set in order to construct a boxplot from several samples from just one data set
#' @return a data frame with the AREA under ROC values, for both full interval and [0,0.05] interval of p values, for every and each one of hypothesis test methods
#' @export ROCAreaData
#' @examples
#' library(compcodeR)
#' library(edgeR)
#' set.seed(123456987)
#' datasCI<-generateSyntheticData(dataset = "EMOresults", n.vars = 15000,samples.per.cond = 6, n.diffexp = 1500,repl.id = 1, seqdepth = 1e7,fraction.upregulated = 0.5,between.group.diffdisp = FALSE,filter.threshold.total = 1,filter.threshold.mediancpm = 0,fraction.non.overdispersed = 0)
#' expressiondata<-datasCI@count.matrix
#' TMMfac<-calcNormFactors.default(expressiondata,method = "TMM")
#' exprT<-t(t(expressiondata)*TMMfac)
#' conditions<-(datasCI@sample.annotations$condition-1)
#' testEMOTMM<-EMObuODlmTest(expr=exprT, condition = conditions,original = exprT,originalm=exprT)
#' edgeR.dgelist = DGEList(counts = expressiondata, group = factor(conditions))
#' edgeR.dgelist = calcNormFactors(edgeR.dgelist, method = "TMM")
#' edgeR.dgelist = estimateCommonDisp(edgeR.dgelist)
#' edgeR.dgelist = estimateTagwiseDisp(edgeR.dgelist, trend = "movingave")
#' edgeR.test = exactTest(edgeR.dgelist)
#' edgeR.pvalues = edgeR.test$table$PValue
#' edgeR.adjpvalues = p.adjust(edgeR.pvalues, method = "BH")
#' testedgeRTMM<-data.frame(FC=edgeR.test$table$logFC, PVALUE=edgeR.test$table$PValue, ADJ.PVAL=edgeR.adjpvalues)
#' dataRoc<-data.frame(EMO=testEMOTMM$ADJ.PVAL,edgeR=testedgeRTMM$ADJ.PVAL, DE=datasCI@variable.annotations$differential.expression)
#' ROCAreasData(dataRoc = dataRoc, Area,resample = 10,DE = 0.1)


ROCAreasData<-function(dataRoc, resample=10, DE=0.1){
  AREAS<-data.frame(matrix(ncol=(dim(dataRoc)[2]-1), nrow=resample))
  colnames(AREAS)<-colnames(dataRoc)[1:(dim(dataRoc)[2]-1)]
  AREAS005<-AREAS
  for(j in 1:resample){
    indexNDE<-sample(rownames(dataRoc[dataRoc[,"DE"]==0,]), size = 1000*(1-DE), replace = FALSE)
    indexDE<-sample(rownames(dataRoc[dataRoc[,"DE"]==1,]), size = 1000*DE, replace = FALSE)
    secuencians<-sort(c(0.05,unique(dataRoc[c(indexNDE, indexDE),(dim(dataRoc)[2]-1)])),decreasing = FALSE)
    FPr<-data.frame(matrix(nrow = length(secuencians), ncol = (dim(dataRoc)[2]-1)))
    TPr<-data.frame(matrix(nrow = length(secuencians), ncol = (dim(dataRoc)[2]-1)))
    areas005<-data.frame(t(as.matrix(rep(0,(dim(dataRoc)[2]-1)))))
    areas<-data.frame(t(as.matrix(rep(0,(dim(dataRoc)[2]-1)))))

    for(i in (1:length(secuencians))){
      FPr[i,]<-apply(dataRoc[c(indexNDE,indexDE),1:(dim(dataRoc)[2]-1)]<=(secuencians[i])&dataRoc[c(indexNDE,indexDE),"DE"]==0, 2,function(x){sum(x,na.rm = TRUE)})/sum(dataRoc[c(indexNDE,indexDE),"DE"]==0,na.rm = TRUE)
      TPr[i,]<-apply(dataRoc[c(indexNDE,indexDE),1:(dim(dataRoc)[2]-1)]<=(secuencians[i])&dataRoc[c(indexNDE,indexDE),"DE"]==1, 2,function(x){sum(x,na.rm = TRUE)})/sum(dataRoc[c(indexNDE,indexDE),"DE"]==1,na.rm = TRUE)
    }

    for(k in 1:(dim(dataRoc)[2]-1)){
      areas005[,k]<- pracma::trapz(y = TPr[FPr[,k]<=0.05,k],x=FPr[FPr[,k]<=0.05,k])
      areas[,k]<- pracma::trapz(y = TPr[,k],x=FPr[,k])
    }
    AREAS005[j,]<-data.frame(areas005)
    AREAS[j,]<-data.frame(areas)
  }

  return(list(AREAS, AREAS005))

}
