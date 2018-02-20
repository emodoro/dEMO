#' @title Confidence interval from log2(p1/p2) or log2(Fold Change)
#' @description This function provides a confidencial interval estimation from
#'   log2(proportions ratio) or log2(FC) where proportions belong to negative
#'   binomial population, such as RNA-seq data sets. Overdispersion estimation
#'   is carried out with \code{\link[edgeR]{estimateTagwiseDisp.default}} function.
#' @details \code{\link{EMObuODlmTest}} function is employed when data sets came
#'   from RNA-seq assays, where each row from expression matrix from RNA-seq
#'   data sets are dependent to the other rows into the same column (sample).
#'   Conversely, data sets store several rows where each of them does not came
#'   from the same experiment, that is, each row has nothing to do with the
#'   other rows, and \code{\link{dEMOOD1sci}} function has to be used.
#' @author Enrique Perez_Riesgo
#' @param expr data.frame, ExpressionSet or matrix after dataset has been
#'   filtered
#' @param condition Binary vector where 0 means control and 1 treatment
#' @param original data.frame, ExpressionSet or matrix before dataset has been
#'   filtered
#' @param alpha significance level for hypothesis test. Default value is 0.05
#' @param method.adj It is the multiple testing correction method. Default value
#'   correspond to Benjamini-Hochberg correction. c("holm", "hochberg",
#'   "hommel", "bonferroni", "BH", "BY","fdr", "none")
#' @return A data.frame is returned, which cointains log2FC, PVALUE, ADJ.PVAL
#'   (adjusted p value) and dEMO.STAT (statistic of dEMO test)
#' @export dEMOODci
#' @examples
#' datasCI <- compcodeR::generateSyntheticData(dataset = dEMOresults,
#' n.vars = 15000,samples.per.cond = 6, n.diffexp = 1500,repl.id = 1,
#' seqdepth = 1e7,fraction.upregulated = 0.5,between.group.diffdisp = FALSE,
#' filter.threshold.total = 1,filter.threshold.mediancpm = 0,
#' fraction.non.overdispersed = 0)
#' expressiondata <- datasCI@count.matrix
#' TMMfac <- edgeR::calcNormFactors.default(expressiondata, method = "TMM")
#' exprT <- t(t(expressiondata)*TMMfac)
#' conditions <- datasCI@sample.annotations$condition
#' testdEMOODci <- dEMOODci(expr = exprT, condition = (conditions-1),original = exprT)

dEMOODci <- function(expr, condition, original = NULL, alpha=0.05,
            method.adj="BH"){
  require(edgeR)
  if(class(expr) == "data.frame"){
    expr <- as.matrix(expr)
    genes <- rownames(expr)
  }
  if(class(expr) == "ExpressionSet"){
    expr <- exprs(expr)
    genes <- featureNames(expr)
  }
  if(class(expr) == "matrix"){
    genes <- dimnames(expr)[[1]]
  }
  if(!is.null(original)){
    if(class(original) == "ExpressionSet"){
      n <- apply(exprs(original), 2, sum)
      originalm <- exprs(original)
    }
    if(class(original) == "matrix"){
      n <- apply(original, 2, sum)
      originalm <- original
    }
    if(class(original) == "data.frame"){
      original <- as.matrix(original)
      n <- apply(original, 2, sum)
      originalm <- original
    }
  }else{
    n <- apply(expr, 2, sum)
  }

  #Overdispersion
  cd <- edgeR::estimateCommonDisp(y = originalm)
  cd1 <- edgeR::estimateCommonDisp(y = originalm[condition == 1])
  cd0 <- edgeR::estimateCommonDisp(y = originalm[condition == 0])
  cdg0 <- edgeR::estimateTagwiseDisp.default(y = expr[, condition == 0],
          dispersion = cd)
  cdg1 <- edgeR::estimateTagwiseDisp.default(y = expr[, condition == 1],
          dispersion = cd)
  cdg <- edgeR::estimateTagwiseDisp.default(y = expr, dispersion = cd)
  r1 <- 1 / cdg1
  r0 <- 1 / cdg0
  rt <- length(condition) / cdg
  rn1 <- 1 / cd1
  rn0 <- 1 / cd0
  rs <- matrix(c(rep(r0, sum(condition == 0)), rep(r1, sum(condition == 1))),
        ncol = length(condition), byrow = FALSE)

  #special cases
  expr[apply(expr[, condition == 0], 1, sum) == 0, 1] <- 1 / 2
  expr[apply(expr[, condition == 1], 1, sum) == 0, (sum(condition == 0) + 1)] <- 1 / 2

  probabilities <- 1 - t(t(expr) / n)
  variances <- t((n - t(expr))) * (1 - probabilities) / (probabilities ^ 2)
  n1 <- sum(n[condition == 1])
  x1 <- apply(expr[, condition == 1], 1, sum)
  A1 <- ((1 / x1) * (1 / log(2))) ^ 2
  uOD1 <- (apply(expr[, condition == 1], 1, sum)) *
          (1 + (apply(expr[, condition == 1], 1, sum)) * (1 / r1))
  variances1 <- apply(variances[, condition == 1], 1, sum)
  urep1 <- apply(expr[, condition == 1], 1, var)
  Fesc1 <- lm(c(1, sum(condition == 1)) ~ c(min(r1), max(r1)))
  Fesc1 <- r1 * Fesc1$coefficients[[2]] + Fesc1$coefficients[[1]]

  u1 <- ((uOD1) / sum(condition == 1)) + (urep1 + variances1) /
        (sum(condition == 1) / Fesc1)

  n0 <- sum(n[condition == 0])
  x0 <- apply(expr[, condition == 0], 1, sum)
  A0 <- ((1 / x0) * (1 / log(2))) ^ 2
  uOD0 <- (apply(expr[, condition == 0], 1, sum)) *
          (1 + (apply(expr[, condition == 0], 1, sum)) * (1 / r0))
  variances0 <- apply(variances[, condition == 0], 1, sum)
  urep0 <- apply(expr[, condition == 0], 1, var)
  Fesc0 <- lm(c(1, sum(condition == 0)) ~ c(min(r0), max(r0)))
  Fesc0 <- r0 * Fesc0$coefficients[[2]] + Fesc0$coefficients[[1]]
  u0 <- ((uOD0) / sum(condition == 0)) + (urep0 + variances0) /
        (sum(condition == 0) / Fesc0)

  FC <- log((x1 / (n1 + sum(condition == 1) * r1)) /
        (x0 / (n0 + sum(condition == 0) * r0)), 2)

  uc <- sqrt(A1 * u1 + A0 * u0)

  k <- qnorm(alpha / 2, lower.tail = FALSE)
  Uc <- k * uc

  CI <- data.frame(FCi = FC - Uc, FCs = FC + Uc)


  Conf.Int.SAM <- data.frame(log2FC = FC, LowBound = CI[,1],
                  UpBound = CI[, 2])
  rownames(Conf.Int.SAM) <- genes

  resultados <- Conf.Int.SAM
  return(resultados)
}
