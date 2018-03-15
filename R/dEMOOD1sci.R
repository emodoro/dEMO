#' @title Regresion coefficients for estimate the Correction Factors for
#'   Overdispersion Estiations
#' @description This function provides those regression coefficients necesaries
#'   for estimate the correction factors for Overdispersion Estimations, which
#'   is used as argument of EMOOD1s function, of data sets by considering each
#'   row as a different experiment, where each row has nothing to do with the
#'   others rows.
#' @details Sometimes, data sets store several rows where each of them does not
#'   came from the same experiment, that is, each row has nothing to do with the
#'   other rows. Conversely, the rows from RNA-seq data sets are dependent to
#'   the other rows into the same column (sample), so other function has to be
#'   apply here, such as EMObuODlmTest().
#' @author Enrique Perez_Riesgo
#' @family One sample functions
#' @seealso \code{\link{dEMOOD1sci}} for use dEMO::Fesc outcomes as argument for
#'   this function and test the differential expression for independents rows.
#' @param data matrix where each row is a different experiment which has nothing
#'   to do with the other rows or experiments.
#' @param librarysize each sample from a row belongs to a negative binomial
#'   distribution with a "n1" or library size, and j is the number of samples
#'   which belong to the same condition that sample 1, the column "j+1" is n1,
#'   column "j+2" is n2, and so on.
#' @param condition Binary vector where 0 means control and 1 treatment
#' @return A data.frame is returned, which cointains tow columns, the first one
#'   with the correction factor from samples which belong to condition "0" and
#'   the second column with the correction fators from samples which belong to
#'   condition "1".
#' @export Fesc
#' @examples
#' #observations
#' obsev <- 100
#' #population parameters
#' p1 <- 0.2
#' p2 <- 0.5
#' r1 <- 10^7
#' r2 <- 10^7
#' #simulations
#' x21 <- rnbinom(prob = p2, obsev, size = r2)
#' n21 <- x21 + r2
#' x22 <- rnbinom(prob = p2, obsev, size = r2)
#' n22 <- x22 + r2
#' x23 <- rnbinom(prob = p2, obsev, size = r2)
#' n23 <- x23 + r2
#' q2p <- x21 + x22 + x23
#' n2s <- n21 + n22 + n23
#' x11 <- rnbinom(prob = p1, obsev, size = r1)
#' n11 <- x11 + r1
#' x12 <- rnbinom(prob = p1, obsev, size = r1)
#' n12 <- x12 + r1
#' x13 <- rnbinom(prob = p1, obsev, size = r1)
#' n13 <- x13 + r1
#' q1p <- x11 + x12 + x13
#' n1s <- n11 + n12 + n13
#' rp1 <- n1s - q1p
#' rp2 <- n2s - q2p
#' #data set with all above simulations
#' expressiondata <- matrix(c(x21, x22, x23, x11, x12, x13, n21, n22,
#' n23, n11, n12, n13), ncol = 12, byrow = FALSE)
#' #correction factors
#' Fescoefs <- Fesc( expressiondata[, 1:6], expressiondata[, 7:12],
#' condition = c(0, 0, 0, 1, 1, 1))


Fesc <- function(data, librarysize, condition){
  cd <- edgeR::estimateCommonDisp( y = data)
  cdg0 <- edgeR::estimateTagwiseDisp.default(y = data[, condition==0],
          dispersion = cd)
  cdg1 <- edgeR::estimateTagwiseDisp.default(y = data[,condition==1],
          dispersion = cd)
  cdg <- edgeR::estimateTagwiseDisp.default(y = data, dispersion = cd)
  r1 <- 1/cdg1
  r0 <- 1/cdg0
  rt <- length(condition)/cdg
  Fesc0 <- lm(c(1, sum(condition == 0))~c(min(r0), max(r0)))
  Fesc1 <- lm(c(1, sum(condition == 1))~c(min(r1), max(r1)))
  Fescoef <- data.frame(Fesc0 = c(Fesc0$coefficients[[1]],
            Fesc0$coefficients[[2]]), Fesc1 = c(Fesc1$coefficients[[1]],
            Fesc1$coefficients[[2]]))
  return(Fescoef)
}


#' @title Confidence interval from log2(p1/p2), feature by feature,
#' @description This function provides a confidencial interval from
#'   log2(proportions ratio) or log2(FC) where proportions belong to negative
#'   binomial population and from data sets by considering each row as a
#'   different experiment, from data sets by considering each row as a different
#'   experiment, where each row has nothing to do with the others rows.
#'   Overdispersion estimation is carried out with
#'   edgeR::estimateTagwiseDisp.default function. dEMO::EMOOD1sci function works
#'   row by row, so if you want to calculate dEMO confidence intervals when rows
#'   are independent experiments, a loop "for" is recomended, as is showed in
#'   the example
#' @details Sometimes, data sets store several rows where each of them does not
#'   came from the same experiment, that is, each row has nothing to do with the
#'   other rows. Conversely, the rows from RNA-seq data sets are dependent to
#'   the other rows into the same column (sample), so other function has to be
#'   apply here, such as \code{\link{EMObuODlmTest}}.
#' @author Enrique Perez_Riesgo
#' @family One sample functions
#' @seealso \code{\link{Fesc}} for use its outcomes as argument of this
#'   function and test the differential expression for independents rows.
#' @param expr matrix where the first set of columns store the counts for each
#'   sample and the second set of columns store the library size for each
#'   sample. For instance, if "n1" is the library size from sample 1 and x1 is
#'   the counts from sample 1, and j is the number of samples which belong to
#'   the same condition that sample 1, the column 1 is x1 and "j+1" is n1, the
#'   column 2 is x2 and "j+2" is n2, and so on.
#' @param condition Binary vector where 0 means control and 1 treatment
#' @param Fescd is the data.frame which stores coefficients calculated by Fesc
#'   function from dEMO package
#' @param cd it is the common disperssion estimated by
#'   \code{\link[edgeR]{estimateCommonDisp}} function from edgeR package
#' @param alpha significance level for hypothesis test. Default value is 0.05
#' @param method.adj It is the multiple testing correction method. Default value
#'   correspond to Benjamini-Hochberg correction. c("holm", "hochberg",
#'   "hommel", "bonferroni", "BH", "BY","fdr", "none")
#' @return A data.frame is returned, which cointains log2FC, FC_LowBound and
#'   FC_UpBound
#' @export dEMOOD1sci
#' @examples
#' #this example is the continuation of example from Fesc function
#' cd <- estimateCommonDisp(y = expressiondata[, 1:6])
#' for(z in 1:dim(expressiondata)[1]){
#' testdEMOODci <- dEMOOD1sci(expr = expressiondata[z, 1:6],
#' Fescd = Fescoefs, condition = c(0, 0, 0, 1, 1, 1),
#' n = expressiondata[z, 7:12], cd = cd, alpha = 0.05, method.adj = "none")
#' }

dEMOOD1sci <- function(expr, condition, Fescd, n=NULL, cd, alpha=0.05,
                       method.adj="BH"){
  if(sum(expr[condition == 0]) == 0){
    expr[1] <- 1/2
    n[1] <- n[1] - 1/2
  }
  if(sum(expr[condition == 1]) == 0){
    expr[(sum(condition == 0) + 1)]<- 1/2
    n[(sum(condition == 0) + 1)] <- n[(sum(condition == 0) + 1)] - 1/2
  }

  #Overdispersion
  cdg0 <- edgeR::estimateTagwiseDisp.default(y = t(expr[condition == 0]),
          dispersion = cd, lib.size = n[condition == 0])
  cdg1 <- edgeR::estimateTagwiseDisp.default(y = t(expr[condition == 1]),
          dispersion = cd, lib.size = n[condition == 1])
  cdg  <- edgeR::estimateTagwiseDisp.default(y = t(expr), dispersion = cd,
          lib.size = n)
  r1   <- 1/cdg1
  r0   <- 1/cdg0
  rt   <- length(condition)/cdg

  probabilities <- 1 - (expr)/n
  variances <- (n - expr)*(1 - probabilities)/(probabilities^2)

  n1   <- sum(n[condition == 1])
  x1   <- sum(expr[condition == 1])
  A1   <- ((1/x1)*(1/log(2)))^2
  uOD1 <- (sum(expr[condition == 1]))*(1 +
          (sum(expr[condition == 1]))*(1/r1))
  variances1 <- sum(variances[condition == 1])
  urep1 <- var(expr[condition == 1])

  n0   <- sum(n[condition == 0])
  x0   <- sum(expr[condition == 0])
  A0   <- ((1/x0)*(1/log(2)))^2
  uOD0 <- (sum(expr[condition == 0]))*(1+
          (sum(expr[condition == 0]))*(1/r0))
  variances0 <- sum(variances[condition == 0])
  urep0 <- var(expr[condition == 0])

  FC <- log((x1/(n1))/(x0/(n0)), 2)

  Fesc1 <- r1*Fescd[2, 2] + Fescd[1, 2]
  Fesc0 <- r0*Fescd[2, 1] + Fescd[1, 1]

  if(is.na(Fesc1)){
    Fesc1=1
  }

  if(is.na(Fesc0)){
    Fesc0=1
  }

  u1 <- ((uOD1)/sum(condition == 1)) + (urep1 + variances1)/
        (sum(condition == 1)/Fesc1)
  u0 <- ((uOD0)/sum(condition == 0)) + (urep0 + variances0)/
    (sum(condition == 0)/Fesc0)
  uc <- sqrt(A1*u1 + A0*u0)

  k  <- qnorm(alpha/(2), lower.tail = FALSE)
  Uc <- k*uc

  CI <- data.frame(FCi = FC - Uc, FCs = FC + Uc)

  Conf.Int.SAM <- data.frame(log2FC = FC, LowBound = CI[1], UpBound = CI[2])

  resultados <- Conf.Int.SAM
  return(resultados)
}
