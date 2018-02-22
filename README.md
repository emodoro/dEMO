# dEMO
dEMO is a R package that implements both dEMO test, a new method which allows carrying out a Differential Expression Analysis 
from RNA-seq data sets, and estimating confidencial intervals for Fold Change.

# Functions
Several functions have find developed in dEMO package for ahead three principal points:

1. **Building confidence intervals for Fold Change** :
    * When data sets store data expression such that each row from expression matrix is independent on the rest, we suggest 
    employ \code{\link[dEMO]{dEMOOD1sci}} function. This situation could pops up when each row from data set cames from a independent 
    study such as data simulation following a negative binomial distribution, each row with a n trials and r succeses. 
    * When data sets came from RNA-seq assays, each row *i* depends on the rest of rows, since each row is a feature (gene) 
    measured from the sample stored into column *j*. Here, we suggest to use \code{\link[dEMO]{dEMOODci}} function.
2. **Differential Expression Analysis** : \code{\link[dEMO]{dEMObuODlmTest}} function allows carrying out *dEMO test method*, which is
    a hypothesis testing where hull hpothesis correspond to $\log_2{Fold Change}=0$ and alternative hypothesis to 
    $\log_2{Fold Change}\ne0$.
3. **Comparing different Differential Expression Analysis Methods**. In orther to testing the performance for several Differential
    Expression Analysis methods, between them, several functions have been developed, which allows:
    * Building ROCs: \code{\link[dEMO]{ROCcurve}}.
    * Getting Area under ROCs: \code{\linck[dEMO]{ROCAreas}}, \code{\link[dEMO]{ROCAreas005}} and \code{\link[dEMO]{ROCAreasData}}.
    * Building FDR curves: \code{\link[dEMO]{FDRcurve}}
    * Testing the capability of controlling Type I error: \code{\link[dEMO]{TIerror}}.

# Install dEMO R package.
Since several vignettes have been developed, we suggest download dEMO repository in local computer and then employ \code{\link[devtools]{build}}.
