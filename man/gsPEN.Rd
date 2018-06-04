\name{gsPEN}
\alias{gsPEN}
\title{
SummaryLasso incorporating multiple traits
}
\description{
	SummaryLasso to model pleiotropy by introducing a group-Lasso type penalty, which is sensitive to select SNPs modestly associated with multiple traits. 
}
\usage{
gsPEN(summaryZ, Nvec, plinkLD, NumIter = 100, breaking = 1, numChrs = 22, ChrIndexBeta = 0, Init_summaryBetas = 0, Zscale = 1, RupperVal = NULL, tuningMatrix = NULL, penalty = c("mixLOG"), OutputForPath = 0, taufactor = c(1/25, 1, 10), llim_length = 10, subtuning = 50, Lambda_limit = c(0.5, 0.9), Lenlam = 5, dfMax = NULL, outputAll = 0)
}
\arguments{
  \item{summaryZ}{
The Z statistics of p SNPs from q GWA studies. A matrix with dimension p x q for p SNPs and q traits. The first column corresponds to the primary trait and the rest columns correspond to the secondary traits.
}
  \item{Nvec}{
A vector of length q for the sample sizes of q GWA studies.
}
  \item{plinkLD}{
.ld file of the LD calculation from plink.
}
  \item{NumIter}{
The number of maximum iteraions for the estimation procedure. 
}
  \item{breaking}{
A binary (0,1) variable to check if there are some certain estimates of coefficients to diverge during the iterations. This may happen when the signs of the correlation coefficinets were estimated incorrectly. The default value is 1.
}
  \item{numChrs}{
The number of chromosomes used in the analysis. Current version of pacakge does not use this argument.
}
  \item{ChrIndexBeta}{
The chromosome index for each SNP. Current version of pacakge does not use this argument.
}
  \item{Init_summaryBetas}{
Can be used to set the initial values of the coefficients for the iterative estimation. 
}
  \item{Zscale}{
A binary (0,1) variable to make the coefficients from different GWA studies with unequal sample sizes comparable.
 The default value is 1.
}
  \item{RupperVal}{
The maximum tolerable magnitude of the estimates of coefficients during the iterations. 
This is to avoid a certain estimates of coefficients to diverge during the iterations. This may happen when the signs of the correlation coefficinets were estimated incorrectly. The default value is 50 times the maximum of coeffcients from the input in absolute values.
}
  \item{tuningMatrix}{
Inputs for the tuning values of the tuning parameters. Default is null and it will be generated automatically. 
}
  \item{penalty}{
Current version of pacakge does not use this argument. 
}
  \item{OutputForPath}{
This is for internal checking. The default value is 0.
}
  \item{taufactor}{
The weights to generate the tuning values for the tuning paramter "tau" and the default is c(1/25, 1, 10) times the median of the p summation of the coefficients for each SNP across q traits.
}
  \item{llim_length}{
The argument to set up the number of tuning values for lambdas between the lower and upper bound. The default value is 10.
}
  \item{subtuning}{
The argument to set up the number of tuning values for lambdas between the lower and upper bound. The default value is 50.
}
  \item{Lambda_limit}{
The quantiles to set up the tuning values of lambda. The default value is c(0.5, 0.9).
}
  \item{Lenlam}{
The number of tuning values for lambda parameter without using the Log penalty. In other words, the initial \code{Lenlam} rows of the tuningMatrix are for summayLasso single trait analysis.
}
  \item{dfMax}{
The upper bound of the number of non-zero estimates of coefficients for the primary trait.
}
  \item{outputAll}{
For internal usage. 
}
}
\details{
Note that the tuning values for the tuning parameters may need to be modified manually when the selected optimal tuning parameters are at the boundary of the inputs.
}
\value{
\item{BetaMatrix}{The output of the coefficients matrix with dimensions (total number of combinations 
	of the tuning values times (pq)). Each column represents the vectorization of the p x q coefficients matrix given a particular combination of the tuning values (stacking its columns into a column vector). 
	} 
\item{Numitervec}{This vector shows the number of iterations to converge for each combination 
	of the tuning values.}
\item{AllTuningMatrix}{This matrix shows all combination of tuning values used in the estimation process. Its dimension is that total number of combinations 
	of the tuning values times total number of tuning parameters.}
}
\references{
This R packages is based on the method introduced in the manuscript "A comprehensive statistical framework for building polygenic risk prediction models based on summary statistics of genome-wide association studies."
}
\author{
Ting-Huei Chen
}



\examples{
library("SummaryLasso") 
data("summaryZ")
data("Nvec")
data("plinkLD")
output = gsPEN(summaryZ=summaryZ, Nvec=Nvec, plinkLD=plinkLD)


## The function is currently defined as
function (summaryZ, Nvec, plinkLD, NumIter = 100, breaking = 1, 
    numChrs = 22, ChrIndexBeta = 0, Init_summaryBetas = 0, Zscale = 1, 
    RupperVal = NULL, tuningMatrix = NULL, penalty = c("mixLOG"), 
    OutputForPath = 0, taufactor = c(1/25, 1, 10), llim_length = 10, 
    subtuning = 50, Lambda_limit = c(0.5, 0.9), Lenlam = 5, dfMax = NULL, 
    outputAll = 0) 
{
    if (Zscale != 1) {
        error("Tuning values set-up for multiple traits analysis requires Zscale=1.")
    }
    Nq = length(Nvec)
    if (Nq == 1) {
        summaryBetas = matrix(summaryZ/sqrt(Nvec), length(summaryZ), 
            Nq)
        SDvec = matrix(1/sqrt(Nvec), length(summaryZ), Nq)
    }
    else {
        summaryBetas = matrix(0, nrow(summaryZ), Nq)
        SDvec = matrix(0, nrow(summaryZ), Nq)
        for (ii in 1:Nq) {
            summaryBetas[, ii] = summaryZ[, ii]/sqrt(Nvec[ii])
            SDvec[, ii] = 1/sqrt(Nvec[ii])
        }
    }
    rownames(summaryBetas) = rownames(summaryZ)
    if (is.null(dfMax)) {
        dfMax = ceiling(0.7 * nrow(summaryZ))
    }
    if (is.null(tuningMatrix)) {
        medianval = median(apply(abs(summaryBetas), 1, sum), 
            na.rm = T)
        tauvec = sort(medianval * taufactor)
        Lambda_limit = quantile(abs(summaryZ[, 1]), Lambda_limit)
        tuningMatrix = Tuning_setup_group_only(tauvec, subtuning, 
            Lambda_limit, Lenlam, llim_length)
    }
    inv_summaryBetas = 0
    count_nonzero = function(xx) {
        length(which(xx != 0))
    }
    counts = apply(summaryBetas, 1, count_nonzero)
    Betaindex = c(1:nrow(summaryBetas)) - 1
    SNPnames = rownames(summaryBetas)
    ldJ = PlinkLD_transform(plinkLD, SNPnames)
    rm(plinkLD)
    JidMatrix = matrix(, nrow(ldJ), 2)
    mat1 = match(ldJ[, 1], SNPnames)
    mat2 = match(ldJ[, 2], SNPnames)
    JidMatrix[, 1] = Betaindex[mat1]
    JidMatrix[, 2] = Betaindex[mat2]
    ldJ[, 1] = JidMatrix[, 1]
    ldJ[, 2] = JidMatrix[, 2]
    od = order(JidMatrix[, 1], JidMatrix[, 2], decreasing = F)
    ldJ = ldJ[od, ]
    wind = which(!Betaindex \%in\% ldJ[, 1])
    IndJ = -1
    if (length(wind) > 0) {
        IndJ = Betaindex[wind]
    }
    Counts = table(ldJ[, 1])
    NumSNP = length(Counts)
    IndexS = c(0, cumsum(Counts)[-NumSNP])
    IndexE = cumsum(Counts) - 1
    IndexMatrix = matrix(, NumSNP, 3)
    IndexMatrix[, 1] = as.numeric(names(Counts))
    nrow_IndexMatrix = NumSNP
    ncol_IndexMatrix = ncol(IndexMatrix)
    IndexMatrix[, 2] = IndexS
    IndexMatrix[, 3] = IndexE
    ldvec = ldJ[, 3]
    ldJ = ldJ[, 2]
    length_ldJ = length(ldJ)
    if (is.null(RupperVal)) {
        RupperVal = ceiling(max(abs(summaryBetas), na.rm = T) * 
            50)
    }
    P = nrow(summaryBetas)
    Q = ncol(summaryBetas)
    if (nrow(tuningMatrix) > 1) {
        NumTuning = nrow(tuningMatrix)
    }
    ncolBetaMatrix = P * Q
    dims = rep(0, 14)
    dims[1] = NumTuning
    dims[2] = P
    dims[3] = NumIter
    dims[4] = breaking
    dims[5] = nrow_IndexMatrix
    dims[6] = ncol_IndexMatrix
    dims[7] = length(wind)
    dims[8] = Zscale
    dims[9] = nrow(tuningMatrix)
    dims[10] = ncol(tuningMatrix)
    dims[11] = Q
    dims[12] = ncolBetaMatrix
    dims[13] = OutputForPath
    dims[14] = dfMax
    Numitervec = rep(0, NumTuning)
    BetaMatrix = matrix(0, NumTuning, ncolBetaMatrix)
    Z = .C("gsPEN", as.double(t(summaryBetas)), as.integer(ldJ), 
        as.integer(dims), Numitervec = as.integer(Numitervec), 
        as.integer(t(IndexMatrix)), as.integer(IndJ), as.double(ldvec), 
        as.double(inv_summaryBetas), as.integer(ChrIndexBeta), 
        as.double(RupperVal), as.double(Init_summaryBetas), as.double(t(SDvec)), 
        as.double(t(tuningMatrix)), BetaMatrix = as.double(t(BetaMatrix)), 
        penalty, PACKAGE = "SummaryLasso")
    BetaMatrix = matrix(Z$BetaMatrix, nrow = NumTuning, ncol = ncolBetaMatrix, 
        byrow = TRUE)
    colnames(tuningMatrix) = c("lam1", "null", "lam2", "tau")
    tuningMatrix = tuningMatrix[, c(1, 3, 4)]
    colnames(BetaMatrix) = paste0(rep(SNPnames, times = Q), ".trait", 
        rep(c(1:Q), each = P))
    Numitervec = Z$Numitervec
    if (outputAll == 0) {
        convergeIndex = which(Numitervec > 0)
        Numitervec = Numitervec[convergeIndex]
        BetaMatrix = BetaMatrix[convergeIndex, ]
        tuningMatrix = tuningMatrix[convergeIndex, ]
    }
    ll = list(BetaMatrix = BetaMatrix, Numitervec = Numitervec, 
        tuningMatrix = tuningMatrix)
    return(ll)
  }
}

