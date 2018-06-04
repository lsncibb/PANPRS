\name{gsfPEN}
\alias{gsfPEN}
\title{
SummaryLasso incorporating multiple traits and functional annotations of SNPs.
}
\description{
	SummaryLasso to model pleiotropy by introducing a group-Lasso type penalty, which is sensitive to select SNPs modestly associated with multiple traits and to incorporate functional annotations of SNPs simultaneously.
}
\usage{
gsfPEN(summaryZ, Nvec, plinkLD, NumIter = 1000, RupperVal = NULL, breaking = 1, numChrs = 22, ChrIndexBeta = 0, Init_summaryBetas = 0, Zscale = 1, tuningMatrix = NULL, penalty = "mixLOG", funcIndex, numfunc, p.Threshold = NULL, p.Thresholdpara = c(0.5, 10^4, 5), taufactor = c(1/25, 1, 10), llim_length = 10, subtuning = 50, Lambda_limit = c(1, 4), Lenlam = 100, lambdavec_func = NULL, lambdavec_func_limit_len = c(1.5, 4), dfMax = NULL, outputAll = 0)
}
\arguments{
  \item{summaryZ}{
The Z statistics of p SNPs from q GWA studies. A matrix with dimension p x q for p SNPs and q traits. The first column corresponds to the primary trait and the rest columns correspond to the secondary traits.
}
  \item{Nvec}{
A vector of length q for the sample sizes of q GWA studies.
}
  \item{plinkLD}{
.ld file obtained from the LD calculation from plink.
}
  \item{NumIter}{
The number of maximum iterations for the estimation procedure. 
}
  \item{RupperVal}{
The maximum tolerable magnitude of the estimates of coefficients during the iterations. 
This is to avoid certain estimates of coefficients to diverge during the iterations. This may happen when the signs of the correlation coefficients were estimated incorrectly. The default value is 50 times the maximum of coefficients from the input in absolute values.
}
  \item{breaking}{
A binary (0,1) variable to check if there are some certain estimates of coefficients to diverge during the iterations. This may happen when the signs of the correlation coefficients were estimated incorrectly. The default value is 1.
}
  \item{numChrs}{
The number of chromosomes used in the analysis. Current version of package does not use this argument. 
}
  \item{ChrIndexBeta}{
The chromosome index for each SNP. Current version of packge does not use this argument. 
}
  \item{Init_summaryBetas}{
Can be used to set the initial values of the coefficients for the iterative estimation. 
}
  \item{Zscale}{
A binary (0,1) variable to make the coefficients from different GWA studies with unequal sample sizes comparable.
 The default value is 1.
}
  \item{tuningMatrix}{
Inputs for the tuning values of the tuning parameters. Default is null and it will be generated automatically.
}
  \item{penalty}{
Current version of pacakge does not use this argument.
}
  \item{funcIndex}{
Inputs for the functional annotations of SNPs. A p x k matrix with (0,1) entry; p is the number of SNPs and 
k is the number of functional annotations. For the element at i-th row, j-th column, the entry 0 means SNP i without j-th functional annotation; entry 1 means otherwise.
}
  \item{numfunc}{
The number of functional annotations.
}
  \item{p.Threshold}{
The p-values threshold to set up the tuning values of the baseline tuning parameter.
}
  \item{p.Thresholdpara}{
When p.Threshold is null, p.Threshold will be generated automatically based on the values of p.Thresholdpara.
The default values are c(0.5, 10^4, 5), where the first element is the maximum of the p-value threshold, 
the second element is the minimum, and the third element is total number of p-value thresholds to be 
generated from the minimum to the maximum.
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
  \item{lambdavec_func}{
The tuning values for the tuning parameters associated with the functional annotations.
}
  \item{lambdavec_func_limit_len}{
When lambdavec_func is null, lambdavec_func will be generated automatically based on the arguments of 
\code{lambdavec_func_limit_len}. The default values are c(1.5, 4). The first element is the maximum of the 
tuning value and the second element is the total number of the tuning values to be 
generated from 0 to the maximum.
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
data("funcIndex")
output = gsfPEN(summaryZ=summaryZ, Nvec=Nvec, plinkLD=plinkLD,funcIndex=funcIndex, numfunc=ncol(funcIndex))

## The function is currently defined as
function (summaryZ, Nvec, plinkLD, NumIter = 1000, RupperVal = NULL, 
    breaking = 1, numChrs = 22, ChrIndexBeta = 0, Init_summaryBetas = 0, 
    Zscale = 1, tuningMatrix = NULL, penalty = "mixLOG", funcIndex, 
    numfunc, lambdavec = c(0, 0.05, 0.075), p.Threshold = NULL, 
    p.Thresholdpara = c(0.5, 10^4, 5), taufactor = c(1/25, 1, 
        10), llim_length = 10, subtuning = 50, Lambda_limit = c(0.5, 0.9),
         Lenlam = 100, lambdavec_func = NULL, lambdavec_func_limit_len = c(1.5, 
        4), dfMax = NULL, outputAll = 0) 
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
    if (is.null(p.Threshold)) {
        p.Threshold = seq(p.Thresholdpara[1], p.Thresholdpara[2], 
            length.out = p.Thresholdpara[3])
    }
    if (any(c(is.null(tuningMatrix), is.null(lambdavec_func)))) {
        medianval = median(apply(abs(summaryBetas), 1, sum), 
            na.rm = T)
        tauvec = sort(medianval * taufactor)
        Lambda_limit = quantile(abs(summaryZ[, 1]), Lambda_limit)
        output = Tuning_setup_group_func(lambdavec_func, lambdavec_func_limit_len, 
            p.Threshold, numfunc, tauvec, subtuning, Lambda_limit, 
            Lenlam)
        funcLambda = output[[1]]
        lambdavec_func = output[[2]]
        tuningMatrix = output[[3]]
    }
    else {
        funcLambda0 = permutations(length(lambdavec_func), numfunc, 
            repeats.allowed = T)
        funcLambda = funcLambda0 - 1
    }
    Betaindex = c(1:nrow(summaryBetas)) - 1
    SNPnames = rownames(summaryBetas)
    ldJ = PlinkLD_transform(plinkLD, SNPnames)
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
    if (nrow(funcIndex) != nrow(summaryBetas)) {
        stop("nrow of summaryBetas and row of funcIndex do not match.")
    }
    Zmatrix = 1 - funcIndex
    sumfuncIndex = apply(Zmatrix, 1, sum)
    IfuncSNP = rep(0, P)
    IfuncSNP[which(sumfuncIndex != 0)] = 1
    nrow_funcLambda = nrow(funcLambda)
    ncol_funcLambda = ncol(funcLambda)
    leng.p.Threshold = length(p.Threshold)
    nrow_Zmatrix = nrow(Zmatrix)
    ncol_Zmatrix = ncol(Zmatrix)
    NumTuning = nrow(tuningMatrix) * leng.p.Threshold * nrow_funcLambda
    nrow_AllTuningMatrix = NumTuning
    ncol_AllTuningMatrix = (numfunc + 1) + 2
    AllTuningMatrix = matrix(0, nrow_AllTuningMatrix, ncol_AllTuningMatrix)
    dims2 = rep(0, 8)
    dims2[1] = NumTuning
    dims2[2] = nrow_funcLambda
    dims2[3] = ncol_funcLambda
    dims2[4] = nrow_Zmatrix
    dims2[5] = ncol_Zmatrix
    dims2[6] = leng.p.Threshold
    dims2[7] = nrow_AllTuningMatrix
    dims2[8] = ncol_AllTuningMatrix
    lambda0vec = abs(-qnorm(p.Threshold/2))
    nrow_BetaMatrix = NumTuning
    ncolBetaMatrix = P * Q
    dims = rep(0, 13)
    dims[1] = P
    dims[2] = NumIter
    dims[3] = breaking
    dims[4] = nrow_IndexMatrix
    dims[5] = ncol_IndexMatrix
    dims[6] = length(wind)
    dims[7] = numChrs
    dims[8] = Zscale
    dims[9] = nrow(tuningMatrix)
    dims[10] = ncol(tuningMatrix)
    dims[11] = Q
    dims[12] = ncolBetaMatrix
    dims[13] = dfMax
    Numitervec = rep(0, NumTuning)
    BetaMatrix = matrix(0, nrow_BetaMatrix, ncolBetaMatrix)
    Z = .C("gsfFunc", as.double(t(summaryBetas)), as.integer(ldJ), 
        as.integer(dims), Numitervec = as.integer(Numitervec), 
        as.integer(t(IndexMatrix)), as.integer(IndJ), as.double(ldvec), 
        as.integer(ChrIndexBeta), as.double(RupperVal), as.double(Init_summaryBetas), 
        as.double(t(SDvec)), as.double(t(tuningMatrix)), BetaMatrix = as.double(t(BetaMatrix)), 
        penalty, as.double(lambda0vec), as.double(t(Zmatrix)), 
        as.integer(dims2), AllTuningMatrix = as.double(t(AllTuningMatrix)), 
        as.double(lambdavec_func), as.integer(t(funcLambda)), 
        as.integer(IfuncSNP), PACKAGE = "SummaryLasso")
    BetaMatrix = matrix(Z$BetaMatrix, nrow = NumTuning, ncol = ncolBetaMatrix, 
        byrow = TRUE)
    colnames(BetaMatrix) = paste0(rep(SNPnames, times = Q), ".trait", 
        rep(c(1:Q), each = P))
    AllTuningMatrix = matrix(Z$AllTuningMatrix, nrow = nrow_AllTuningMatrix, 
        ncol = ncol_AllTuningMatrix, byrow = TRUE)
    Numitervec = Z$Numitervec
    colnames(AllTuningMatrix) = c("lam0", paste0("lam.f", c(1:numfunc)), 
        "lam2", "tau")
    if (outputAll == 0) {
        convergeIndex = which(Numitervec > 0)
        Numitervec = Numitervec[convergeIndex]
        BetaMatrix = BetaMatrix[convergeIndex, ]
        AllTuningMatrix = AllTuningMatrix[convergeIndex, ]
    }
    ll = list(BetaMatrix = BetaMatrix, Numitervec = Numitervec, 
        AllTuningMatrix = AllTuningMatrix)
    return(ll)
  }
}

