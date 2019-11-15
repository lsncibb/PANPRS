# PANPRS




### Input for PANPRS incorporating multiple traits and functional annotations of SNPs.

summaryZ, 
The Z statistics of p SNPs from q GWA studies. A matrix with dimension p x q for p SNPs and q traits. The first column corresponds to the primary trait and the rest columns correspond to the secondary traits.


Nvec,
A vector of length q for the sample sizes of q GWA studies.

plinkLD,
LD matrix information.

NumIter,
The number of maximum iterations for the estimation procedure. 


funcIndex,
Inputs for the functional annotations of SNPs. A p x k matrix with (0,1) entry; p is the number of SNPs and 
k is the number of functional annotations. For the element at i-th row, j-th column, the entry 0 means SNP i without j-th functional annotation; entry 1 means otherwise.

numfunc,
The number of functional annotations.


dfMax
The upper bound of the number of non-zero estimates of coefficients for the primary trait.



### Usage:
The current version only work on Unix, Linux and Mac System, R(>=3.4.3), R packages "gtools" and "permutations" and GCC(>=4.4.7) are required.

Modify the parameters in the gsfPEN.R, and execute it.



## Example
Please find it in the R package.


## Future extensions



## Contact
* Lei Song, lei.song@nih.gov
* Ting-huei Chen, ting-huei.chen@mat.ulaval.ca

