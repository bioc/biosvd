##' Excludes specified eigenfeatures/eigenassays from the original data
##'
##' The function excludes eigenfeatures/eigenassays from the data that correspond to steady-state intensity, steady-scale variance, experimental artifacts and/or noise as specified by the user. In case no eigenfeatures are specified, the eigenfeature(s) corresponding to steady-state/steady-scale is/are removed. Filtering out steady-state expression/intensity corresponds to centering the expression/intensity patterns at steady-state expression/intensity level (arithmetic mean of expression/intensity ~ 0). Filtering out steady-scale variance corresponds to normalization by the steady scale of expression/intensity variance (geometric mean of variance ~ 1).
##' @title Excludes specified eigenfeatures/eigenassays from the data
##' @param x object of class eigensystem
##' @param excludeEigenfeatures vector of positive integers representing eigenfeatures to be excluded from the data
##' @return Object of class eigensystem
##' @import Biobase grid BiocGenerics
##' @importFrom methods callGeneric show
##' @export
##' @docType methods
##' @rdname exclude-methods
##' @aliases exclude exclude,Eigensystem-method
##' @section Methods:
##' \describe{
##' \item{\code{signature(x = "Eigensystem")}}{
##' }
##' }
##' @family "eigensystem"
##' @examples
##' ## Metabolomics starvation data obtained from http://genomics-pubs.princeton.edu/StarvationMetabolomics/Download.shtml
##' data(StarvationData)
##'
##' ## Computes the eigensystem for the actual data
##' eigensystem <- compute(StarvationData)
##' ## Excludes the eigenfeature(s) representing steady-state expression/intensity as defined in compute
##' exclude(eigensystem)
##' ## Excludes user-specified eigenfeatures 1, 4 and 5
##' exclude(eigensystem, excludeEigenfeatures=c(1,4,5))
##'
##' ## Computes the eigensystem for the variance in the data
##' eigensystem <- compute(StarvationData, apply="variance")
##' ## Excludes the eigenfeature(s) representing steady-scale variance as defined in compute
##' eigensystem <- exclude(eigensystem)
##' ## Excludes none of the eigenfeatures and recalculates the eigensystem for the actual data
##' eigensystem <- exclude(eigensystem, excludeEigenfeatures=0)
##' @references Alter O, Brown PO and Botstein D. Singular value decomposition for genome-wide expression data processing and modeling. Proc Natl Acad Sci U.S.A. 97(18), 10101-10106 (2000).
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
##'
setMethod("exclude", "Eigensystem", 
          function(x, excludeEigenfeatures=NULL)
{
  eigensystem <- x
  .checkExcludeEigenfeatures(excludeEigenfeatures, eigensystem)
  
  if (is.numeric(excludeEigenfeatures)) excludeEigenfeatures(eigensystem) <- excludeEigenfeatures
  if (excludeEigenfeatures(eigensystem)[1]!=0) eigenexpressions(eigensystem)[excludeEigenfeatures(eigensystem)] <- 0
  matrix <- eigenassays(eigensystem) %*% diag(eigenexpressions(eigensystem)) %*% eigenfeatures(eigensystem)
  if (apply(eigensystem)=="variance") matrix <- signMatrix(eigensystem) * sqrt(exp(matrix))
  
  rownames(matrix) <- rownames(matrix(eigensystem))
  colnames(matrix) <- colnames(matrix(eigensystem))
  matrix(eigensystem) <- matrix
  eigensystem <- compute(eigensystem)
})