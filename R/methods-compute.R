.compute <- function(object, apply=c("data","variance")) {
  apply <- match.arg(apply)
  result <- .transformObject(object)
  assay.matrix <- result[[1]]
  feature.matrix <- result[[2]]
  matrix <- result[[3]]
  sign.matrix <- sign(matrix)
  matrix <- .transformMatrix(matrix, apply)
  
  svd.results <- svd(matrix)
  eigenassays <- svd.results$u
  eigenexpressions <- svd.results$d
  eigenfeatures <- t(svd.results$v)
  colnames(eigenfeatures) <- colnames(matrix)
  rownames(eigenfeatures) <- c(1:min(nrow(matrix),ncol(matrix)))
  fractions <- eigenexpressions^2/sum(eigenexpressions^2)
  names(fractions) <- c(1:min(nrow(matrix),ncol(matrix)))
  
  new("Eigensystem",
      matrix = matrix,
      signMatrix = sign.matrix,
      assayMatrix = assay.matrix,
      featureMatrix = feature.matrix,
      eigenassays = eigenassays,
      eigenexpressions = eigenexpressions,
      eigenfeatures = eigenfeatures,
      assaycorrelations = diag(eigenexpressions) %*% eigenfeatures,
      featurecorrelations = t(eigenassays %*% diag(eigenexpressions)),
      fractions = fractions,
      entropy = round(100*(-sum(fractions*log(fractions))/log(ncol(eigenassays))))/100,
      apply = apply,
      excludeEigenfeatures = if(fractions[2] > 0.15 && sum(fractions[1:2]) > 0.85) {c(1,2)} else {1}
      )
}


##' Computes the eigensystem for a feature by assay matrix by applying Singular Value Decomposition.
##'
##' Function compute decomposes the input data set from the feature x assay space to the reduced diagonalized "eigenfeatures x eigenassays" space, with the eigenfeatures and eigenassays unique orthonormal superpositions of the features and assays, respectively. This approach allows filtering out eigenfeatures and eigenassays that are inferred to represent noise or experimental artifacts, either at the expression/intensity level or the variance level.
##' The function can be applied to an object of class matrix, data.frame, ExpressionSet, or eigensystem.
##' @title Compute the eigensystem for a feature by assay matrix
##' @param object object of class matrix, data.frame, Expressionset or Eigensystem, containing the feature x assay expression or intensity data
##' @param apply the actual expression/intensity data (data) or variance in the data (variance) to which to apply the function (default = data)
##' @return Object of class Eigensystem
##' @import Biobase grid BiocGenerics
##' @importFrom methods callGeneric show
##' @export
##' @docType methods
##' @rdname compute-methods
##' @aliases compute compute,Eigensystem-method compute,data.frame-method compute,matrix-method compute,ExpressionSet-method
##' @section Methods:
##' \describe{
##' \item{\code{signature(x = "Eigensystem")}}{
##' }
##' \item{\code{signature(x = "data.frame")}}{
##' }
##' \item{\code{signature(x = "matrix")}}{
##' }
##' \item{\code{signature(x = "ExpressionSet")}}{
##' }
##' }
##' @family "eigensystem"
##' @seealso Eigensystem-class
##' @examples
##' ## Metabolomics starvation data obtained from http://genomics-pubs.princeton.edu/StarvationMetabolomics/Download.shtml
##' data(StarvationData)
##'
##' ## Computes the eigensystem for the actual expression/intensity data
##' eigensystem <- compute(StarvationData)
##' ## Computes the eigensystem for the variance in the data
##' eigensystem <- compute(StarvationData, apply="variance")
##' @references Alter O, Brown PO and Botstein D. Singular value decomposition for genome-wide expression data processing and modeling. Proc Natl Acad Sci U.S.A. 97(18), 10101-10106 (2000).
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
setMethod("compute", "Eigensystem", .compute)
setMethod("compute", "data.frame", .compute)
setMethod("compute", "ExpressionSet", .compute)
setMethod("compute", "matrix", .compute)

