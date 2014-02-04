##' Sorts the data according to two eigenfeatures and eigenassays
##'
##' The function sorts the data according to two certain eigenfeatures and eigenassays. This gives a global picture of the dynamics of expression/intensities, in which individual features and assays are classified in groups of similar regulation and function or similar cellular state and biological phenotype.
##' @title Sorts the data according to two eigenfeatures and eigenassays
##' @param x object of class eigensystem
##' @param decreasing FALSE
##' @param eigenfeature.xaxis first eigenfeature to sort against (default 2)
##' @param eigenfeature.yaxis second eigenfeature to sort against (default 1)
##' @param colorIdFeatures vector of integers or a string variable in the feature data of the original data set in case from class ExpressionSet, with feature annotation used for coloring purposes on the polar plot
##' @return Object of class eigensystem
##' @import Biobase grid BiocGenerics
##' @importFrom methods callGeneric show
##' @export
##' @docType methods
##' @rdname sort-methods
##' @aliases sort sort,Eigensystem-method
##' @section Methods:
##' \describe{
##' \item{\code{signature(x = "Eigensystem")}}{
##' }
##' }
##' @examples
##' ## Data obtained from http://genomics-pubs.princeton.edu/StarvationMetabolomics/Download.shtm
##' data(StarvationData)
##'
##' ## Computes the eigensystem for the actual data
##' eigensystem <- compute(StarvationData)
##' ## Excludes the eigenfeatures representing steady-state intensity
##' eigensystem <- exclude(eigensystem)
##'
##' ## Sort the data according to eigenfeature 1 and 2
##' eigensystem.sorted <- sort(eigensystem)
##' ## Sort the data according to eigenfeature 3 and 4
##' eigensystem.sorted <- sort(eigensystem, eigenfeature.xaxis=4, eigenfeature.yaxis=3)
##' ## Visualization of the data after sorting to eigenfeature 3 and 4
##' plot(eigensystem.sorted)
##' @references Alter O, Brown PO and Botstein D. Singular value decomposition for genome-wide expression data processing and modeling. Proc Natl Acad Sci U.S.A. 97(18), 10101-10106 (2000).
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
##'
setMethod("sort", "Eigensystem",
          function(x, decreasing=FALSE, eigenfeature.xaxis=2,
                   eigenfeature.yaxis=1,
                   colorIdFeatures=rep(1,nrow(matrix(eigensystem))))
{
  eigensystem <- x
  if (isTRUE(decreasing))
    stop("'decreasing' argument not yet supported")
  .checkEigenfeatureAxes(eigenfeature.xaxis, eigenfeature.yaxis, eigensystem)
  colorIdFeatures <- .checkColorID(colorIdFeatures, eigensystem, "features")

  zerophase <- atan(assaycorrelations(eigensystem)[eigenfeature.yaxis,1]/assaycorrelations(eigensystem)[eigenfeature.xaxis,1])/pi
  nr.matrix <- nrow(matrix(eigensystem))
  coordinates.features <- base::matrix(0,nrow=nr.matrix,ncol=2)
  amplitude.scale.factor <- sqrt(featurecorrelations(eigensystem)[eigenfeature.xaxis,]^2+featurecorrelations(eigensystem)[eigenfeature.yaxis,]^2)
  for (z in c(1:nr.matrix)) coordinates.features[z,] <- c(featurecorrelations(eigensystem)[eigenfeature.xaxis,z]/amplitude.scale.factor[z], featurecorrelations(eigensystem)[eigenfeature.yaxis,z]/amplitude.scale.factor[z])
  angular.distance <- atan(featurecorrelations(eigensystem)[eigenfeature.yaxis,]/featurecorrelations(eigensystem)[eigenfeature.xaxis,])/pi
  squared.distance <- sqrt(coordinates.features[,1]^2+coordinates.features[,2]^2)
  sort.matrix <- cbind(coordinates.features, angular.distance, squared.distance, matrix(eigensystem))
  rownames(sort.matrix) <- rownames(matrix(eigensystem))
  sort.matrix <- sort.matrix[order(sort.matrix[,1],sort.matrix[,2]),]
  negative1 <- max(which(sort.matrix[,1]<0))
  positive1 <- min(which(sort.matrix[,1]>0))
  sort.matrix <- sort.matrix[,-1]

  if (is.infinite(negative1)) {
    sort.matrix <- sort.matrix[order(sort.matrix[,1],sort.matrix[,2]),]
    sort.matrix <- sort.matrix[,-1]
    sort.matrix[,1] <- sort.matrix[,1]+1-zerophase
  } else if (is.infinite(positive1)) {
    sort.matrix <- sort.matrix[order(sort.matrix[,2],sort.matrix[,1]),]
    sort.matrix <- sort.matrix[,-1]
    sort.matrix[,1] <- sort.matrix[,1]-zerophase
  } else {
    sort.matrix <- rbind(sort.matrix[order(sort.matrix[1:negative1,2],sort.matrix[1:negative1,1]),], sort.matrix[positive1-1+order(sort.matrix[positive1:nr.matrix,1],sort.matrix[positive1:nr.matrix,2]),])
    sort.matrix <- sort.matrix[,-1]
    sort.matrix[1:negative1,1] <- sort.matrix[1:negative1,1]-zerophase
    sort.matrix[positive1:nr.matrix,1] <- sort.matrix[positive1:nr.matrix,1]+1-zerophase
  }

  negative2 <- max(which(sort.matrix[,1]<0))
  positive2 <- min(which(sort.matrix[,1]>0))
  rownames.before <- rownames(sort.matrix)
  if (is.infinite(negative2) || is.infinite(positive2)) {
    rownames.after <- rownames.before
  } else {
    rownames.after <- c(rownames.before[positive2:nr.matrix], rownames.before[1:negative2])
    sort.matrix <- rbind(sort.matrix[positive2:nr.matrix,], sort.matrix[1:negative2,])
  } 
  rownames(sort.matrix) <- rownames.after

  positive3 <- max(which(sort.matrix[,1]>0))
  negative3 <- min(which(sort.matrix[,1]<0))
  if (is.finite(negative3) && is.finite(positive3)) {
    sort.matrix[negative3:nr.matrix,1] <- sort.matrix[negative3:nr.matrix,1] + ceiling((sort.matrix[positive3,1]-sort.matrix[negative3,1])*10)/10
  } 
  matrix <- sort.matrix[,3:(ncol(matrix(eigensystem))+2)]
  eigenassays <- eigenfeatures(eigensystem) %*% t(matrix)
  eigenassays <- t(eigenassays/base::matrix(eigenexpressions(eigensystem),nrow=nrow(eigenfeatures(eigensystem)),ncol=nr.matrix))
  assaycorrelations <- diag(eigenexpressions(eigensystem)) %*% eigenfeatures(eigensystem)
  featurecorrelations <- t(eigenassays %*% diag(eigenexpressions(eigensystem)))

  Eigensystem(matrix = matrix,
      signMatrix = signMatrix(eigensystem),
      assayMatrix = assayMatrix(eigensystem),
      featureMatrix = featureMatrix(eigensystem),
      eigenassays = eigenassays,
      eigenexpressions = eigenexpressions(eigensystem),
      eigenfeatures = eigenfeatures(eigensystem),
      assaycorrelations = assaycorrelations,
      featurecorrelations = featurecorrelations,
      fractions = fractions(eigensystem),
      entropy = entropy(eigensystem),
      apply = apply(eigensystem),
      excludeEigenfeatures = excludeEigenfeatures(eigensystem),
      colorIdFeatures = colorIdFeatures[match(rownames(matrix),rownames(matrix(eigensystem)))]
      )
})
