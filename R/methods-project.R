##' Returns the rectangular and polar coordinates of the eigensystem projection onto one or two eigenfeatures and eigenassays
##'
##' The function allows the sorting of the data according to one or two specified eigenfeatures and eigenassays. This gives a global picture of the dynamics of expression/intensities, in which individual features and assays are classified in groups of similar regulation and function or similar cellular state and biological phenotype.
##' @title Projects the data onto one or two eigenfeatures and eigenassays
##' @param x object of class Eigensystem
##' @param axes numerical vector specifying eigenfeatures to project onto (default c(2,1))
##' @param type string specifying the dimensions to return the coordinates for, either features or assays (default features)
##' @return data.frame
##' @import Biobase grid BiocGenerics
##' @export
##' @docType methods
##' @rdname project-methods
##' @aliases project project,Eigensystem-method
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
##' ## Find the projection of the data onto eigenfeature 1 and 2
##' projection <- project(eigensystem)
##' ## Project the data onto eigenfeature 3 and 4
##' projection <- project(eigensystem, axes=c(4,3))
##' @references Alter O, Brown PO and Botstein D. Singular value decomposition for genome-wide expression data processing and modeling. Proc Natl Acad Sci U.S.A. 97(18), 10101-10106 (2000).
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
##'
setMethod("project", "Eigensystem",
          function(x, axes=c(2,1), type="features") {
            eigensystem <- x
            axes <- .checkEigenfeatureAxes(eigensystem, axes)
            type <- match.arg(type, c("features","assays"), several.ok=FALSE)

            if (type %in% "features") {
              zerophase <- atan(assaycorrelations(eigensystem)[axes[2],1]/assaycorrelations(eigensystem)[axes[1],1])/pi
              cartesian.coords <- t(featurecorrelations(eigensystem)[axes,])/sqrt(diag(matrix(eigensystem) %*% t(matrix(eigensystem))))
              scaled.coords <- t(featurecorrelations(eigensystem)[axes,])/
                sqrt(featurecorrelations(eigensystem)[axes[1],]^2 + featurecorrelations(eigensystem)[axes[2],]^2)
              polar.coords <- cbind(atan(featurecorrelations(eigensystem)[axes[2],]/featurecorrelations(eigensystem)[axes[1],])/pi,
                                    sqrt(cartesian.coords[,1]^2 + cartesian.coords[,2]^2))
              colnames(polar.coords) <- c("theta", "r")
              
              polar.coords[,"theta"] <- polar.coords[,"theta"] - zerophase
              polar.coords[scaled.coords[,1] < 0, "theta"] <- polar.coords[scaled.coords[,1] < 0, "theta"] + 1
              
              coordinates <- cbind(cartesian.coords, polar.coords)
              colnames(coordinates) <- c(paste0("eigenfeature", axes), "theta", "r")
              rownames(coordinates) <- rownames(matrix(eigensystem))
              
              as.data.frame(coordinates)
            } else if (type %in% "assays") {
              zerophase <- atan(featurecorrelations(eigensystem)[axes[2],1]/featurecorrelations(eigensystem)[axes[1],1])/pi
              cartesian.coords <- t(assaycorrelations(eigensystem)[axes,])/sqrt(diag(t(matrix(eigensystem)) %*% matrix(eigensystem)))
              scaled.coords <- t(assaycorrelations(eigensystem)[axes,])/
                sqrt(assaycorrelations(eigensystem)[axes[1],]^2 + assaycorrelations(eigensystem)[axes[2],]^2)
              polar.coords <- cbind(atan(assaycorrelations(eigensystem)[axes[2],]/assaycorrelations(eigensystem)[axes[1],])/pi,
                                    sqrt(cartesian.coords[,1]^2 + cartesian.coords[,2]^2))
              colnames(polar.coords) <- c("theta", "r")
              
              polar.coords[,"theta"] <- polar.coords[,"theta"] - zerophase
              polar.coords[scaled.coords[,1] < 0, "theta"] <- polar.coords[scaled.coords[,1] < 0, "theta"] + 1
              
              coordinates <- cbind(cartesian.coords, polar.coords)
              colnames(coordinates) <- c(paste0("eigenfeature", axes), "theta", "r")
              rownames(coordinates) <- colnames(matrix(eigensystem))
              
              as.data.frame(coordinates)
            }
          })
