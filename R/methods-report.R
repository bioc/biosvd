##' Generates a summarizing report of the eigensystem
##'
##' The function generates a txt file of the eigensystem, containing the list of features with their coordinates, radius and phase in the polar plot according to two selected eigenfeatures.
##' 
##' @title Creates a report of the eigensystem
##' @param x object of class Eigensystem
##' @param y object of class EigensystemPlotParam
##' @param ... Additional arguments that can be passed on
##' @return NULL
##' @importFrom methods callGeneric show
##' @export
##' @docType methods
##' @rdname report-methods
##' @aliases report report,Eigensystem,EigensystemPlotParam-method
##' @section Methods:
##' \describe{
##' \item{\code{signature(x = "Eigensystem", y = "EigensystemPlotParam")}}{
##' }
##' }
##' @family "eigensystem"
##' @examples
##' ## Metabolomics starvation data obtained from http://genomics-pubs.princeton.edu/StarvationMetabolomics/Download.shtml
##' data(StarvationData)
##'
##' ## Computes the eigensystem for the actual data
##' eigensystem <- compute(StarvationData)
##' ## Exclude the eigenfeatures representing steady-state intensity
##' eigensystem <- exclude(eigensystem)
##' ## Computes the eigensystem on the variance in the data after filtering out stead-state intensity
##' eigensystem <- compute(eigensystem, apply="variance")
##' ## No exclusion of eigenfeatures representing steady-scale variance
##' eigensystem <- exclude(eigensystem, excludeEigenfeatures=0)
##'
##' ## Generate report for eigenfeatures 1 and 2
##' params <- new("EigensystemPlotParam")
##' if (.Platform$OS.type %in% "windows") path(params) <- getwd()
##' report(eigensystem, params)
##' ## Generate report for eigenfeatures 2 and 3
##' whichPolarAxes(params) <- c(3,2)
##' report(eigensystem, params)
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
setMethod("report", signature(x="Eigensystem", y="EigensystemPlotParam"),
          function(x, y, ...) {
            eigensystem <- x
            params <- y
            polarAxes <- .checkEigenfeatureAxes(eigensystem, whichPolarAxes(params))
            
            polarplot.info <- cbind(rownames(matrix(eigensystem)),project(eigensystem, polarAxes))
            colnames(polarplot.info)=c("Feature",paste("Coordinate to eigenfeature ",polarAxes[2],sep=""),paste("Coordinate to eigenfeature ",polarAxes[1],sep=""),"Radius","Phase")
            if (class(featureMatrix(eigensystem))!="NULL" && ncol(featureMatrix(eigensystem)) > 0) polarplot.info <- cbind(polarplot.info, featureMatrix(eigensystem))
            
            write.table(polarplot.info, file=paste(prefix(params),"report.eigenfeature",polarAxes[1],"vs",polarAxes[2],"txt",sep="."), row.names=FALSE)
          })

