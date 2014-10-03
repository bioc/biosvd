setClassUnion("numericORfactor", c("numeric", "factor"))


##' The Eigensystem class
##' 
##' \code{Eigensystem} is a list-based class for storing the results of applying Singular Value Decomposition (SVD) to a feature by assay data set. Objects are normally created by \code{\link{compute,Eigensystem-method}}.
##' 
##' The \code{Eigensystem} class stores the original data and all SVD-derived information obtained with \code{compute}.
##' Data in the \code{Eigensystem} are organized into different slots,
##' \sQuote{matrix}, \sQuote{signMatrix}, \sQuote{assayMatrix}, \sQuote{featureMatrix},
##' \sQuote{eigenassays}, \sQuote{eigenexpressions}, \sQuote{eigenfeatures},
##' \sQuote{assaycorrelations}, \sQuote{featurecorrelations}, \sQuote{fractions},
##' \sQuote{entropy}, \sQuote{apply}, \sQuote{excludeEigenfeatures}, and \sQuote{colorIdFeatures}.
##' Brief descriptions of these slots are provided below.
##'
##' @name Eigensystem-class
##' @rdname Eigensystem-class
##' @aliases Eigensystem matrix signMatrix assayMatrix featureMatrix eigenassays eigenexpressions eigenfeatures assaycorrelations featurecorrelations fractions entropy apply excludeEigenfeatures colorIdFeatures show,Eigensystem-method
##' @exportClass Eigensystem 
##' @section Slots:
##' \code{Eigensystem} objects contain the following slots
##' \describe{
##' 	\item{\code{matrix}:}{\code{matrix} containing the feature by assay data without missing values}
##' 	\item{\code{signMatrix}:}{\code{matrix} containing the sign of each element in matrix}
##' 	\item{\code{assayMatrix}:}{\code{matrix} containing additional information about assays, with rows as assays and columns as additional variables}
##' 	\item{\code{featureMatrix}:}{\code{matrix} containing additional information about features, with rows as features and columns as additional variables}
##' 	\item{\code{eigenassays}:}{\code{matrix} containing the feature by eigenassay data, with each column in eigenassays corresponding to a left singular vector, representing genome-wide expression, proteome-wide abundance or metabolome-wide intensity in the corresponding eigenassay}
##' 	\item{\code{eigenexpressions}:}{\code{numeric} vector containing the eigenexpression fraction of each {eigenfeature, eigenassay}-pair, constituting the diagonal elements of the diagonal matrix connecting the left and right singular values; the diagonal matrix reflects the decoupling and decorrelation of the data, with expression of each eigenfeature restricted to the corresponding eigenassay}
##' 	\item{\code{eigenfeatures}:}{\code{matrix} containing the eigenfeatures by assay data, with each row corresponding to a right singular vector, representing the expression, abundance or intensity of the corresponding eigenfeature across all assays}
##' 	\item{\code{assaycorrelations}:}{\code{matrix} containing the correlation between the eigenassays as rows and the assays as columns}
##' 	\item{\code{featurecorrelations}:}{\code{matrix} containing the correlation between the eigenfeatures as rows and features as columns}
##' 	\item{\code{fractions}:}{\code{numeric} vector containing the eigenexpression fraction for each {eigenfeature, eigenassay}-pair, defined as the relative fraction of overall expression that each eigenfeature and eigenassay capture}
##' 	\item{\code{entropy}:}{\code{numeric} value between 0 and 1 giving the Shannon entropy as measure for data complexity, with an entropy of 0 corresponding to an ordered and redundant data set with all expression captured by a single {eigenfeature, eigenassay}-pair, and an entropy of 1 corresponding to a disordered and random data set with all {eigenfeature, eigenassay}-pairs equally expressed}
##' 	\item{\code{apply}:}{\code{character} containing whether the eigensystem should be computed for the actual data or the variance in the data}
##' 	\item{\code{excludeEigenfeatures}:}{\code{numeric} vector containing eigenfeature 1 and 2 in case they capture >85\% of the data with eigenfeature 2 capturing at least 15\%, otherwise numeric value containing eigenfeature 1}
##' 	\item{\code{colorIdFeatures}:}{\code{numeric} vector or \code{factor} containing annotation information on the features}
##' }
##'
##' @section Accessors:
##' \describe{
##'    \item{}{\code{matrix(x)}, \code{matrix(x) <- value}}
##'    \item{}{\code{signMatrix(x)}, \code{signMatrix(x) <- value}}
##'    \item{}{\code{assayMatrix(x)}, \code{assayMatrix(x) <- value}}
##'    \item{}{\code{featureMatrix(x)}, \code{featureMatrix(x) <- value}}
##'    \item{}{\code{eigenassays(x)}, \code{eigenassays(x) <- value}}
##'    \item{}{\code{eigenexpressions(x)}, \code{eigenexpressions(x) <- value}}
##'    \item{}{\code{eigenfeatures(x)}, \code{eigenfeatures(x) <- value}}
##'    \item{}{\code{assaycorrelations(x)}, \code{assaycorrelations(x) <- value}}
##'    \item{}{\code{featurecorrelations(x)}, \code{featurecorrelations(x) <- value}}
##'    \item{}{\code{fractions(x)}, \code{fractions(x) <- value}}
##'    \item{}{\code{entropy(x)}, \code{entropy(x) <- value}}
##'    \item{}{\code{apply(x)}, \code{apply(x) <- value}}
##'    \item{}{\code{excludeEigenfeatures(x)}, \code{excludeEigenfeatures(x) <- value}}
##'    \item{}{\code{colorIdFeatures(x)}, \code{colorIdFeatures(x) <- value}}
##'  }
##'
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
##' @examples
##' ## Metabolomics starvation data obtained from http://genomics-pubs.princeton.edu/StarvationMetabolomics/Download.shtml
##' data(StarvationData)
##'
##' ## An object from class Eigensystem is obtained with the compute method
##' eigensystem <- compute(StarvationData)
##' 
##' ## Obtain entropy
##' entropy(eigensystem)
##' @seealso \link{compute,Eigensystem-method}
##' @references Alter O, Brown PO and Botstein D. Singular value decomposition for genome-wide expression data processing and modeling. Proc Natl Acad Sci U.S.A. 97(18), 10101-10106 (2000).
##' @keywords methods classes
setClass(
         Class="Eigensystem",
         representation=representation(
           matrix = "matrix",
           signMatrix = "matrix",
           assayMatrix = "matrix",
           featureMatrix = "matrix",
           eigenassays = "matrix",
           eigenexpressions = "numeric",
           eigenfeatures = "matrix",
           assaycorrelations = "matrix",
           featurecorrelations = "matrix",
           fractions = "numeric",
           entropy = "numeric",
           apply = "character",
           excludeEigenfeatures = "numeric",
           colorIdFeatures = "numericORfactor"
           ),
         validity=function(object) {
         },
         prototype=prototype(
           matrix = matrix(),
           signMatrix = matrix(),
           assayMatrix = matrix(),
           featureMatrix = matrix(),
           eigenassays = matrix(),
           eigenexpressions = numeric(),
           eigenfeatures = matrix(),
           assaycorrelations = matrix(),
           featurecorrelations = matrix(),
           fractions = numeric(),
           entropy = numeric(),
           apply = character(),
           excludeEigenfeatures = numeric(),
           colorIdFeatures = numeric()
           )
         )

##' The EigensystemPlotParam class
##' 
##' \code{EigensystemPlotParam} is a list-based class for storing the parameters needed to specify plot features used by \code{link{plot,plot-method}} in plotting Eigensystem diagnostics, projections and transformations.
##' 
##' The \code{EigensystemPlotParam} class stores the list of desired plots, color palettes and keys, figure file and directory names and other necessary parameters.
##' Data in the \code{EigensystemPlotParam} may be organized into slots:
##' \sQuote{plots}, \sQuote{palette}, \sQuote{whichAssays}, \sQuote{whichFeatures},
##' \sQuote{whichEigenassays}, \sQuote{whichEigenfeatures}, \sQuote{whichPolarAxes},
##' \sQuote{assayColorMap}, \sQuote{featureColorMap}, \sQuote{contrast},
##' \sQuote{negativeValues}, \sQuote{path}, \sQuote{prefix},
##' \sQuote{filename}, \sQuote{figure},
##' brief descriptions of which follow.
##'
##' @name EigensystemPlotParam-class
##' @rdname EigensystemPlotParam-class
##' @aliases EigensystemPlotParam plots palette whichAssays whichFeatures whichEigenassays whichEigenfeatures whichPolarAxes assayColorMap
##' featureColorMap contrast negativeValues path prefix filenames figure show,EigensystemPlotParam-method
##' @exportClass EigensystemPlotParam 
##' @section Slots:
##' \code{EigensystemPlotParam} objects contain the following slots
##' \describe{
##' 	\item{\code{plots}:}{\code{character vector} indicating one or more plot choices from: "eigenfeatureHeatmap", "eigenassayHeatmap", "sortedHeatmap",
##'           "fraction","scree","zoomedFraction", "lines", "allLines", "eigenfeaturePolar", "eigenassayPolar". Defaults to all. }
##' 	\item{\code{palette}:}{\code{function} defining the palette to be used for heatmaps. Default is a Blue-Yellow color ramp. }
##' 	\item{\code{whichAssays}:}{\code{numeric vector} listing which assays are to be included in the plot(s). Default is all. }
##' 	\item{\code{whichFeatures}:}{\code{numeric vector} which features to include in plots. Default is all. }
##' 	\item{\code{whichEigenassays}:}{\code{numeric vector}  which eigenassays to include in plots. }
##' 	\item{\code{whichEigenfeatures}:}{\code{numeric vector} which eigenfeatures to include in plots. Defaults to first four. }
##'     \item{\code{whichPolarAxes}:}{\code{numeric vector} which two eigenassays/eigenfeatures to include in polar plots. Default is first two. }
##' 	\item{\code{assayColorMap}:}{ assayColorMap and featureColorMap are optional \code{list}s of colors corresponding to the
##'           levels of these annotations for assays and features. The elements of each list are named by the annotation column they
##'           correspond to. Each element is itself a named vector of colors, named by the level of the annotation it reflects (see example). }
##' 	\item{\code{featureColorMap}:}{\code{list} }
##' 	\item{\code{contrast}:}{\code{numeric} value specifying the contrast to use in heatmaps. }
##' 	\item{\code{negativeValues}:}{\code{logical} indicating whether scaling of values for heatmap should result in range that includes negative values. Default is TRUE.}
##' 	\item{\code{path}:}{\code{character} specifying the path of files for figure output. Default is current working directory.}
##' 	\item{\code{prefix}:}{\code{character} specifying an optional prefix to add to filenames. Default is biosvd.}
##' 	\item{\code{filenames}:}{\code{character} optional name for the files containing the plots. Names correspond to elements of 'plots' vector.}
##' 	\item{\code{figure}:}{\code{logical} specifying whether to ouptut plots into files. Default is FALSE. }
##' }
##' 
##' @section Accessors:
##' \describe{
##'    \item{}{\code{plots(x)}, \code{plots(x) <- value}}
##'    \item{}{\code{palette(x)}, \code{palette(x) <- value}}
##'    \item{}{\code{whichAssays(x)}, \code{whichAssays(x) <- value}}
##'    \item{}{\code{whichFeatures(x)}, \code{whichFeatures(x) <- value}}
##'    \item{}{\code{whichEigenassays(x)}, \code{whichEigenassays(x) <- value}}
##'    \item{}{\code{whichEigenfeatures(x)}, \code{whichEigenfeatures(x) <- value}}
##'    \item{}{\code{whichPolarAxes(x)}, \code{whichPolarAxes(x) <- value}}
##'    \item{}{\code{assayColorMap(x)}, \code{assayColorMap(x) <- value}}
##'    \item{}{\code{featureColorMap(x)}, \code{featureColorMap(x) <- value}}
##'    \item{}{\code{contrast(x)}, \code{contrast(x) <- value}}
##'    \item{}{\code{negativeValues(x)}, \code{negativeValues(x) <- value}}
##'    \item{}{\code{path(x)}, \code{path(x) <- value}}
##'    \item{}{\code{prefix(x)}, \code{prefix(x) <- value}}
##'    \item{}{\code{filenames(x)}, \code{filenames(x) <- value}}
##'    \item{}{\code{figure(x)}, \code{figure(x) <- value}}
##'  }
##'
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
##' @examples
##' data(YeastData_alpha)
##'
##' params <- new("EigensystemPlotParam")
##' cellcycle.col.map <- c("orange2", "darkgreen", "blue2", "magenta2", "red2")
##' names(cellcycle.col.map) <- c("S", "G2", "M", "M/G1", "G1")
##' assayColorMap(params) <- list(Cell.cycle.stage=cellcycle.col.map)
##' featureColorMap(params) <- list(Cell.cycle.stage=NA)
##' 
##' @references Alter O, Brown PO and Botstein D. Singular value decomposition for genome-wide expression data processing and modeling. Proc Natl Acad Sci U.S.A. 97(18), 10101-10106 (2000).
##' @keywords methods classes
setClass(
         Class="EigensystemPlotParam",
         representation=representation(
           plots="character",
           palette="function",
           whichAssays="numeric",
           whichFeatures="numeric",
           whichEigenassays="numeric",
           whichEigenfeatures="numeric",
           whichPolarAxes="numeric",
           assayColorMap="list",
           featureColorMap="list",
           contrast="numeric",
           negativeValues="logical",
           path="character",
           prefix="character",
           filenames="character",
           figure="logical"
           ),
         validity=function(object) {
           },
         prototype=prototype(
           plots=c("eigenfeatureHeatmap","eigenassayHeatmap","sortedHeatmap","fraction","scree","zoomedFraction","lines","allLines","eigenfeaturePolar","eigenassayPolar"),
           palette=colorRampPalette(c("blue","yellow"), space="rgb"),
           whichAssays=numeric(),
           whichFeatures=numeric(),
           whichEigenassays=numeric(),
           whichEigenfeatures=c(1:4),
           whichPolarAxes=c(2,1),
           assayColorMap=list(),
           featureColorMap=list(),
           contrast=3,
           negativeValues=TRUE,
           path=getwd(),
           prefix="biosvd",
           filenames=character(),
           figure=FALSE
           )
         )




