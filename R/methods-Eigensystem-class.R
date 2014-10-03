### =========================================================================
### Eigensystem class methods 
### =========================================================================

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

##' @export
matrix <- function(x) slot(x, "matrix") 
##' @export
signMatrix <- function(x) slot(x, "signMatrix") 
##' @export
assayMatrix <- function(x) slot(x, "assayMatrix") 
##' @export
featureMatrix <- function(x) slot(x, "featureMatrix") 
##' @export
eigenassays <- function(x) slot(x, "eigenassays") 
##' @export
eigenexpressions <- function(x) slot(x, "eigenexpressions") 
##' @export
eigenfeatures <- function(x) slot(x, "eigenfeatures") 
##' @export
assaycorrelations <- function(x) slot(x, "assaycorrelations") 
##' @export
featurecorrelations <- function(x) slot(x, "featurecorrelations") 
##' @export
fractions <- function(x) slot(x, "fractions") 
##' @export
entropy <- function(x) slot(x, "entropy") 
##' @export
apply <- function(x) slot(x, "apply") 
##' @export
excludeEigenfeatures <- function(x) slot(x, "excludeEigenfeatures") 
##' @export
colorIdFeatures <- function(x) slot(x, "colorIdFeatures") 


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Setters
###

"matrix<-" <- function(x, value) 
{
  slot(x, "matrix") <- value
  x
}
"signMatrix<-" <- function(x, value) 
{
  slot(x, "signMatrix") <- value
  x
}
"assayMatrix<-" <- function(x, value) 
{
  slot(x, "assayMatrix") <- value
  x
}
"featureMatrix<-" <- function(x, value) 
{
  slot(x, "featureMatrix") <- value
  x
}
"eigenassays<-" <- function(x, value) 
{
  slot(x, "eigenassays") <- value
  x
}
"eigenexpressions<-" <- function(x, value) 
{
  slot(x, "eigenexpressions") <- value
  x
}
"eigenfeatures<-" <- function(x, value) 
{
  slot(x, "eigenfeatures") <- value
  x
}
"assaycorrelations<-" <- function(x, value) 
{
  slot(x, "assaycorrelations") <- value
  x
}
"featurecorrelations<-" <- function(x, value) 
{
  slot(x, "featurecorrelations") <- value
  x
}
"fractions<-" <- function(x, value) 
{
  slot(x, "fractions") <- value
  x
}
"entropy<-" <- function(x, value) 
{
  slot(x, "entropy") <- value
  x
}
"apply<-" <- function(x, value) 
{
  slot(x, "apply") <- value
  x
}
"excludeEigenfeatures<-" <- function(x, value) 
{
  slot(x, "excludeEigenfeatures") <- value
  x
}
"colorIdFeatures<-" <- function(x, value) 
{
  slot(x, "colorIdFeatures") <- value
  x
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

setMethod(show, "Eigensystem",
          function(object)
          {
            cat("class:", class(object), "\n")
            slots <- c("matrix", "signMatrix", "assayMatrix", "featureMatrix",
                       "eigenassays", "eigenexpressions", "eigenfeatures",
                       "assaycorrelations", "featurecorrelations", "fractions", "entropy",
                       "apply", "excludeEigenfeatures", "colorIdFeatures")
            for (i in slots) {
              elt <- slot(object, i)
              if (class(elt) %in% "matrix")
              	cat(i,": dimensions = ",dim(elt),"\n")
              else if (class(elt) %in% c("numeric","character","factor"))
              	cat(i,": ",elt,"\n")
             }
          }
          )



### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions
###


##' Transforms object input
##' 
##' Transforms object input
##' @title Transforms object input
##' @param object 
##' @return List with assay.matrix, feature.matrix, and matrix
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
.transformObject <- function(object) {

  assay.matrix <- feature.matrix <- matrix <- NULL
  if ("Eigensystem" %in% class(object)) {
    assay.matrix <- assayMatrix(object)
    feature.matrix <- featureMatrix(object)
    matrix <- matrix(object)
  } else if ("ExpressionSet" %in% class(object)) {
    assay.matrix <- as.matrix(pData(object))
    feature.matrix <- as.matrix(fData(object))
    matrix <- exprs(object)
  } else if (is.matrix(object)) {
    matrix <- object
  } else if (is.data.frame(object)) {
    matrix <- as.matrix(object)
  } else {
    stop("Input data needs to be of class matrix, data frame, ExpressionSet or eigensystem")
  }
  
  if (any(is.na(matrix))) {
    stop("Data matrix contains missing values")
  }
  
  return(list(assay.matrix, feature.matrix, matrix))
}


##' Transforms matrix
##' 
##' Transforms matrix to its variance in case apply equals variance
##' @title Transforms matrix
##' @param matrix
##' @param apply
##' @return matrix
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
.transformMatrix <- function(matrix, apply) {
  
  if (apply=="variance") {
    matrix[which(matrix==0)] <- 1e-10
    matrix <- log(matrix^2)
  }
  return(matrix)  
}


##' range check for eigenfeature axes
##' 
##' Checks that a vector of integers is between 1 and max number of eigenfeatures (inclusive) 
##' @title Checks validity of the eigenfeature xaxis and yaxis
##' @param eigensystem
##' @param axes
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
.checkEigenfeatureAxes <- function(eigensystem, axes) {
  neigenfeatures <- dim(eigenfeatures(eigensystem))[1]

  axesI <- as.integer(axes)
  if(any(axesI < 1) | any(axesI > neigenfeatures) | any(axes != axesI) | any(is.na(axesI))) 
    stop("Axes must be integers ranging between 1 and the total number of eigenfeatures inclusive")
  else
    axesI
}



##' Checks validity of color ID for assays or features
##' 
##' Checks validity of color ID for assays or features
##' @title Checks validity of color ID for assays or features
##' @param color.id
##' @param eigensystem
##' @param ids
##' @return color.id
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
.checkColorID <- function(color.id, eigensystem, ids=c("assays","features")) {
  
  ids <- match.arg(ids)
  if (ids=="assays") {
    direction <- "columns"
    eset <- "pheno"
    nbElements <- ncol(matrix(eigensystem))
    matrix <- assayMatrix(eigensystem)
  } else if (ids=="features") {
    direction <- "rows"
    eset <- "feature"
    nbElements <- nrow(matrix(eigensystem))
    matrix <- featureMatrix(eigensystem)
  }
  
  if (is.character(color.id) && length(color.id)==1) {
    if (!any(colnames(matrix)==color.id)) {
      stop(paste0("Color ID for the ",ids," does not match one of the variables in the ",eset," data of the original ExpressionSet data"))
    }
    color.id <- factor(matrix[,color.id])
  } else if (is.character(color.id)) {
    stop("Only one variable should be given for the color IDs")
  } else if (class(color.id) %in% c("numeric","factor") && length(color.id) != nbElements) {
    stop(paste0("Color IDs for the ",ids," do not match the number of ",direction," in the input data"))
  } else if (!class(color.id) %in% c("numeric","factor")) {
    stop(paste0("Color IDs for the ",ids," should be a vector of integers or a factor variable from the ",eset," data of the original ExpressionSet data"))
  }

  return(color.id)
}


##' Checks validity of contrast
##' 
##' Checks validity of contrast
##' @title Checks validity of contrast
##' @param contrast
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
.checkContrast <- function(contrast) {
  if (!is.numeric(contrast)) {
    stop("The contrast for heatmap visualization should be a positive integer")
  } else if (length(contrast)>1) {
    stop("Only one positive integer should be provided for the contrast for heatmap visualization")
  } else if (length(contrast)==1 && contrast < 0) {
    stop("The contrast for heatmap visualization should be a positive integer")
  } else {
    contrast
  }
}



##' Checks validity of prefix
##' 
##' Checks validity of prefix
##' @title Checks validity of prefix
##' @param prefix
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
.checkPrefix <- function(prefix) {
  if (!is.character(prefix)) {
    stop("The prefix for the plot and html file names should be a character")
  } else if (length(prefix)>1) {
    stop("Only one string should be provided for the prefix for the plot and html file names")
  } else {
    prefix
  }
}

##' Checks validity of eigenfeatures for exclusion
##' 
##' Checks validity of eigenfeatures for exclusion
##' @title Checks validity of eigenfeatures for exclusion
##' @param exclude.eigenfeatures
##' @param eigensystem object of class eigensystem
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
.checkExcludeEigenfeatures <- function(excludeEigenfeatures, eigensystem) {

  if (is.numeric(excludeEigenfeatures) && min(excludeEigenfeatures)<0) {
    stop("ExcludeEigenfeatures should be a vector of positive integers")
  } else if (is.numeric(excludeEigenfeatures) && min(excludeEigenfeatures)==0 && length(excludeEigenfeatures)>1) {
    stop("ExcludeEigenfeatures should be 0 when no eigenfeatures need to be excluded or a vector of positive integers ranging from 1 to the number of eigenfeatures in the eigensystem")
  } else if (is.numeric(excludeEigenfeatures) && max(excludeEigenfeatures)>ncol(matrix(eigensystem))) {
    stop("ExcludeEigenfeatures should be a vector of integers ranging from 1 to the number of eigenfeatures in the eigensystem")
  } else if (is.numeric(excludeEigenfeatures) && length(excludeEigenfeatures)>ncol(matrix(eigensystem))) {
    stop("ExcludeEigenfeatures should be a vector of maximum length equal to the number of eigenfeatures in the eigensystem")
  }
}


##' Checks dimensions of provided eigensystem and eigensystemPlotParam objects
##' 
##' Checks dimensions of provided eigensystem and eigensystemPlotParam objects
##' @title Checks dimensions of provided eigensystem and eigensystemPlotParam objects
##' @param eigensystem object of class Eigensystem
##' @param eigensystemPlotParams object of class EigensystemPlotParam
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
.checkDimensions <- function(eigensystem, eigensystemPlotParams) {
  nassays <- ncol(matrix(eigensystem))
  nfeatures <- nrow(matrix(eigensystem))
  whichAssays <- whichAssays(eigensystemPlotParams)
  whichFeatures <- whichFeatures(eigensystemPlotParams)
  whichEigenassays <- whichEigenassays(eigensystemPlotParams)
  whichEigenfeatures <- whichEigenfeatures(eigensystemPlotParams)
  if(length(whichFeatures) == 0) {
    whichFeatures <- seq(1, nfeatures) 
  }
  if(length(whichEigenfeatures) == 0) {
    whichEigenFeatures <- seq(1, nfeatures)
  }
  if(length(whichAssays) == 0) {
    whichAssays <- seq(1, nassays)
  }
  if(length(whichEigenassays) == 0) {
    whichEigenAssays <- seq(1, nassays)
  }
  if(length(whichFeatures) > nfeatures | length(whichEigenfeatures) > nfeatures |
     length(whichAssays) > nassays | length(whichEigenassays) > nassays) {
    stop("[EigensystemPlotParam: validation] too many features or assays being requested for plot.")
  } else {
    list(nassays=nassays, nfeatures=nfeatures,
         whichFeatures=whichFeatures, whichEigenfeatures=whichEigenfeatures,
         whichAssays=whichAssays, whichEigenassays=whichEigenassays,
         eigensystemrank=dim(eigenfeatures(eigensystem))[1])
  }
}


##' Checks validity of provided color map
##' 
##' Checks validity of provided color map
##' @title Checks validity of provided color map for assays and features
##' @param eigensystem object of class Eigensystem
##' @param eigensystemPlotParams object of class EigensystemPlotParam
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
.checkColorMaps <- function(eigensystem, eigensystemPlotParams) {
  # check or generate assayColorMap
  assayColorMap <- assayColorMap(eigensystemPlotParams)
  variables <- as.list(names(assayColorMap))
  assay.color.map <- lapply(variables,
                            function(variable) {
                              if(!variable %in% colnames(assayMatrix(eigensystem))) {
                                stop(paste("Color map defined for assay annotation that does not exist:",variable))
                              }
                              if(all(is.na(assayColorMap[[variable]]))) {
                                # create color map
                                assay.levels <- levels(as.factor(assayMatrix(eigensystem)[,variable]))
                                assay.colors <- rainbow(length(assay.levels))
                                names(assay.colors) <- assay.levels
                                assay.colors
                              } else { 
                                if(!any(names(assayColorMap[[variable]])  %in% assayMatrix(eigensystem)[,variable])) {
                                  stop(paste("Color map defined for assay annotations that do not exist:",
                                             variable,
                                             names(assayColorMap[[variable]][!names(assayColorMap[[variable]])  %in% assayMatrix(eigensystem)[,variable]])))
                                } else {
                                  assayColorMap[[variable]]
                                }
                              }
                            })
  names(assay.color.map) <- variables

  # check or generate featureColorMap
  featureColorMap <- featureColorMap(eigensystemPlotParams)
  variables <- as.list(names(featureColorMap))
  feature.color.map <- lapply(variables,
                              function(variable) {
                                if(!variable %in% colnames(featureMatrix(eigensystem))) {
                                  stop(paste("Color map defined for feature annotation that does not exist:",variable))
                                }
                                if(all(is.na(featureColorMap[[variable]]))) {
                                  feature.levels <- levels(as.factor(featureMatrix(eigensystem)[,variable]))
                                  feature.colors <- rainbow(length(feature.levels))
                                  names(feature.colors) <- feature.levels
                                  feature.colors
                                } else { 
                                  if(!any(names(featureColorMap[[variable]])  %in% featureMatrix(eigensystem)[,variable])) {
                                    stop(paste("Color map defined for feature annotations that do not exist:",
                                               variable,
                                               names(featureColorMap[[variable]][!names(featureColorMap[[variable]])  %in% featureMatrix(eigensystem)[,variable]])))
                                  } else {
                                    featureColorMap[[variable]]
                                  }
                                }
                              })
  names(feature.color.map) <- variables
  list(assayColorMap=assay.color.map, featureColorMap=feature.color.map)
}

##' Checks validity of figure
##' 
##' Checks validity of figure
##' @title Checks validity of figure
##' @param figure
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
.checkFigure <- function(figure) { 
  if (!is.logical(figure)) {
    stop("Variable figure should be a boolean (TRUE or FALSE) indicating whether figures should be shown or saved as pdf")
  } else {
    figure
  }
}

