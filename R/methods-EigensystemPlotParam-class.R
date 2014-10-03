### =========================================================================
### EigensystemPlotParam class methods 
### =========================================================================

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###
##' @export
plots <- function(x) slot(x, "plots") 
##' @export
palette <- function(x) slot(x, "palette") 
##' @export
whichAssays <- function(x) slot(x, "whichAssays") 
##' @export
whichFeatures <- function(x) slot(x, "whichFeatures") 
##' @export
whichEigenassays <- function(x) slot(x, "whichEigenassays") 
##' @export
whichEigenfeatures <- function(x) slot(x, "whichEigenfeatures") 
##' @export
whichPolarAxes <- function(x) slot(x, "whichPolarAxes") 
##' @export
assayColorMap <- function(x) slot(x, "assayColorMap") 
##' @export
featureColorMap <- function(x) slot(x, "featureColorMap") 
##' @export
contrast <- function(x) slot(x, "contrast") 
##' @export
negativeValues <- function(x) slot(x, "negativeValues") 
##' @export
path <- function(x) slot(x, "path") 
##' @export
prefix <- function(x) slot(x, "prefix") 
##' @export
filenames <- function(x) slot(x, "filenames") 
##' @export
figure <- function(x) slot(x, "figure") 

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Setters
###

##' Set plots for object of EigensystemPlotParam
##' 
##' Set plots for object of EigensystemPlotParam
##' @title Set plots for object of EigensystemPlotParam
##' @param x object of class Eigensystem
##' @param value new value for slot plots
##' @return object of class Eigensystem
##' @export
##' @seealso EigensystemPlotParam-class
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
"plots<-" <- function(x, value)
{
  slot(x, "plots") <- value
  x
} 

##' Set palette for object of EigensystemPlotParam
##' 
##' Set palette for object of EigensystemPlotParam
##' @title Set palette for object of EigensystemPlotParam
##' @param x object of class Eigensystem
##' @param value new value for slot palette
##' @return object of class Eigensystem
##' @export
##' @seealso EigensystemPlotParam-class
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
"palette<-" <- function(x, value)
{
  slot(x, "palette") <- value
  x
} 

##' Set whichAssays for object of EigensystemPlotParam
##' 
##' Set whichAssays for object of EigensystemPlotParam
##' @title Set whichAssays for object of EigensystemPlotParam
##' @param x object of class Eigensystem
##' @param value new value for slot whichAssays
##' @return object of class Eigensystem
##' @export
##' @seealso EigensystemPlotParam-class
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
"whichAssays<-" <- function(x, value)
{
  slot(x, "whichAssays") <- value
  x
} 

##' Set whichFeatures for object of EigensystemPlotParam
##' 
##' Set whichFeatures for object of EigensystemPlotParam
##' @title Set whichFeatures for object of EigensystemPlotParam
##' @param x object of class Eigensystem
##' @param value new value for slot whichFeatures
##' @return object of class Eigensystem
##' @export
##' @seealso EigensystemPlotParam-class
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
"whichFeatures<-" <- function(x, value)
{
  slot(x, "whichFeatures") <- value
  x
} 

##' Set whichEigenassays for object of EigensystemPlotParam
##' 
##' Set whichEigenassays for object of EigensystemPlotParam
##' @title Set whichEigenassays for object of EigensystemPlotParam
##' @param x object of class Eigensystem
##' @param value new value for slot whichEigenassays
##' @return object of class Eigensystem
##' @export
##' @seealso EigensystemPlotParam-class
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
"whichEigenassays<-" <- function(x, value)
{
  slot(x, "whichEigenassays") <- value
  x
} 

##' Set whichEigenfeatures for object of EigensystemPlotParam
##' 
##' Set whichEigenfeatures for object of EigensystemPlotParam
##' @title Set whichEigenfeatures for object of EigensystemPlotParam
##' @param x object of class Eigensystem
##' @param value new value for slot whichEigenfeatures
##' @return object of class Eigensystem
##' @export
##' @seealso EigensystemPlotParam-class
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
"whichEigenfeatures<-" <- function(x, value)
{
  slot(x, "whichEigenfeatures") <- value
  x
} 

##' Set whichPolarAxes for object of EigensystemPlotParam
##' 
##' Set whichPolarAxes for object of EigensystemPlotParam
##' @title Set whichPolarAxes for object of EigensystemPlotParam
##' @param x object of class Eigensystem
##' @param value new value for slot whichPolarAxes
##' @return object of class Eigensystem
##' @export
##' @seealso EigensystemPlotParam-class
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
"whichPolarAxes<-" <- function(x, value) 
{
  slot(x, "whichPolarAxes") <- value
  x
}

##' Set assayColorMap for object of EigensystemPlotParam
##' 
##' Set assayColorMap for object of EigensystemPlotParam
##' @title Set assayColorMap for object of EigensystemPlotParam
##' @param x object of class Eigensystem
##' @param value new value for slot assayColorMap
##' @return object of class Eigensystem
##' @export
##' @seealso EigensystemPlotParam-class
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
"assayColorMap<-" <- function(x, value)
{
  slot(x, "assayColorMap") <- value
  x
} 

##' Set featureColorMap for object of EigensystemPlotParam
##' 
##' Set featureColorMap for object of EigensystemPlotParam
##' @title Set featureColorMap for object of EigensystemPlotParam
##' @param x object of class Eigensystem
##' @param value new value for slot featureColorMap
##' @return object of class Eigensystem
##' @export
##' @seealso EigensystemPlotParam-class
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
"featureColorMap<-" <- function(x, value)
{
  slot(x, "featureColorMap") <- value
  x
} 

##' Set contrast for object of EigensystemPlotParam
##' 
##' Set contrast for object of EigensystemPlotParam
##' @title Set contrast for object of EigensystemPlotParam
##' @param x object of class Eigensystem
##' @param value new value for slot contrast
##' @return object of class Eigensystem
##' @export
##' @seealso EigensystemPlotParam-class
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
"contrast<-" <- function(x, value)
{
  slot(x, "contrast") <- value
  x
} 

##' Set negativeValues for object of EigensystemPlotParam
##' 
##' Set negativeValues for object of EigensystemPlotParam
##' @title Set negativeValues for object of EigensystemPlotParam
##' @param x object of class Eigensystem
##' @param value new value for slot negativeValues
##' @return object of class Eigensystem
##' @export
##' @seealso EigensystemPlotParam-class
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
"negativeValues<-" <- function(x, value)
{
  slot(x, "negativeValues") <- value
  x
} 

##' Set path for object of EigensystemPlotParam
##' 
##' Set path for object of EigensystemPlotParam
##' @title Set path for object of EigensystemPlotParam
##' @param x object of class Eigensystem
##' @param value new value for slot path
##' @return object of class Eigensystem
##' @export
##' @seealso EigensystemPlotParam-class
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
"path<-" <- function(x, value)
{
  slot(x, "path") <- value
  x
} 

##' Set prefix for object of EigensystemPlotParam
##' 
##' Set prefix for object of EigensystemPlotParam
##' @title Set prefix for object of EigensystemPlotParam
##' @param x object of class Eigensystem
##' @param value new value for slot prefix
##' @return object of class Eigensystem
##' @export
##' @seealso EigensystemPlotParam-class
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
"prefix<-" <- function(x, value)
{
  slot(x, "prefix") <- value
  x
} 

##' Set filenames for object of EigensystemPlotParam
##' 
##' Set filenames for object of EigensystemPlotParam
##' @title Set filenames for object of EigensystemPlotParam
##' @param x object of class Eigensystem
##' @param value new value for slot filenames
##' @return object of class Eigensystem
##' @export
##' @seealso EigensystemPlotParam-class
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
"filenames<-" <- function(x, value)
{
  slot(x, "filenames") <- value
  x
} 

##' Set figure for object of EigensystemPlotParam
##' 
##' Set figure for object of EigensystemPlotParam
##' @title Set figure for object of EigensystemPlotParam
##' @param x object of class Eigensystem
##' @param value new value for slot figure
##' @return object of class Eigensystem
##' @export
##' @seealso EigensystemPlotParam-class
##' @author Anneleen Daemen \email{daemen.anneleen@@gene.com}, Matthew Brauer \email{brauer.matthew@@gene.com}
"figure<-" <- function(x, value)
{
  slot(x, "figure") <- value
  x
} 


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###
setMethod(show, "EigensystemPlotParam",
          function(object)
          {
            cat("class:", class(object), "\n")
            slots <- c("plots","palette","whichAssays",
                       "whichFeatures","whichEigenassays","whichEigenfeatures","whichPolarAxes",
                       "assayColorMap","featureColorMap","contrast",
                       "negativeValues","path","prefix",
                       "filenames","figure")
            for (i in slots) {
              elt <- slot(object, i)
              if(class(elt) == "list")
                cat(i,": ",unlist(elt),"\n")
              else if (!class(elt) %in% "function")
              	cat(i,": ",elt,"\n")
             }
          }
          )
