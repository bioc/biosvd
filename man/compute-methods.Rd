\docType{methods}
\name{compute,Eigensystem-method}
\alias{compute}
\alias{compute,Eigensystem-method}
\alias{compute,ExpressionSet-method}
\alias{compute,data.frame-method}
\alias{compute,matrix-method}
\title{Compute the eigensystem for a feature by assay matrix}
\usage{
\S4method{compute}{Eigensystem}(object, apply = c("data", "variance"))
}
\arguments{
  \item{object}{object of class matrix, data.frame,
  Expressionset or Eigensystem, containing the feature x
  assay expression or intensity data}

  \item{apply}{the actual expression/intensity data (data)
  or variance in the data (variance) to which to apply the
  function (default = data)}
}
\value{
Object of class Eigensystem
}
\description{
Computes the eigensystem for a feature by assay matrix by
applying Singular Value Decomposition.
}
\details{
Function compute decomposes the input data set from the
feature x assay space to the reduced diagonalized
"eigenfeatures x eigenassays" space, with the eigenfeatures
and eigenassays unique orthonormal superpositions of the
features and assays, respectively. This approach allows
filtering out eigenfeatures and eigenassays that are
inferred to represent noise or experimental artifacts,
either at the expression/intensity level or the variance
level. The function can be applied to an object of class
matrix, data.frame, ExpressionSet, or eigensystem.
}
\section{Methods}{
  \describe{ \item{\code{signature(x = "Eigensystem")}}{ }
  \item{\code{signature(x = "data.frame")}}{ }
  \item{\code{signature(x = "matrix")}}{ }
  \item{\code{signature(x = "ExpressionSet")}}{ } }
}
\examples{
## Metabolomics starvation data obtained from http://genomics-pubs.princeton.edu/StarvationMetabolomics/Download.shtml
data(StarvationData)

## Computes the eigensystem for the actual expression/intensity data
eigensystem <- compute(StarvationData)
## Computes the eigensystem for the variance in the data
eigensystem <- compute(StarvationData, apply="variance")
}
\author{
Anneleen Daemen \email{daemen.anneleen@gene.com}, Matthew
Brauer \email{brauer.matthew@gene.com}
}
\references{
Alter O, Brown PO and Botstein D. Singular value
decomposition for genome-wide expression data processing
and modeling. Proc Natl Acad Sci U.S.A. 97(18), 10101-10106
(2000).
}
\seealso{
Eigensystem-class

Other X.eigensystem.: \code{\link{exclude}},
\code{\link{exclude,Eigensystem-method}},
\code{\link{exclude,Eigensystem-method}};
\code{\link{plot}},
\code{\link{plot,Eigensystem,ANY-method}},
\code{\link{plot,Eigensystem,ANY-method}};
\code{\link{report}},
\code{\link{report,Eigensystem-method}},
\code{\link{report,Eigensystem-method}}
}

