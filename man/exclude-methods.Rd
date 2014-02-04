\docType{methods}
\name{exclude,Eigensystem-method}
\alias{exclude}
\alias{exclude,Eigensystem-method}
\title{Excludes specified eigenfeatures/eigenassays from the data}
\usage{
\S4method{exclude}{Eigensystem}(x, excludeEigenfeatures = NULL)
}
\arguments{
  \item{x}{object of class eigensystem}

  \item{excludeEigenfeatures}{vector of positive integers
  representing eigenfeatures to be excluded from the data}
}
\value{
Object of class eigensystem
}
\description{
Excludes specified eigenfeatures/eigenassays from the
original data
}
\details{
The function excludes eigenfeatures/eigenassays from the
data that correspond to steady-state intensity,
steady-scale variance, experimental artifacts and/or noise
as specified by the user. In case no eigenfeatures are
specified, the eigenfeature(s) corresponding to
steady-state/steady-scale is/are removed. Filtering out
steady-state expression/intensity corresponds to centering
the expression/intensity patterns at steady-state
expression/intensity level (arithmetic mean of
expression/intensity ~ 0). Filtering out steady-scale
variance corresponds to normalization by the steady scale
of expression/intensity variance (geometric mean of
variance ~ 1).
}
\section{Methods}{
  \describe{ \item{\code{signature(x = "Eigensystem")}}{ }
  }
}
\examples{
## Metabolomics starvation data obtained from http://genomics-pubs.princeton.edu/StarvationMetabolomics/Download.shtml
data(StarvationData)

## Computes the eigensystem for the actual data
eigensystem <- compute(StarvationData)
## Excludes the eigenfeature(s) representing steady-state expression/intensity as defined in compute
exclude(eigensystem)
## Excludes user-specified eigenfeatures 1, 4 and 5
exclude(eigensystem, excludeEigenfeatures=c(1,4,5))

## Computes the eigensystem for the variance in the data
eigensystem <- compute(StarvationData, apply="variance")
## Excludes the eigenfeature(s) representing steady-scale variance as defined in compute
eigensystem <- exclude(eigensystem)
## Excludes none of the eigenfeatures and recalculates the eigensystem for the actual data
eigensystem <- exclude(eigensystem, excludeEigenfeatures=0)
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
Other X.eigensystem.: \code{\link{compute}},
\code{\link{compute,data.frame-method}},
\code{\link{compute,Eigensystem-method}},
\code{\link{compute,Eigensystem-method}},
\code{\link{compute,ExpressionSet-method}},
\code{\link{compute,matrix-method}}; \code{\link{plot}},
\code{\link{plot,Eigensystem,ANY-method}},
\code{\link{plot,Eigensystem,ANY-method}};
\code{\link{report}},
\code{\link{report,Eigensystem-method}},
\code{\link{report,Eigensystem-method}}
}

