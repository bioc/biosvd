\name{YeastData_alpha}
\docType{data}
\alias{YeastData}
\title{Yeast Cell Cycle Expression}
\usage{
  data(YeastData_alpha)
  YeastData
}
\description{
To catalog genes in Saccharomyces cerevisiae whose transcript levels vary periodically within the cell cycle, mRNA levels in samples from yeast cultures were synchronized in G1 phase with alpha factor arrest. After release of the alpha factor, cells were sampled every 7 minutes over a timespan of 140 minutes, during which the cells synchronously completed two cell cycles.
}
\format{
An ExpressionSet with gene x sample expression data, gene annotation and sample information. Expression data comprise the un-logtransformed ratio of gene expression to reference mRNA from an asynchronous yeast culture. For the samples cell cycle phase is known. For 800 cell cycle-regulated genes, the phase in which these genes reach their peak expression was determined based on published timing of the expression of known cell cycle-regulated genes.
}
\source{
Spellman PT, Sherlock G, Zhang MQ, et al. Comprehensive identification of cell cycle-regulated genes of the Yeast Saccharomyces cerevisiae by microarray hybridization. Mol Biol Cell 9, 3273-3297 (1998).
}
\references{
  Alter O, Brown PO and Botstein D. Singular value
  decomposition for genome-wide expression data processing
  and modeling. Proc Natl Acad Sci U.S.A. 97(18),
  10101-10106 (2000).
}
\keyword{datasets}
\author{
  Anneleen Daemen, Matthew Brauer
}
