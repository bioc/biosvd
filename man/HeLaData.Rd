\name{HeLaData_exp_DoubleThym_2}
\docType{data}
\alias{HeLaData}
\title{Human HeLa Cell Cycle Expression}
\usage{
  data(HeLaData_exp_DoubleThym_2)
  HeLaData
}
\description{
	Cells were arrested at the beginning of S phase by using a double thymidine block. Upon release from the thymidine block, cells were sampled every 1-2 hours for 44 hours during which the cells completed three cell cycles.
}
\format{
	An ExpressionSet with gene x sample expression data, gene annotation and sample information. Expression data comprise the un-logtransformed ratio of gene expression to reference mRNA from an asynchronous HeLa culture. For the samples cell cycle phase is known. For >850 genes that were identied by Whitfield et al to be periodically expressed during the cell cycle, the phase was determined based on correlation with genes known to be expressed in each cell cycle phase (e.g. cyclin E1 at the G1/S boundary, RAD51 in S phase, and TOP2A in G2).
}
\source{
	Whitfield ML, Sherlock G, Saldanha AJ, et al. Identification of genes periodically expressed in the human cell cycle and their expression in tumors. Mol Biol Cell 13, 1977-2000 (2002).
}
\keyword{datasets}
\author{
  Anneleen Daemen, Matthew Brauer
}
