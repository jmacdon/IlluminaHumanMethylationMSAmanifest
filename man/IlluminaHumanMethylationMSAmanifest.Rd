\name{IlluminaHumanMethylationMSAmanifest}
\alias{IlluminaHumanMethylationMSAmanifest}
\docType{data}
\title{
  Manifest object for the IlluminaHumanMethylationMSA array
}
\description{
  This object contains the array design for Illumina's Human Methylation
  MSA microarray.
}
\usage{data(IlluminaHumanMethylationMSAmanifest)}
\details{
  This package is based on the file
  \url{https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/infiniummethylationscreening/MSA-48v1-0_20102838_A1.csv}

  The script used to generate this package can be found in
  \code{system.files("scripts", "manifest.R", package = "IlluminaHumanMethylationMSAmanifest")}
}
\format{
  An object of class \code{IlluminaMethylationManifest}
}
\seealso{
  The description of the class in the minfi package,
  \code{\link[minfi]{IlluminaMethylationManifest-class}}.
}
\examples{
data(IlluminaHumanMethylationMSAmanifest)
}
\keyword{datasets}
