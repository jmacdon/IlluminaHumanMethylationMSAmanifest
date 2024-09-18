# package for MSA Infinium DNA Methylation array and minfi compatibility


```
library(BiocManager)
install_github("jmacdon/IlluminaHumanMethylationMSAmanifest")

```

## Run

```
#read intensities
MSA_example = read.metharray.exp(...)

#Annotation is assigned by hand presently
annotation(MSA_example) <- list(array = "IlluminaHumanMethylationMSA",
	                            annotation = "20a1.hg38")

```
