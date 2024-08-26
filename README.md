# # Workaround to make MSA Infinium DNA methylation array compatible with minfi package

## Install from GitHub or from .tar.gz


```
library(devtools)
install_github('https://github.com/Illumina/IlluminaHumanMethylationMSAmanifest')
install.packages('IlluminaHumanMethylationMSAmanifest_0.1.0.tar.gz',repos=NULL,type='source')
```

## Run

```
#read intensities
MSA_example = read.metharray.exp(base='MSA_IDATs')

#assign annotation
annotation(MSA_example)["array"] = "IlluminaHumanMethylationMSA"
annotation(MSA_example)["annotation"] = "20a1.hg38"

#read beta into matrix
betas_msa_example= as.data.frame(getBeta(MSA_example))
```
