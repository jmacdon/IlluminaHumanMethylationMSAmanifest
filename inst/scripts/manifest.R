## run from home dir
setwd("~/")
library(minfi)
wheretoput <- "~/ogopogo/data3/Rpacks/IlluminaHumanMethylationMSAmanifest/data/IlluminaHumanMethylationMSAmanifest.rda"

stopifnot(file.exists("MSA-48v1-0_20102838_A1.csv"))

z <- minfi:::read.manifest.EPIC("MSA-48v1-0_20102838_A1.csv")
zz <- z$manifestList

IlluminaHumanMethylationMSAmanifest <- IlluminaMethylationManifest(TypeI = zz$TypeI,
                                                                   TypeII = zz$TypeII,
                                                                   TypeControl = zz$TypeControl,
                                                                   TypeSnpI = zz$TypeSnpI,
                                                                   TypeSnpII = zz$TypeSnpII,
                                                                   annotation = "IlluminaHumanMethylationMSA")

stopifnot(validObject(IlluminaHumanMethylationMSAmanifest))
save(IlluminaHumanMethylationMSAmanifest, compress = "xz",
     file = wheretoput)
