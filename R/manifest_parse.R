

read_custom_manifest=function(manifest_file=NULL,saveRDA=FALSE,name=NULL){

#init, checks & read file
  require(vroom)
  require(minfi)
  require(illuminaio)
  require(devtools)
  require(dplyr)

  if (is.null(manifest_file)){
    stop('manifest file option is empty, provide compatible .csv')
  }
  if(saveRDA&is.null(name)){
    stop('define the name for manifest .rda')
  }

#file <- "./custom_manifests/ExampleManifest-NA-NA-GRCh38.csv"
file=manifest_file

e1 <- vroom(file)

control.line <- grep("Controls",e1$Illumina)+1
assay.line <- grep("\\[Assay",e1$Illumina)+1

if(isEmpty(assay.line)){
  assay.line=0
  stop('Assay section is empty')
  }
if(isEmpty(control.line)){
  assay.line=0
  stop('Control line is empty. If control lines should be included in the manifest, check the input')
}

rm(e1)
gc()

manifest <- vroom(file,skip=assay.line,n_max =control.line-assay.line-2)

# check remaining duplicates & update empty values
duplicated_probes <- manifest$IlmnID[duplicated(manifest$IlmnID)]

remove_duplicated <- function(probe){
  sel <- which(manifest$Name == probe)
  if(any(manifest[sel,"Manifest_probe_match"])){
    sel_ID <- manifest$IlmnID[sel[which(manifest[sel,"Manifest_probe_match"]%>%pull())]][1]
  }else{
    sel_ID <- manifest$IlmnID[sel][1]
  }
  rmv_ID <- setdiff(manifest$IlmnID[sel],sel_ID)
  return(rmv_ID)
}

rmv_IlmnID <- unlist(sapply(duplicated_probes,remove_duplicated))

if(!is.null(rmv_IlmnID)){
manifest <- manifest[-which(manifest$IlmnID %in% rmv_IlmnID),]
}

if(any(duplicated(manifest$IlmnID))){
  warning('Found diplicated probes:',manifest$IlmnID[duplicated(manifest$IlmnID)])
}

manifest$AddressA_ID <- gsub("^0*", "", manifest$AddressA_ID)
manifest$AddressB_ID <- gsub("^0*", "", manifest$AddressB_ID)

manifest$AddressA_ID[is.na(manifest$AddressA_ID)] <- ""
manifest$AddressB_ID[is.na(manifest$AddressB_ID)] <- ""


#parse fields

TypeI <- manifest[manifest$Infinium_Design_Type == "I",
                  c("IlmnID", "AddressA_ID", "AddressB_ID", "Color_Channel", "Next_Base",
#                    c("Name", "AddressA_ID", "AddressB_ID", "Color_Channel", "Next_Base",
                    "AlleleA_ProbeSeq", "AlleleB_ProbeSeq")]

names(TypeI)[c(2, 3, 4, 5, 6, 7)] <- c("AddressA", "AddressB", "Color", "NextBase", "ProbeSeqA", "ProbeSeqB")

TypeI <- as(TypeI, "DataFrame")
TypeI$ProbeSeqA <- DNAStringSet(TypeI$ProbeSeqA)
TypeI$ProbeSeqB <- DNAStringSet(TypeI$ProbeSeqB)
TypeI$NextBase <- DNAStringSet(TypeI$NextBase)
TypeI$nCpG <- as.integer(
  oligonucleotideFrequency(TypeI$ProbeSeqB, width = 2)[, "CG"] - 1L)
TypeI$nCpG[TypeI$nCpG < 0] <- 0L

#Filter out nv & rs
TypeSnpI <- TypeI[grep("^rs", TypeI$IlmnID), ]
colnames(TypeSnpI)[1]='Name'
TypeI <- TypeI[!grepl("^rs", TypeI$IlmnID), ]
TypeI <- TypeI[!grepl("^nv", TypeI$IlmnID), ]
colnames(TypeI)[1]='Name'
#TypeSnpI <- TypeI[grep("^rs", TypeI$Name), ]
#TypeI <- TypeI[!grep("^rs", TypeI$Name), ]
#TypeI <- TypeI[!grep("^nv", TypeI$Name), ]

TypeII <- manifest[
  manifest$Infinium_Design_Type == "II",
  c("IlmnID", "AddressA_ID", "AlleleA_ProbeSeq")]
#  c("Name", "AddressA_ID", "AlleleA_ProbeSeq")]
names(TypeII)[c(2,3)] <- c("AddressA", "ProbeSeqA")
TypeII <- as(TypeII, "DataFrame")
TypeII$ProbeSeqA <- DNAStringSet(TypeII$ProbeSeqA)
TypeII$nCpG <- as.integer(letterFrequency(TypeII$ProbeSeqA, letters = "R"))
TypeII$nCpG[TypeII$nCpG < 0] <- 0L

#Filter out nv & rs
TypeSnpII <- TypeII[grep("^rs", TypeII$IlmnID), ]
colnames(TypeSnpII)[1]='Name'
TypeII <- TypeII[!grepl("^rs", TypeII$IlmnID), ]
TypeII <- TypeII[!grepl("^nv", TypeII$IlmnID), ]
colnames(TypeII)[1]='Name'
#TypeSnpII <- TypeII[grep("^rs", TypeII$Name), ]
#TypeII <- TypeII[-grep("^rs", TypeII$Name), ]
#TypeII <- TypeII[-grep("^nv", TypeII$Name), ]

controls <- read.table(
  file = file,
  skip = control.line,
  sep = ",",
  comment.char = "",
  quote = "",
  colClasses = c(rep("character", 5)))[, 1:5]

TypeControl <- controls[, 1:4]
names(TypeControl) <- c("Address", "Type", "Color", "ExtendedType")
TypeControl <- as(TypeControl, "DataFrame")

## assemble manifest & save
maniTmp <- list(
  manifestList = list(
    TypeI = TypeI,
    TypeII = TypeII,
    TypeControl = TypeControl,
    TypeSnpI = TypeSnpI,
    TypeSnpII = TypeSnpII),
  manifest = manifest,
  controls = controls)


maniList <- maniTmp$manifestList

IlluminaHumanMethylationCUSTOMmanifest <- IlluminaMethylationManifest(TypeI = maniList$TypeI,
                                                                      TypeII = maniList$TypeII,
                                                                      TypeControl = maniList$TypeControl,
                                                                      TypeSnpI = maniList$TypeSnpI,
                                                                      TypeSnpII = maniList$TypeSnpII,
                                                                      annotation = "IlluminaHumanMethylationCUSTOM")

return(IlluminaHumanMethylationCUSTOMmanifest)

if(saveRDA){
  save(IlluminaHumanMethylationCUSTOMmanifest,file=name)
}

}
