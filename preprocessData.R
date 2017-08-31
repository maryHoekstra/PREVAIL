# correct for batch effects and filter 
# read in annotation data 

load("prevailAug_eset.rda")
prevail.eset <- prevailJul3.eset[,1:188] # last 21 samples are healthy controls

# get patient codes
groupDesignation <- read.csv("PREVAIL_randomization_20170303.csv")
IDs <- read.csv("PREVAIL_patient_20170303.csv")
patientIDs <- prevail.eset$patient_id
seqID_indices <- match(patientIDs,IDs$patient_id)
patientSeqIDs <- as.character(IDs[seqID_indices,"seq_id"])
patientCodeIndices <- match(patientSeqIDs,groupDesignation$seq_id)
patientCodes <- as.character(groupDesignation[patientCodeIndices,"code"])

filtered.eset$code <- patientCodes

# load packages
library(Biobase)
library(genefilter)
library(sva)
library(AnnotationDbi)
library(dplyr)

# correct for batch effects
phenoData <- pData(prevail.eset)
batch <- phenoData$Batch
modcombat <- model.matrix(~1, data=phenoData) # no adjustment variables, just fit an intercept term
exprsData <- exprs(prevail.eset)
correctedExprs <- ComBat(dat=exprsData, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

# create new corrected expression set
corrected.eset <- prevail.eset
exprs(corrected.eset) <- correctedExprs

# take the top 50% most variable probes
# in many tissues, only 40% of genes are expressed
# IQR is robust to outliers 
filtered.eset <- varFilter(corrected.eset, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)
filteredExprs <- exprs(filtered.eset)

# read in annotation data
getAnnotation <- function() {
  geo_sheet <- read.delim('/Users/maryhoekstra/GPL15207-14509.txt', header=TRUE, sep='\t', skip=24)
  probes <- dplyr::select(geo_sheet, ID, Gene.Title, Gene.Symbol, Chromosomal.Location, Entrez.Gene, SwissProt, OMIM)
  return(probes)
}
annoData <- getAnnotation()




