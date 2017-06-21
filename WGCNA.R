# based on WGCNA tutorial at https://labs.genetics.ucla.edu

library(WGCNA)
options(stringsAsFactors = FALSE);

# transpose expression data so rows are samples and columns are genes
datExprs <- t(collapsedExprs)

## try subsetting the matrix and see about the scale free topology
# use all samples from group A
groupA.eset <- filtered.eset[ , filtered.eset$code=="A"]
groupA_patientIDs <- colnames(exprs(groupA.eset))
groupA_sampleIndices <- match(groupA_patientIDs,rownames(datExprs))
datExprs_A <- datExprs[groupA_sampleIndices,]

# use all samples from group B
groupB.eset <- filtered.eset[ , filtered.eset$code=="B"]
groupB_patientIDs <- colnames(exprs(groupB.eset))
groupB_sampleIndices <- match(groupB_patientIDs,rownames(datExprs))
datExprs_B <- datExprs[groupB_sampleIndices,]

# use only samples at day 3 for both groups
day3.eset <- filtered.eset[,filtered.eset$Day=="day 3"]
day3_samples <- colnames(exprs(day3.eset))
day3_sampleIndices <- match(day3_samples,rownames(datExprs))
datExprs_day3 <- datExprs[day3_sampleIndices,]

# use only samples from day 1 in group B
day1B.eset <- groupB.eset[,groupB.eset$Day=="day 1"]
day1B_samples <- colnames(exprs(day1B.eset))
day1B_sampleIndices <- match(day1B_samples,rownames(datExprs))
datExprs_day1B <- datExprs[day1B_sampleIndices,]

# use only samples from day 7 in group B
day14B.eset <- groupB.eset[,groupB.eset$Day=="day 14"]
day14B_samples <- colnames(exprs(day14B.eset))
day14B_sampleIndices <- match(day14B_samples,rownames(datExprs))
datExprs_day14B <- datExprs[day14B_sampleIndices,]

# choose set of soft-thresholding powers
candidatePowers = c(c(1:10), seq(from = 12, to=32, by=2))

# call network topology analysis function
softThreshold = pickSoftThreshold(datExprs_A,networkType="signed",corFnc="bicor", powerVector = candidatePowers, verbose = 5)

# plot results
sizeGrWindow(10,7)
par(mfrow = c(1,2))
cex1 = 0.9

# scale-free topology fit index as a function of the soft-thresholding power
plot(softThreshold$fitIndices[,1], -sign(softThreshold$fitIndices[,3])*softThreshold$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(softThreshold$fitIndices[,1], -sign(softThreshold$fitIndices[,3])*softThreshold$fitIndices[,2],
     labels=candidatePowers,cex=cex1,col="red")

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# mean connectivity as a function of the soft-thresholding power
plot(softThreshold$fitIndices[,1], softThreshold$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(softThreshold$fitIndices[,1], softThreshold$fitIndices[,5], labels=candidatePowers, cex=cex1,col="red")

# do blockwise analysis since we have more than 8000 probes on a 4GB laptop
getModules <- function(datExprs,sfPower) {
  bwnet = blockwiseModules(datExprs, maxBlockSize = 6074,
                           power = sfPower, networkType= "signed", TOMType = "signed", minModuleSize = 30,
                           reassignThreshold = 0, mergeCutHeight = 0.25,
                           numericLabels = TRUE,
                           corType = "bicor",
                           saveTOMs = TRUE,
                           saveTOMFileBase = "TOM-blockwise",
                           verbose = 3)
return (bwnet)
}

bwnetA <- getModules(datExprs_A,12)
