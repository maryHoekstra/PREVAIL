# based on WGCNA tutorial at https://labs.genetics.ucla.edu

library(WGCNA)
library(cluster)
options(stringsAsFactors = FALSE)

datExprs <- t(collapsedExprs)
samples <- rownames(datExprs)

# create expression set with only early and late samples
early <- c("day 1")
late <- c("day 14", "day 21", "day 28")
earlyLate.eset <- getESet(filtered.eset,early,late)

# divide early and late samples
early.eset <- earlyLate.eset[,earlyLate.eset$Day=="day 1"]
late.eset <- earlyLate.eset[,earlyLate.eset$Day=="day 14" | earlyLate.eset$Day=="day 21" | earlyLate.eset$Day=="day 28"]

# divide groups A and B
early_groupA.eset <- early.eset[ , early.eset$code=="A"]
late_groupA.eset <- late.eset[, late.eset$code=="A"]
early_groupB.eset <- early.eset[, early.eset$code=="B"]
late_groupB.eset <- late.eset[, late.eset$code=="B"]

earlyIndices_A <- match(colnames(exprs(early_groupA.eset)),samples)
earlyExprs_A <- datExprs[earlyIndices_A,]
lateIndices_A <- match(colnames(exprs(late_groupA.eset)),samples)
lateExprs_A <- datExprs[lateIndices_A,]
earlyIndices_B <- match(colnames(exprs(early_groupB.eset)),samples)
earlyExprs_B <- datExprs[earlyIndices_B,]
lateIndices_B <- match(colnames(exprs(late_groupB.eset)),samples)
lateExprs_B <- datExprs[lateIndices_B,]

lateExprs <- datExprs[c(lateIndices_A,lateIndices_B),]

# one method for removing outliers from groups A and B
sampleTree_late = hclust(dist(lateExprs), method = "average")
plot(sampleTree_late, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 65, col = "red");
clust = cutreeStatic(sampleTree_late, cutHeight = 65, minSize = 10)
table(clust)
keepSamples = (clust==1)
lateExprs = lateExprs[keepSamples, ]

# another method for removing outliers based on squared Euclidean distance 
# note that we transpose the data 
A=adjacency(t(lateExprs_A),type="distance")
# this calculates the whole network connectivity 
k=as.numeric(apply(A,2,sum))-1
# standardized connectivity 
Z.k=scale(k)
# Designate samples as outlying
# if their Z.k value is below the threshold 
thresholdZ.k=-2.5 # often -2.5
# the color vector indicates outlyingness (red) 
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
# calculate the cluster tree using flashClust or hclust
sampleTree = hclust(as.dist(1-A), method = "average")
# Plot the sample dendrogram and the colors underneath. 
plotDendroAndColors(sampleTree,groupLabels=names(outlierColor), colors=outlierColor,main="Sample dendrogram")

# remove sample 270755 in B
lateExprs_B <- lateExprs_B[-12,]
# remove sample 270820 from A
lateExprs_A <- lateExprs_A[-12,]

# choose set of soft-thresholding powers
candidatePowers = c(c(1:10), seq(from = 12, to=32, by=2))

# call network topology analysis function
softThreshold = pickSoftThreshold(lateExprs_B,networkType="signed",corFnc="bicor", powerVector = candidatePowers, verbose = 5)

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

bwnetA <- getModules(datExprs = earlyExprs_A,sfPower = 26)
bwnetB <- getModules(datExprs = earlyExprs_B,sfPower = 26)
