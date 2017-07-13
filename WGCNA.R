# based on WGCNA tutorial at https://labs.genetics.ucla.edu

library(WGCNA)
library(cluster)
library(ggplot2)
library(reshape2)
options(stringsAsFactors = FALSE)

datExprs <- t(collapsedExprs)
samples <- rownames(datExprs)

# create expression set with only early and late samples
early <- c("day 1")
late <- c("day 14")
earlyLate.eset <- getESet(filtered.eset,early,late)
groupA.eset = earlyLate.eset[,earlyLate.eset$code=="A"]
groupB.eset = earlyLate.eset[,earlyLate.eset$code=="B"]
timepointLabels <- sub("day 14","Late",groupB.eset$Day)
timepointLabels <- sub("day 1","Early",timepointLabels)
timepointLabels <- sub("day 7","Late",timepointLabels)

timepointLabels <- sub("day 14","Late",earlyLate.eset$Day)
timepointLabels <- sub("day 1","Early",timepointLabels)
timepointLabels <- sub("day 7","Late",timepointLabels)
timepointLabels <- sub("day 3","NA",timepointLabels)
timepointLabels <- sub("day 21","NA",timepointLabels)
timepointLabels <- sub("day 28","NA",timepointLabels)

# divide early and late samples
early.eset <- earlyLate.eset[,earlyLate.eset$Day=="day 1"]
late.eset <- earlyLate.eset[,earlyLate.eset$Day=="day 14"]
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
earlyExprs <- datExprs[c(earlyIndices_A,earlyIndices_B),]
earlyLateA <- datExprs[c(earlyIndices_A,lateIndices_A),]
earlyLateB <- datExprs[c(earlyIndices_B,lateIndices_B),]

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
removeOutliers <- function(exprsMatrix) {
  outlierIDs <- numeric(0)
  outliers <- TRUE
  # loop until no more outliers
  while (outliers) {
    sampleIDs <- rownames(exprsMatrix)
    A=adjacency(t(exprsMatrix),type="distance")
    # this calculates the whole network connectivity 
    k=as.numeric(apply(A,2,sum))-1
    # standardized connectivity 
    Z.k=scale(k)
    # Designate samples as outlying
    # if their Z.k value is below the threshold 
    thresholdZ.k=-2.5
    # the color vector indicates outlyingness (red) 
    outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
    # calculate the cluster tree using flashClust or hclust
    sampleTree = hclust(as.dist(1-A), method = "average")
    # Plot the sample dendrogram and the colors underneath. 
    plotDendroAndColors(sampleTree,groupLabels=names(outlierColor), colors=outlierColor,main="Sample dendrogram")

    outlierIndices <- which(outlierColor=="red")
    outlierIDs <- c(outlierIDs,sampleIDs[outlierIndices])
    exprsMatrix <- exprsMatrix[-outlierIndices,]
    if (length(outlierIndices)<1) 
      outliers=FALSE
  }
    return(outlierIDs)
}

outlierIDs_earlyA <- removeOutliers(earlyExprs_A)
outlierIDs_earlyB <- removeOutliers(earlyExprs_B) # none

outlierIDs_lateA <- removeOutliers(lateExprs_A)
outlierIDs_lateB <- removeOutliers(lateExprs_B)

outlierIDs_early <- removeOutliers(earlyExprs)
outlierIDs_late <- removeOutliers(lateExprs)

allOutliers <- removeOutliers(datExprs)
outlierIndices <- which((rownames(datExprs) %in% allOutliers))
datExprs <- datExprs[-outlierIndices,]


outlierIndices_earlyA <- which((rownames(earlyExprs_A) %in% outlierIDs_early))
outlierIndices_earlyB <- which((rownames(earlyExprs_B) %in% outlierIDs_early))
outlierIndices_lateA <- which((rownames(lateExprs_A) %in% outlierIDs_late))
outlierIndices_lateB <- which((rownames(lateExprs_B) %in% outlierIDs_late))

outlierIndices_earlyA <- which((rownames(earlyExprs_A) %in% outlierIDs_earlyA))
outlierIndices_lateA <- which((rownames(lateExprs_A) %in% outlierIDs_lateA))
outlierIndices_lateB <- which((rownames(lateExprs_B) %in% outlierIDs_lateB))

earlyExprs_A <- earlyExprs_A[-outlierIndices_earlyA,]
earlyExprs_B <- earlyExprs_B[-outlierIndices_B,]

lateExprs_A <- lateExprs_A[-outlierIndices_lateA,]
lateExprs_B <- lateExprs_B[-outlierIndices_lateB,]

# choose set of soft-thresholding powers
candidatePowers = c(c(1:10), seq(from = 12, to=32, by=2))

# call network topology analysis function
softThreshold = pickSoftThreshold(earlyExprs,networkType="signed",corFnc="bicor", powerVector = candidatePowers, verbose = 5)

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
                           verbose = 3)
return (bwnet)
}

### comparing expression of module eigengenes

net <- getModules(datExprs,12) 
moduleLabelsAutomatic=net$colors
# Convert labels to colors for plotting
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic)
# A data frame with module eigengenes can be obtained as follows 
MEsAutomatic=net$MEs


# subset into groups A and B and Early vs Late
allDays <- filtered.eset$Day
allCodes <- filtered.eset$code
B <- which(allCodes=="B")
A <- which(allCodes=="A")
early <- which(allDays=="day 1")
late <- which(allDays=="day 14")
inds <- intersect(A,c(early,late))
inds <- intersect(B,c(early,late))
MEsAutomatic <- MEsAutomatic[inds,]

groupA.eset <- baseline.eset[,baseline.eset$code=="A"]
groupA_indices <- match(colnames(exprs(groupA.eset)),samples)
datExprsA <- datExprs[groupA_indices,]
MEs_A <- MEsAutomatic[groupA_indices,]

groupB.eset <- filtered.eset[,filtered.eset$code=="B"]
groupB_indices <- match(colnames(exprs(groupB.eset)),samples)
datExprsB <- datExprs[groupB_indices,]
MEs_B <- MEsAutomatic[groupB_indices,]


MEsAutomatic$Timepoint <- factor(timepointLabels)
df <- melt(MEsAutomatic)
ggplot(data=df) + geom_boxplot(aes(x=Timepoint,y=value)) + facet_wrap(~variable,scales = "free") + ggtitle("Module Eigenegene Expression - Early vs. Late in Group A") + scale_fill_brewer(palette = "Accent")

library(genefilter)
# perform a t-test using the colttests function
ctt <- colttests(data.matrix(MEsAutomatic),MEsAutomatic$Timepoint,tstatOnly = FALSE)
barplot(ctt$p.value[-20],names.arg = rownames(ctt)[-20],ylim = c(0,1),main = "P-values - Early vs Late (Group B)",xlab = "Module Eigengenes",ylab = "P-values")
abline(h=0.05,col="red")


sizeGrWindow(10,7)
par(mfrow = c(1,2))
cex1 = 0.9
#this is the body weight
#weight = as.data.frame(datTraits$weight_g)
#names(weight)="weight"
# Next use this trait to define a gene significance variable 
GS.weight=as.numeric(cor(datExprFemale,weight,use="p"))
# This translates the numeric values into colors 
GS.weightColor=numbers2colors(GS.weight,signed=T)
blocknumber=3
datColors=data.frame(moduleColorsAutomatic)[net$blockGenes[[blocknumber]],]
# Plot the dendrogram and the module colors underneath 
plotDendroAndColors(net$dendrograms[[blocknumber]],colors=datColors, groupLabels=c("Module colors"),dendroLabels=FALSE, hang=0.03,addGuide=TRUE,guideHang=0.05)

