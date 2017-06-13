# based on WGCNA tutorial at https://labs.genetics.ucla.edu

options(stringsAsFactors = FALSE);

# transpose expression data so rows are samples and columns are genes
datExprs <- t(collapsedExprs)

# choose set of soft-thresholding powers
candidatePowers = c(c(1:10), seq(from = 12, to=32, by=2))

# call network topology analysis function
softThreshold = pickSoftThreshold(datExprs,networkType="signed",corFnc="bicor", powerVector = candidatePowers, verbose = 5)

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
bwnet = blockwiseModules(datExprs, maxBlockSize = 6074,
                         power = 12, networkType= "signed", TOMType = "signed", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE,
                         corType = "bicor",
                         saveTOMs = TRUE,
                         saveTOMFileBase = "TOM-blockwise",
                         verbose = 3)

# analysis of modules 
# define numbers of genes and samples
nGenes = ncol(datExprs)
nSamples = nrow(datExprs)

# recalculate MEs with color labels
MEs0 = moduleEigengenes(datExprs, moduleColours)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, phenoData, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
