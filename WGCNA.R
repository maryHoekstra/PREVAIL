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
softThreshold = pickSoftThreshold(datExprs_day14B,networkType="signed",corFnc="bicor", powerVector = candidatePowers, verbose = 5)

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

##

# do blockwise analysis since we have more than 8000 probes on a 4GB laptop
bwnet = blockwiseModules(datExprs, maxBlockSize = 6074,
                         power = 30, networkType= "signed", TOMType = "signed", minModuleSize = 30,
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

# write gene list for each module to a file 
moduleColours <- bwnet$colors
genes <- colnames(datExprs)

for (i in 1:max(moduleColours_30)) {
  moduleGenes <- (moduleColours_30==i)
  moduleEntrezIDs <- genes[moduleGenes]
  write(moduleEntrezIDs,file=paste("/Users/maryhoekstra/Desktop/module",i,"_30.txt"))
}
  
# recalculate MEs with color labels
MEs0 = moduleEigengenes(datExprs, moduleColours)$eigengenes
MEs = orderMEs(MEs0)
# examine correlations between modules and clinical traits
# look at code (group) and day
moduleTraitCor = cor(MEs, phenoData, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(selectPhenoData),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

## consensus analysis
# separate data from group A and B
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

# we work with two sets:
nSets = 2
# for easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("Group A", "Group B")
shortLabels = c("A", "B")
# form multi-set expression data: columns starting from 9 contain actual expression data.
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = as.data.frame(datExprs_A))
#names(multiExpr[[1]]$data) = femData$substanceBXH;
#rownames(multiExpr[[1]]$data) = names(femData)[-c(1:8)];
multiExpr[[2]] = list(data = as.data.frame(datExprs_B))
#names(multiExpr[[2]]$data) = maleData$substanceBXH;
#rownames(multiExpr[[2]]$data) = names(maleData)[-c(1:8)];
# Check that the data has the correct format for many functions operating on multiple sets:
exprSize = checkSets(multiExpr)

# check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK

# cluster samples on euclidean distance
sampleTrees = list()
for (set in 1:nSets)
{
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}

bnet = blockwiseConsensusModules(
  multiExpr, maxBlockSize = 6074, power = 12, minModuleSize = 30,
  deepSplit = 2, networkType= "signed", TOMType = "signed",
  corType = "bicor",
  pamRespectsDendro = FALSE,
  mergeCutHeight = 0.25, numericLabels = TRUE,
  minKMEtoStay = 0,
  saveTOMs = TRUE, verbose = 5)

moduleLabels <- bnet$colors
moduleColours = labels2colors(moduleLabels)
# recalculate consMEs to give them color names
consMEsC = multiSetMEs(multiExpr, universalColors = moduleColours)
# order eigengenes by consesus hierarchical clustering:
MET = consensusOrderMEs(consMEsC)

# perform differential analysis
sizeGrWindow(8,10);
#pdf(file = "Plots/EigengeneNetworks.pdf", width= 8, height = 10);
par(cex = 0.9)
plotEigengeneNetworks(MET, setLabels, marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1),
                      zlimPreservation = c(0.5, 1), xLabelsAngle = 90)
dev.off()

# choose the least preserved modules by hand?
magentaModuleGenes <- (moduleColours=="magenta")
magentaModuleIDs <- genes[magentaModuleGenes]
write(magentaModuleIDs,file=paste("/Users/maryhoekstra/Desktop/magentaModule.txt"))
midnightBlueModuleGenes <- (moduleColours=="midnightblue")
midnightBlueModuleIDs <- genes[midnightBlueModuleGenes]
write(midnightBlueModuleIDs,file=paste("/Users/maryhoekstra/Desktop/midnightBlueModule.txt"))

# try module preservation analysis between group A and B
# the gene expression data: columns are genes, rows are arrays (samples)
# separate data from group A and B
groupA.eset <- filtered.eset[ , filtered.eset$code=="A"]
groupA_patientIDs <- colnames(exprs(groupA.eset))
groupA_sampleIndices <- match(groupA_patientIDs,rownames(datExprs))
datExprs_A <- datExprs[groupA_sampleIndices,]

# use all samples from group B
groupB.eset <- filtered.eset[ , filtered.eset$code=="B"]
groupB_patientIDs <- colnames(exprs(groupB.eset))
groupB_sampleIndices <- match(groupB_patientIDs,rownames(datExprs))
datExprs_B <- datExprs[groupB_sampleIndices,]

# do blockwise analysis since we have more than 8000 probes on a 4GB laptop
bwnet = blockwiseModules(datExprs_A, maxBlockSize = 6074,
                         power = 28, networkType= "signed", TOMType = "signed", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE,
                         corType = "bicor",
                         saveTOMs = TRUE,
                         saveTOMFileBase = "TOM-blockwise",
                         verbose = 3)
groupALabels <- bwnet$colors
numModulesA <- max(groupALabels)
groupAColours = labels2colors(groupALabels)

setLabels = c("A", "B")
multiExpr = list(A = list(data = datExprs_A), Male = list(data = datExprs_B))
multiColour = list(A = groupAColours)

system.time( {
  mp = modulePreservation(multiExpr, multiColour,
                          referenceNetworks = 1,
                          nPermutations = 200,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3)
} )
# Save the results
save(mp, file = "modulePreservation.RData")

ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1])
# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1]

# leave grey and gold modules out (jk)
plotMods = !(logical(length(modColors)))
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
#pdf(fi="Plots/BxHLiverFemaleOnly-modulePreservation-Zsummary-medianRank.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
  # For Zsummary, add threshold lines
  if (p==2) {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  } }
# If plotting into a file, close it
#dev.off();

# let's see what the least preserved module is
midnightBlueModuleGenes <- (groupAColours=="midnightblue")
genes <- colnames(datExprs_A)
midnightBlueModuleIDs <- genes[midnightBlueModuleGenes]
write(midnightBlueModuleIDs,file=paste("/Users/maryhoekstra/Desktop/midnightBlueModule_A.txt"))
