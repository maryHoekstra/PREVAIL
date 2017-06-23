# module preservation analysis
# see if modules in a reference set are preserved in a test set

allGenes <- colnames(datExprs)

# call getModules on reference dataset
bwnet <- getModules(refExprs,refPower)
refLabels <- bwnet$colors
refColours <- labels2colors(refLabels)

findPreservation <- function(refExprs,refPower,refColours,testExprs) {
  
  setLabels = c("Reference","Test")
  multiExpr = list(Reference = list(data = refExprs), Test = list(data = testExprs))
  multiColour = list(Reference = refColours)
  
  mp = modulePreservation(multiExpr, multiColour,
                               referenceNetworks = 1,
                               nPermutations = 200,
                               randomSeed = 1,
                               quickCor = 0,
                               verbose = 3)
  
  return (mp)
}

mp_A <- findPreservation(refExprs = earlyExprs_A,testExprs = lateExprs_A,refPower = 26)
mp_B <- findPreservation(refExprs = earlyExprs_,testExprs = lateExprs_B,refPower = 26)
mp_BvsA_late <- findPreservation(refExprs = lateExprs_B,testExprs = lateExprs_A,refPower = 16)
# save the results
save(mp_earlyVsLate_A, file = "modulePreservation_EvsL_A.RData")
save(mp_earlyVsLate_B, file = "modulePreservation_EvsL_B.RData")
save(mp_BvsA_late, file = "modulePreservation_BvsA_late.RData")

# plot statistics
ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1])
# compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )
# module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1]

# leave grey and gold modules out
plotMods = !(modColors %in% c("grey","gold"))
# text labels for points
text = modColors[plotMods]
# auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# start the plot
sizeGrWindow(10, 5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE)
  max = max(plotData[, p], na.rm = TRUE)
  # adjust ploting ranges appropriately
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
  # for Zsummary, add threshold lines
  if (p==2) {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  } }

getModuleGenes <- function(geneList,geneColours,moduleColour) {
  moduleGenes <- (geneColours==moduleColour)
  moduleIDs <- geneList[moduleGenes]
  return (moduleIDs)
}

writeFile <- function(IDs,fileName) {
  write(IDs,file=paste("/Users/maryhoekstra/Desktop/",fileName,".txt",sep=""))  
}

nonpreservedModules <- c("lightcyan")
for (i in 1:length(nonpreservedModules)) {
  moduleIDs <- getModuleGenes(geneList=allGenes,geneColours = refColours,moduleColour = nonpreservedModules[i])
  writeFile(moduleIDs,paste(nonpreservedModules[i],"_BvsA_late",sep = ""))
}

