# module preservation analysis
# see if modules in a reference set are preserved in a test set

findPreservation <- function(refExprs,testExprs,refPower) {
  # call getModules on reference dataset
  bwnet <- getModules(refExprs,refPower)
  refLabels <- bwnet$colors
  refColours <- labels2colors(refLabels)
  
  setLabels = c("Reference","Test")
  multiExpr = list(Reference = list(data = refExprs), Test = list(data = testExprs))
  multiColour = list(Reference = refColours)
  
  mp = modulePreservation(multiExpr, multiColour,
                               referenceNetworks = 1,
                               nPermutations = 200,
                               randomSeed = 1,
                               quickCor = 0,
                               verbose = 3)
  ref = 1
  test = 2
  statsObs = cbind(mp_B$quality$observed[[ref]][[test]][, -1], mp_B$preservation$observed[[ref]][[test]][, -1])
  statsZ = cbind(mp_B$quality$Z[[ref]][[test]][, -1], mp_B$preservation$Z[[ref]][[test]][, -1])
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
  
  return (mp)
}

getModuleGenes <- function(geneList,mpList,moduleColour) {
  colours <- labels2colors(mpList$colors)
  moduleGenes <- (colours==moduleColour)
  moduleIDs <- geneList[moduleGenes]
  return (moduleIDs)
}

writeFile <- function(IDs,fileName) {
  write(IDs,file=paste("/Users/maryhoekstra/Desktop/",fileName,".txt",sep=""))  
}

mp_AvsB <- findPreservation(refExprs = datExprs_A,testExprs = datExprs_B,refPower = 12)
# save the results
save(mp_AvsB, file = "modulePreservation_AvsB.RData")
moduleIDs <- getModuleGenes(geneList=colnames(datExprs_A),mpList = mp_AvsB,moduleColour = "midnightBlue")
writeFile(moduleIDs,"midnightB_AvsB")

