datExprs <- t(collapsedExprs)
samples <- rownames(datExprs)
early <- c("day 1")
late <- c("day 14")
earlyLate.eset <- getESet(filtered.eset,early,late)
groupA.eset = earlyLate.eset[,earlyLate.eset$code=="A"]
groupB.eset = earlyLate.eset[,earlyLate.eset$code=="B"]
timepointLabels <- sub("day 14","Late",groupA.eset$Day)
timepointLabels <- sub("day 1","Early",timepointLabels)

early.eset <- earlyLate.eset[,earlyLate.eset$Day=="day 1"]
late.eset <- earlyLate.eset[,earlyLate.eset$Day=="day 14"]
early_groupB.eset <- early.eset[, early.eset$code=="B"]
late_groupB.eset <- late.eset[, late.eset$code=="B"]
early_groupA.eset <- early.eset[, early.eset$code=="A"]
late_groupA.eset <- late.eset[, late.eset$code=="A"]

earlyIndices_B <- match(colnames(exprs(early_groupB.eset)),samples)
lateIndices_B <- match(colnames(exprs(late_groupB.eset)),samples)
earlyIndices_A <- match(colnames(exprs(early_groupA.eset)),samples)
lateIndices_A <- match(colnames(exprs(late_groupA.eset)),samples)

earlyLateB <- datExprs[c(earlyIndices_B,lateIndices_B),]
earlyLateA <- datExprs[c(earlyIndices_A,lateIndices_A),]

designVec <- sub("Early",1,timepointLabels)
designVec <- as.numeric(sub("Late",0,designVec))
designMat <- cbind(designVec,!designVec)
colnames(designMat) <- c("Early","Late")
patientIDs <- groupA.eset$patient_id
moduleGenes <- getModuleGenes(geneList=allGenes,geneColours = moduleLabelsAutomatic,moduleColour = 14)
which_IDs <- match(moduleGenes,colnames(earlyLateA))
moduleExprs <- earlyLateB[,which_IDs]


corfit <- duplicateCorrelation(t(moduleExprs),designMat,block=patientIDs)
print(corfit$consensus)
fit <- lmFit(t(moduleExprs),designMat,block=patientIDs,correlation=corfit$consensus)
cont.matrix <- makeContrasts(contrasts = "Early-Late",levels=designMat)
fit2 <- contrasts.fit(fit,cont.matrix)
fit2 <- eBayes(fit2)
diffs <- topTable(fit2, adjust='holm',number=10000)

numUpregulated <- nrow(diffs[diffs$logFC<0,])
numDownregulated <- nrow(diffs[diffs$logFC>0,])
percentUpregulated <- numUpregulated/ncol(moduleExprs) * 100

# determine which genes are upregulated or downregulated based on fold changes
earlySamples <- moduleExprs[which(timepointLabels=="Early"),]
lateSamples <- moduleExprs[which(timepointLabels=="Late"),]
earlyMeans <- colMeans(earlySamples)
lateMeans <- colMeans(lateSamples)
quotientMeans <- lateMeans/earlyMeans
logFoldChanges <- log2(quotientMeans)
# when late > early, the log fold change is positive since log2(+) is also positive
numUpregulated <- length(logFoldChanges[logFoldChanges>0])
numDownregulated <- length(logFoldChanges[logFoldChanges<0])
percentUpregulated <- numUpregulated/length(logFoldChanges) * 100

plot(logFoldChanges)

# just look at module eigengenes
earlySamples <- MEsAutomatic[which(timepointLabels=="Early"),]
lateSamples <- MEsAutomatic[which(timepointLabels=="Late"),]
earlyMeans <- colMeans(earlySamples)
lateMeans <- colMeans(lateSamples)
quotientMeans <- lateMeans/earlyMeans
logFoldChanges <- log2(quotientMeans)
# when late > early, the log fold change is positive since log2(+) is also positive
numUpregulated <- length(logFoldChanges[logFoldChanges>0])
numDownregulated <- length(logFoldChanges[logFoldChanges<0])
percentUpregulated <- numUpregulated/length(logFoldChanges) * 100