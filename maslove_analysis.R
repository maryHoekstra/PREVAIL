splitGeneList <-function(x) {
  return(strsplit(as.character(x), '///')[[1]])
}

# LOAD 'prevail.eset'
clean.eset <- varFilter(prevail.eset, var.cutoff = 0.5)
# early = Day 1 late = Day 14,21,28
clean.eset <- clean.eset[ , clean.eset$code=='A']

EarlyLate.eset <- clean.eset[ , clean.eset$Day == "day 1" | clean.eset$Day == 'day 14' | clean.eset$Day == 'day 21' | clean.eset$Day == 'day 28']
early <- as.integer(EarlyLate.eset$Day == 'day 1')
late <- as.integer(EarlyLate.eset$Day == 'day 14' | EarlyLate.eset$Day == 'day 21' | EarlyLate.eset$Day == 'day 28')
pheno <- pData(EarlyLate.eset)
edata <- exprs(EarlyLate.eset)
batch = EarlyLate.eset$Batch
modcombat = model.matrix(~1, data=pheno)
combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

design <- cbind(early, late)
colnames(design) <- c('early', 'late')
fit <- lmFit(combat_edata, design)
cont.matrix <- makeContrasts(EvsL=early-late, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
diffs <- topTable(fit2, p.value=0.05, adjust='holm', number=5000)
diffs$ID <- as.factor(rownames(diffs))
diffs.anno <- merge(diffs, select(anno, ID, Entrez.Gene))

# get annotations and analyze further
geneEntrez <- select(diffs.anno, Entrez.Gene)
genelist <- unlist(sapply(geneEntrez[,1], splitGeneList))
genelist2 <- genelist[genelist != "---"]
genelist3 <- gsub(' ', '', genelist2)
