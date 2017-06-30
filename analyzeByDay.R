# measure differential expression between groups A and B on a specific day

# create expression set with only samples from that day
analyzeByDay <- function(filteredESet,dayString) {
  
  subsetted.eset <- filteredESet[ , filteredESet$Day==dayString]
  
  # create vectors indicating whether a sample belongs to group A or B
  groupA <- as.integer(subsetted.eset$code=='A')
  groupB <- as.integer(subsetted.eset$code=='B') 
  
  # create design matrix with the two groups we are interested in 
  design <- cbind(groupA, groupB)
  colnames(design) <- c('groupA', 'groupB')
  AvsB <- "groupA-groupB"
  
  # fit a linear model using the design
  fit <- lmFit(subsetted.eset, design)
  
  # create a contrast matrix 
  cont.matrix <- makeContrasts(contrasts=AvsB, levels=design)
  
  # fit a linear model with specified contrasts
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  # adjust p-values based on Holm's method, which controls family-wide error rate
  diffs <- topTable(fit2, p.value=0.05, adjust='holm', number=5000)
  
  return(diffs)
}

# determine success of randomization by measuring differential expression between groups A and B at baseline
diffs <- analyzeByDay(filtered.eset, "day 28")
