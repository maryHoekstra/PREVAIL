# measure differential expression between groups A and B on a specific day

# create expression set with only samples from that day
analyzeByDay <- function(filteredESet,dayString) {
  
  subsetted.eset <- filtered.eset[ , filtered.eset$Day==dayString]
  
  # create vectors indicating whether a sample belongs to group A or B
  groupA <- as.integer(subsetted.eset$code=='A')
  groupB <- as.integer(subsetted.eset$code=='B') 
  
  # create design matrix with the two groups we are interested in 
  design <- cbind(groupA, groupB)
  colnames(design) <- c('groupA', 'groupB')
  AvsB <- "groupA-groupB"
  
  diffs <- runLimma(exprs(subsetted.eset),design,AvsB)
  
  return(diffs)
}

# determine success of randomization by measuring differential expression between groups A and B at baseline
diffs <- analyzeByDay(filtered.eset, "day 1")
