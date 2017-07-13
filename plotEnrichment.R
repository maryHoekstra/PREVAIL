# visualize annnotation results

library(gdata)

GO_BP_module14 = read.xls ("/Users/maryhoekstra/Desktop/GO_BP_module14.xlsx", sheet = 1, header = TRUE)
allCounts <- sum(GO_BP_module14$Count)
logPVals <- -log(GO_BP_module14$Benjamini)
meanLogPVal <- mean(logPVals)
maxLogPVal <- max(logPVals)
yMax <- maxLogPVal * 1.1
topTerms <- GO_BP_module14$Term[1:7]
strippedTerms <- unlist(lapply(topTerms,function(string) sub(".*~", "", string)))
barplot(height=logPVals[1:7], 
        main = "Highly Enriched GO Terms",
        xlab = "GO Term", ylab = "Enrichment Score (-log(pval))",
        names.arg = strippedTerms,
        cex.names = 0.7,
        col = "darkblue",
        ylim = c(0,yMax))
abline(h=meanLogPVal,col="red")
legend(x="topright",legend="Mean Enrichment",lwd=1,col="red")
