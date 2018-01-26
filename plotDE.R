# PLOT DE
# graph number of differentially expressed genes between baseline and the other time points for each group

# last modified: January 15th, 2018

numDE_days7and14_A <- findNumberDE("day 1",c("day 7","day 14"),"A")
numDE_days7and14_B <- findNumberDE("day 1",c("day 7","day 14"),"B")

# find number of DE genes from day 1 to a given day for all days
findNumberDE <- function(earlyDay,lateDay,code) {
  earlyVsLateESet <- getESet(filtered.eset,earlyDay,lateDay)
  design <- createDesign(earlyVsLateESet,earlyDay,lateDay)
  diffs <- getDiffs(earlyVsLateESet,groupCode=code,design)
  geneList_up <- subsetByRegulation(diffs,"up")
  geneList_down <- subsetByRegulation(diffs,"down")
  numDE <- length(geneList_up) + length(geneList_down)
  
  return(numDE)
}

earlyDay <- "day 1"
code="A"
numDE_A_day3 <- findNumberDE(earlyDay,"day 3",code) #65
numDE_A_day7 <- findNumberDE(earlyDay,"day 7",code) #6
numDE_A_day14 <- findNumberDE(earlyDay,"day 14",code) #43
numDE_A_day21 <- findNumberDE(earlyDay,"day 21",code) #22
numDE_A_day28 <- findNumberDE(earlyDay,"day 28",code) #2

earlyDay <- "day 1"
code="B"
numDE_B_day3 <- findNumberDE(earlyDay,"day 3",code) #342
numDE_B_day7 <- findNumberDE(earlyDay,"day 7",code) #581
numDE_B_day14 <- findNumberDE(earlyDay,"day 14",code) #588
numDE_B_day21 <- findNumberDE(earlyDay,"day 21",code) #170
numDE_B_day28 <- findNumberDE(earlyDay,"day 28",code) #7

numDE_A <- c(0,numDE_A_day3,numDE_A_day7,numDE_A_day14,numDE_A_day21,numDE_A_day28)
numDE_B <- c(0,numDE_B_day3,numDE_B_day7,numDE_B_day14,numDE_B_day21,numDE_B_day28)

# plot number of DE genes                  
plot(x = numDE_B,type = "l",
     main = "Number of differentially expressed genes over time",
     xlab = "Number of days from baseline",ylab = "Number of differentially expressed genes",
     axes = FALSE,
     col = "blue",
     ylim = c(0,600))
axis(2,at=seq(-250,800,100),labels=seq(-250,800,100))
axis(1, pos=0,at=c(1,2,3,4,5,6),labels=c("0","3","7","14","21","28"))

#lines(x = numDE_B_down, col = "blue")
lines(x = numDE_A, col = "red")
#lines(x = numDE_A_down, col = "red")

legend(x = "topright",legend = c("Placebo group", "Lactoferrin group"),lty = c(1,1),col = c("red","blue"))
