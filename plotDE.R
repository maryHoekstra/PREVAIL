# graph number of differentially expressed genes between baseline and the other time points for each group
             

earlyDays <- c("day 1")
lateDays <- c("day 28")
numDE_B_up <- numeric(6)
numDE_B_down <- numeric(6)
earlyVsLateESet <- getESet(filtered.eset,earlyDays,lateDays)
design <- createDesign(earlyVsLateESet)
diffsB_day28 <- getDiffs(earlyVsLateESet,groupCode="B",design)
geneListB_up <- subsetByRegulation(diffs,"up")
geneListB_down <- subsetByRegulation(diffs,"down")
numDE_B_up[6] <- length(geneListB_up)
numDE_B_down[6] <- length(geneListB_down)  

numDE_A <- c(0,40,480,44,20,2)
numDE_B <- c(0,193,338,605,184,7)
numDE_A_up <- c(0,30,321,38,20,0)
numDE_A_down <- c(0,-10,-191,-7,0,-2)
numDE_B_up <- c(0,153,307,547,163,7)
numDE_B_down <- c(0,-45,-47,-78,-28,0)
                  
plot(x = numDE_B_up,type = "l",
     main = "Number of differentially expressed genes over time",
     xlab = "Day",ylab = "Number of differentially expressed genes",
     axes = FALSE,
     col = "blue",
     ylim = c(-250,600))
axis(2,at=seq(-250,600,100),labels=seq(-250,600,100))
axis(1, pos=0,at=c(1,2,3,4,5,6),labels=c("1","3","7","14","21","28"))

lines(x = numDE_B_down, col = "blue")
lines(x = numDE_A_up, col = "red")
lines(x = numDE_A_down, col = "red")
legend(x = "topright",legend = c("Group A", "Group B"),lty = c(1,1),col = c("red","blue"))
