# perform a two-sampled t-test
nA = ncol(groupAExprs)
nB = ncol(groupBExprs)
meansA <- rowMeans(groupAExprs)
meansB <- rowMeans(groupBExprs)
stdevA <- rowSds(groupAExprs)
stdevB <- rowSds(groupBExprs)

ttest = (meansA-meansB)/sqrt(stdevA^2/nA + stdevB^2/nB)
DOF = (nA+nB)-2
pvals = 2*(1-pt(abs(ttest),DOF))

# create factor using group codes 
groupFactor <- factor(treatment_group)

# perform a t-test using the rowttests function
rtt <- rowttests(baselineExprs,groupFactor,tstatOnly = FALSE)