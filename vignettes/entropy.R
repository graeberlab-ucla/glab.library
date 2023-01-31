# install.packages("DescTools")
require(DescTools) #Entropy
# Entropy(x, y = NULL, base = 2, …)
# MutInf(x, y, base = 2, …)

# 2nd option
# #install.packages("entropy")
# require(entropy)
# # freqs.empirical(y)
# # entropy.empirical(y, unit=c("log", "log2", "log10"))
# # KL.empirical(y1, y2, unit=c("log", "log2", "log10"))
# # chi2.empirical(y1, y2, unit=c("log", "log2", "log10"))
# # mi.empirical(y2d, unit=c("log", "log2", "log10"))
# # chi2indep.empirical(y2d, unit=c("log", "log2", "log10"))

# examples

entropy.sample <- as.data.frame(sapply(rnaseq.log.pcode[,-1], function(x) Entropy((x))))



