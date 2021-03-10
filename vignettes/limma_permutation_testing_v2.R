#### LIMMA function ####
# input should be raw counts, that are then normalized by limma's function voom similar to log transformation

run.limma <- function(exp_matrix, labels, contrasts_vector){
  design <- model.matrix(~ 0 + labels)
  colnames(design) <- levels(factor(labels))
  v <- voom(exp_matrix, design, plot = F)
  lm <- lmFit(v, design)
  contrasts <- makeContrasts(contrasts = contrasts_vector, levels = lm$design)
  fit <- contrasts.fit(lm, contrasts)
  EB <- eBayes(fit)
  tt <- topTable(EB, number = length(lm$sigma))
  return(tt)
}

#### DESeq2 function ####
# input should be raw counts, that are then normalized by DESeq  

run.deseq <- function(exp_matrix, labels, contrasts_vector){
  dds <-DESeqDataSetFromMatrix(countData = round(exp_matrix), colData = data.frame(type = labels), design = ~ type)
  seq <- DESeq(dds)
  res <- results(seq, contrast = c("type", contrasts_vector))
  res_filtered <- res[!is.na(res$padj),]
  res_filtered <- res_filtered[res_filtered$padj < 0.05,]
  res_ordered <- res_filtered[order(res_filtered$padj),]
  return(res_ordered)
}

#### Permutations ####

# expression matrix (n genes x m samples)
#CHANGE
exp_matrix <- norm_xg_matrix[,1:10]
# number of permutations
#CHANGE
num_perms <- 20
# label vector for true groups
#CHANGE
group_labels <- c(rep("one", 5), rep("two", 5))

group_size <- sort(table(group_labels))[1]
group_names <- names(table(group_labels))
#
max_permutations <- choose(length(group_labels), group_size)
if(max_permutations < num_perms){
  print("Invalid number of permutations requested")
}
# initialize matrix to hold p-value results for each gene in each permutation(n genes x m permutations)
perm_stats <- matrix(nrow = nrow(exp_matrix), ncol = num_perms)
# intialize matrix for holding permuted labels
perm_log <- matrix(nrow = ncol(exp_matrix), ncol = num_perms+1)
perm_log[,1] <- group_labels
# permute the labels and run differential expression
pb <- txtProgressBar(min = 0, max = num_perms, initial = 0, style = 3)
for(i in 1:num_perms){
  temp_inds <- sample(1:length(group_labels), length(group_labels), replace = F)
  temp_labels <- group_labels[temp_inds]
  while(max(colSums(temp_labels == perm_log), na.rm = T) == nrow(perm_log)){
    temp_inds <- sample(1:length(group_labels), length(group_labels), replace = F)
    temp_labels <- group_labels[temp_inds]
  }
  perm_log[,i+1] <- temp_labels
  t_res <- run.limma(exp_matrix, temp_labels, paste0(group_names[1], "-", group_names[2]))
  perm_stats[,i] <- t_res$adj.P.Val
  setTxtProgressBar(pb, i)
}

# resetting p-values of NA to 1
perm_stats[is.na(perm_stats)] <- 1
# removing permutations where no gene is statistically significant
no_sig_perms <- -which(colSums(perm_stats) == nrow(perm_stats))
if(length(no_sig_perms)){
  print(paste0(length(no_sig_perms), " permutations with zero significantly different genes removed from analysis"))
  perm_stats <- perm_stats[,-no_sig_perms]
}

#### Actual comparison ####

t_res <- run.limma(exp_matrix, group_labels, paste0(group_names[1], "-", group_names[2]))
p_values <- t_res$adj.P.Val

#### Counting ####

# set range of p-values to evaluate
# CHANGE
cutoffs <- seq(0.75, 0.05, by = -0.05)

# initialize matrix to hold count information for each permutation
perm_counts <- matrix(nrow = length(cutoffs), ncol = ncol(perm_stats))
# count number of genes below p-value threshold in each permutation
for(i in 1:ncol(perm_stats)){
  count_stat <- sapply(1:length(cutoffs), function(j){
    sum(perm_stats[,i] < cutoffs[j])
  })
  perm_counts[,i] <- count_stat
}
# count number of genes below p-value threshold in "real" stratification
true_count <- sapply(1:length(cutoffs), function(i){
  sum(p_values < cutoffs[i])
})

# set percentile lines to appear on graph
# CHANGE
quantiles <- c(0.95, 0.9, 0.5, 0.99)

# evaluate number of genes below p-value at specific percentiles for all permutations
quantile_mat <- sapply(1:length(quantiles), function(a){
  count <- sapply(1:length(cutoffs), function(i){
    quantile(perm_counts[i,], quantiles[a])
  })
})

#### Plotting ####

plot(quantile_mat[,1], type = "l", lty = 3, col = "red", lwd = 2, xaxt= "n", 
     xlab = "Adjusted p-value", ylab = "# of genes below threshold")
axis(1, at = seq(1, 15, by = 2), cutoffs[seq(1, 15, by = 2)])
x <- seq(1:15)
y1 <- quantile_mat[,3]
y2 <- quantile_mat[,2]
polygon(c(x,rev(x)),c(y2,rev(y1)),col=rgb(1,0,0,0.35), border = NA)

y1 <- quantile_mat[,1]
polygon(c(x,rev(x)),c(y2,rev(y1)),col=rgb(1,0,0,0.2), border = NA)

y2 <- quantile_mat[,4]
polygon(c(x,rev(x)),c(y2,rev(y1)),col=rgb(1,0,0,0.1), border = NA)

lines(quantile_mat[,2], lty = 4, col = "red", lwd = 2)
lines(quantile_mat[,3], col = "red")
points(quantile_mat[,3], pch = 15, col = "red")
points(true_count, pch = 16, col = "black")
lines(true_count, col = "black")
lines(quantile_mat[,4], lty= 2, col = "red", lwd = 2)
legend("topright", legend = c("Met Status", "1st percentile",
                              "5th percentile", "10th percentile", "median"),
       lty = c(1, 2,3,4,1), col = c("black", rep("red", 4)), pch = c(16, rep(NA, 3), 15))
