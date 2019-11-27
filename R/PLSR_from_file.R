#' Partial Least Squares Regression from file
#' Last mod: 6/17/18, ksheu
#'
#' Writes out Xscores, Xloadings, and PVE.
#' Plot scores output with plot_pls()
#'
#' @param file file for X matrix
#' @param sample.names vector of names of the samples
#' @param sample.type vector of sample groups
#' @param y.response numeric vector of response values in same order as samples
#' @param comps number of components to compute
#' @param scale default=F
#' @param labels label the plot, default = F
#' @param comp.x,comp.y comps to display
#' @param title title of the plot
#' @param fread defaylt = F
#'
#' @export
#'

PLSR_from_file = function(file, sample.names, sample.type, y.response, title = "PLSR",comps = 5, scale = F, comp.x = "comp.1", comp.y = "comp.2", labels = F, fread = F){

  requireNamespace(mixOmics)
  requireNamespace(ggplot2)

  if (fread == T) {
    data = fread(file)
    data = data[rowSums((data[, -1, with = F] == 0)) < ncol(data[-1]),
                ]
    t.data = t(data[, -1, with = F])
    colnames(t.data) = data$gene
    y.response = (data.frame(y.response)[match(rownames(t.data), sample.names), ])
    y.response = as.matrix(y.response)

    pls.fit = pls(X = t.data, Y = y.response, scale = scale, ncomp = comps)
    print(pls.fit$explained_variance$X)

    # #do CV
    # set.seed(1) # for reproducibility here, only when the `cpus' argument is not used
    # perf.pls = perf(pls.fit, validation = "loo", #folds = 5,
    #                    progressBar = TRUE, auc = TRUE, nrepeat = 1)
    # print(perf.pls$error.rate)  # error rates
    # plot(perf.pls, col = color.mixo(1:3), sd = TRUE, legend.position = "horizontal")

    #write out
    x.variates = data.frame(pls.fit$variates$X)
    x.loadings = data.frame(pls.fit$loadings$X)
    x.exp_variance = data.frame(pls.fit$explained_variance$X)
    variates.X = cbind(Score = rownames(pls.fit$variates$X), x.variates)
    loadings.X = cbind(Loading = rownames(pls.fit$loadings$X), x.loadings)
    rownames(x.exp_variance) = paste0("comp.",seq(1,nrow(x.exp_variance)))

  }
  else {
  data = read.table(file, sep='\t',header=T,stringsAsFactors=FALSE, quote = "")
  data = data[rowSums((data[, -1] == 0)) < ncol(data[-1]), ] #remove genes with no variance
  rownames(data) = make.names(data[, 1], unique=TRUE)
  t.data = data.frame(t(data[, -1]))
  y.response = (data.frame(y.response)[match(rownames(t.data), sample.names), ])
  y.response = as.matrix(y.response)

  pls.fit = pls(X = t.data, Y = y.response, scale = scale, ncomp = comps)
  print(pls.fit$explained_variance$X)

  # #do CV
  # set.seed(1) # for reproducibility here, only when the `cpus' argument is not used
  # perf.pls = perf(pls.fit, validation = "Mfold", folds = 5,
  #                 progressBar = TRUE, auc = TRUE, nrepeat = 10)
  # print(perf.pls$error.rate)  # error rates
  # plot(perf.pls, col = color.mixo(1:3), sd = TRUE, legend.position = "horizontal")
  #
  #write out
  x.variates = data.frame(pls.fit$variates$X)
  x.loadings = data.frame(pls.fit$loadings$X)
  x.exp_variance = data.frame(pls.fit$explained_variance$X)
  variates.X = cbind(Score = rownames(pls.fit$variates$X), x.variates)
  loadings.X = cbind(Loading = rownames(pls.fit$loadings$X), x.loadings)
  rownames(x.exp_variance) = paste0("comp.",seq(1,nrow(x.exp_variance)))
  }
  write.table(as.data.frame(variates.X), paste0(gsub(".txt", "", file), "_PLSR_Xscores.txt"), sep = "\t", row.names = F, quote = F)
  write.table(as.data.frame(loadings.X), paste0(gsub(".txt", "", file), "_PLSR_Xloadings.txt"), sep = "\t", row.names = F, quote = F)
  write.table(as.data.frame(x.exp_variance), paste0(gsub(".txt", "", file), "_PLSR_Xpve.txt"), sep = "\t", row.names = T, quote = F)

  #plot it
  # variates.X$type = sample.type[match(rownames(variates.X), sample.names)]
  # pcx.y <- ggplot(variates.X, aes_string(x = comp.x, y = comp.y)) + geom_point(size = I(2), aes(color = factor(type))) +
  #   theme(legend.position = "right", plot.title = element_text(size = 30), legend.text = element_text(size = 22),
  #   legend.title = element_text(size = 20), axis.title = element_text(size = 30), legend.background = element_rect(),
  #   axis.text.x = element_text(margin = margin(b = -2)), axis.text.y = element_text(margin = margin(l = -14))) +
  #   labs(title = title) + theme_bw() +
  #   if (labels == TRUE) {
  #     geom_text(data = variates.X, mapping = aes(label =(rownames(variates.X))), check_overlap = TRUE, size = 3)
  #   }
  # pcx.y
}
