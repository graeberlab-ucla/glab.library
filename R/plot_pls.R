#' Plot PLS
#'
#' Plots PLS from scores file (output of PLSR_from_file)
#'
#' @param file File containing scores matrix
#' @param info.name Vector of sample names
#' @param info.type Vector of sample types in the same order
#' @param title Title of the plot
#' @param labels default=T
#' @param PCx,PCy PCs to display
#' @param ellipse Construct confidence region based on groups in info.type, default = T
#' @param conf default = 0.95
#' @param fliph default = F
#' @param flipv default = F
#'
#' @importFrom ggplot2 ggplot aes aes_string element_rect element_text geom_point geom_text labs margin theme theme_bw
#'
#' @export
#'
plot_pls <- function(file, info.name, info.type, title = "", labels = TRUE, PCx = "comp1", PCy = "comp2", ellipse = F, conf = 0.95,
                     fliph = F, flipv = F) {
  # Input: PLSR scores file to be plotted
  ## process PLS output and adds groupings

  table <- read.table(file, header = TRUE)
  table$type <- info.type[match(table$Score, info.name)]

  if (fliph == T) {
    table[, PCx] <- table[, PCx] * -1
  }
  if (flipv == T) {
    table[, PCy] <- table[, PCy] * -1
  }


  exp_var <- read.delim(paste0(gsub("scores.txt", "", file), "pve.txt"), row.names = 1)
  exp_var$pve <- unlist(round(exp_var[, 1] * 100, digits = 2))
  rownames(exp_var) <- paste0("comp", seq(1, nrow(exp_var)))

  pcx.y <- ggplot(table, aes_string(x = PCx, y = PCy)) +
    geom_point(size = I(3), aes(color = factor(type))) +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 30),
      legend.text = element_text(size = 22),
      legend.title = element_text(size = 20),
      axis.title = element_text(size = 30),
      legend.background = element_rect(),
      axis.text.x = element_text(margin = margin(b = -2)),
      axis.text.y = element_text(margin = margin(l = -14))
    ) +
    guides(color = guide_legend(title = "Type")) +
    labs(
      title = title,
      x = paste0(PCx, " (", exp_var$pve[match(PCx, rownames(exp_var))], "%)"),
      y = paste0(PCy, " (", exp_var$pve[match(PCy, rownames(exp_var))], "%)")
    ) +
    theme_bw(base_size = 18) +
    if (labels == TRUE) {
      ggrepel::geom_text_repel(data = table, mapping = aes(label = Score), size = 3)
    }


  if (ellipse == TRUE) {
    plot(table[, c(PCx, PCy)], main = title)
    ord <- vegan::ordiellipse(table[, c(PCx, PCy)], table$type, kind = "sd", conf = conf)

    cov_ellipse <- function(cov, center = c(0, 0), scale = 1, npoints = 100) {
      theta <- (0:npoints) * 2 * pi / npoints
      Circle <- cbind(cos(theta), sin(theta))
      t(center + scale * t(Circle %*% chol(cov)))
    }

    df_ell <- data.frame(matrix(ncol = 0, nrow = 0))
    for (g in (table$type)) {
      df_ell <- rbind(df_ell, cbind(as.data.frame(with(
        table[table$type == g, ],
        cov_ellipse(ord[[g]]$cov, ord[[g]]$center, ord[[g]]$scale)
      )),
      type = g
      ))
    }

    pcx.y2 <- pcx.y +
      geom_path(data = df_ell, aes(x = df_ell[, PCx], y = df_ell[, PCy], colour = type), size = 1, linetype = 1)
    pcx.y2
  } else {
    pcx.y
  }
}
