
R PCA tutorial set

README

PCA and PLSR scripts, with example input and output.

Recommend running these scripts in RStudio, or a similar R GUI (graphical user interface)


These scripts are basically a wrapper for the R packace 'prcomp'.
	https://stat.ethz.ch/R-manual/R-devel/library/stats/html/prcomp.html


-------------------------------
 
#For help information on each parameter either (i) see the header portion of each .R file, or (ii) when installed as packages in R/RStudio, use help("PCA_from_file") or help('plot_pca') in R or RStudio. 

#Open the four *.R files in RStudio. 

#'source' the four files to load the algorithms into memory. In more advanced applications they can be installed as packages. 

#execute the following commands from the console of the Rstudio GUI:

#set working directory.  Replace "/path/to/working/dir" with the actual working directory where the files are located
setwd("/path/to/working/dir")

#read annotation file----
human.info = read.delim("human.info.rsem.expression.txt")
	
# for PCA run----
PCA_from_file("Beltran_2016_rsem_genes_upper_norm_counts_coding_log2.txt", center = T, scale = F)
plot_pca("Beltran_2016_rsem_genes_upper_norm_counts_coding_log2_prcomp_scores.txt",
    human.info$sample, human.info$type, labels = F, ellipse = T, PCx = "PC1", PCy = "PC2")

# for PLSR run----
PLSR_from_file("Beltran_2016_rsem_genes_upper_norm_counts_coding_log2.txt", human.info$sample, human.info$type,
    (ifelse(human.info$type=="NEPC", 1, 0)), comps = 5, scale = F)
plot_pls("Beltran_2016_rsem_genes_upper_norm_counts_coding_log2_PLSR_Xscores.txt",
    human.info$sample, human.info$type, title = "PLSR", labels =F, PCx = "comp.1", PCy = "comp.2", ellipse = T, conf = 0.90)




#tutorial set initiated by Katerine Sheu