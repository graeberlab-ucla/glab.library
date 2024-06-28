#' run.enrichR
#' 
#' runs a full enrichR analysis
#' 
#' CRAN
#' install.packages("enrichR")  
#' Github
#' library(devtools)
#' install_github("wjawaid/enrichR")
#'   
#' 
#' @param topN list of the genes for the enrichmnet analysis
#' @return enricher.df
#' 
#' @importFrom enrichR enrichr
#' 
#' @export
#' 

# topN = top500


run.enrichR <- function(topN){ # topN = top500 ;  topN = bottom500
  require(enrichR)
  require(dplyr)
  # listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes
  websiteLive <- TRUE
  dbs <- listEnrichrDbs()
  if (is.null(dbs)) websiteLive <- FALSE
  if (0 && websiteLive) head(dbs)
  
  topN = as.vector(topN)
  
  if (0) { #likely not needed, was added during a debugging stage
    # if (sum(grepl("^C.*orf", topN))) { }
    file_background = "/Users/tgraeber/Library/CloudStorage/Dropbox/glab/data/enrichr human background geneset.txt"
    background = read.delim(file_background, header = F , sep = "\t", stringsAsFactors = F)
    
    # topN_missing = topN[!(topN %in% background$V1)]
    topN = topN[topN %in% background$V1]
  }

  
  #dbs.run <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015")
  #dbs.run <- c("ChEA 2022")
  dbs.run <- dbs$libraryName
  
  if (websiteLive) {
    #enriched <- enrichr(c("Runx1", "Gfi1", "Gfi1b", "Spi1", "Gata1", "Kdr"), dbs.run)
    enriched <- enrichr(topN, dbs.run)
  }
  
  if (0) {
    if (websiteLive) enriched[["GO_Biological_Process_2015"]]
    if (websiteLive) plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
  }
  
  # df1 <- enriched[["GO_Biological_Process_2015"]]
  # df2 <- enriched[["GO_Cellular_Component_2015"]]
  
  #remove empty results
  enriched2 <- enriched[sapply(enriched, function(x) dim(x)[1]) > 0]
  
  enricher.df = dplyr::bind_rows(enriched2, .id = "Database")
  
  # debugging
  # for (i in 500:500) { #i=2
  #   print(i)
  #   #topN = top500[1:i]
  #   topN = bottom500[1:i]
  #   enriched <- enrichr(topN, dbs.run)
  #   enricher.df = dplyr::bind_rows(enriched, .id = "Database")
  #   x = dim(enricher.df)
  #   print(x)
  # }
  
  # if (0) {write.table(enricher.df, "enricher_proeteomics_extreme_top500.txt", sep = "\t", col.names = T, row.names = F, quote = F)}
  # if (0) {write.table(enricher.df, "enricher M249 VSR FA CRISPR extreme.rank top500.txt", sep = "\t", col.names = T, row.names = F, quote = F)}
  # if (0) {write.table(enricher.df, "enricher M249 VSR FA CRISPR and proteomics extreme.rank top500.txt", sep = "\t", col.names = T, row.names = F, quote = F)}
  # if (0) {write.table(enricher.df, "enricher M249 VSR FA geneexp extreme.rank top500.txt", sep = "\t", col.names = T, row.names = F, quote = F)}
  # if (0) {write.table(enricher.df, "enricher M249 VSR FA CRISPR and geneexp extreme.rank top500.txt", sep = "\t", col.names = T, row.names = F, quote = F)}
  # if (0) {write.table(enricher.df, "enricher M249 VSR FA CRISPR geneexp proteomics extreme.rank top500.txt", sep = "\t", col.names = T, row.names = F, quote = F)}
  
  enricher.df$logP = -log10(as.numeric(enricher.df$P.value))
  enricher.df$logP[is.na(enricher.df$logP)] = -1
  
  enricher.df$notes = ""
  
  # names(enricher.df)
  
  
  file.erichr.dbs.trim = "/Users/tomgraeber/Dropbox/glab/data/enrichr/enrichr datasets culling Oct 27 2023.txt"
  erichr.dbs.trim = read.delim(file.erichr.dbs.trim, sep = "\t", stringsAsFactors = F)
  erichr.dbs.trim.cut = erichr.dbs.trim[erichr.dbs.trim$trim == 0,] # trim and trim2 are different stringencies of culling
  erichr.dbs.trim2.cut = erichr.dbs.trim[erichr.dbs.trim$trim2 == 0,] # trim and trim2 are different stringencies of culling
  
  enricher.df$trim = ifelse(enricher.df$Database %in% erichr.dbs.trim.cut$Database,1,0)
  enricher.df$trim2 = ifelse(enricher.df$Database %in% erichr.dbs.trim2.cut$Database,1,0)
  
  enricher.df <- enricher.df %>% dplyr::select(Database,Term,notes,starts_with("trim"),logP,everything())
  enricher.df <- enricher.df %>% arrange(-logP)

    # enricher.df.trim = enricher.df
  # 
  # if (1) {
  #   file.erichr.dbs.trim = "/Users/tomgraeber/Dropbox/glab/data/enrichr/enrichr datasets culling Oct 27 2023.txt"
  #   erichr.dbs.trim = read.delim(file.erichr.dbs.trim, sep = "\t", stringsAsFactors = F)
  #   erichr.dbs.trim.cut = erichr.dbs.trim[erichr.dbs.trim$X == 0,]
  #   
  #   enricher.df.trim = enricher.df.trim[!(enricher.df.trim$Database %in% erichr.dbs.trim.cut$Database),]
  #   #enricher.df.trim = enricher.df.trim[!grepl("NIH_Funded_PIs_", enricher.df.trim$Database), ]
  #   
  # }
  
  return(enricher.df)
  
}
