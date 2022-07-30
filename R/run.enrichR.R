#' run.enrichR
#' 
#' runs a full enrichR analhysis
#' 
#' @param topN list of the genes for the enrichmnet analysis
#' @return enricher.df
#' 
#' @importFrom enrichR enrichr
#' 
#' @export
#' 

run.enrichR <- function(topN){ #topN = top500
  #library(enrichR)
  # listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes
  websiteLive <- TRUE
  dbs <- listEnrichrDbs()
  if (is.null(dbs)) websiteLive <- FALSE
  if (0 && websiteLive) head(dbs)
  
  #
  dbs.run <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015")
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
  
  enricher.df = bind_rows(enriched, .id = "Database")
  # if (0) {write.table(enricher.df, "enricher_proeteomics_extreme_top500.txt", sep = "\t", col.names = T, row.names = F, quote = F)}
  # if (0) {write.table(enricher.df, "enricher M249 VSR FA CRISPR extreme.rank top500.txt", sep = "\t", col.names = T, row.names = F, quote = F)}
  # if (0) {write.table(enricher.df, "enricher M249 VSR FA CRISPR and proteomics extreme.rank top500.txt", sep = "\t", col.names = T, row.names = F, quote = F)}
  # if (0) {write.table(enricher.df, "enricher M249 VSR FA geneexp extreme.rank top500.txt", sep = "\t", col.names = T, row.names = F, quote = F)}
  # if (0) {write.table(enricher.df, "enricher M249 VSR FA CRISPR and geneexp extreme.rank top500.txt", sep = "\t", col.names = T, row.names = F, quote = F)}
  # if (0) {write.table(enricher.df, "enricher M249 VSR FA CRISPR geneexp proteomics extreme.rank top500.txt", sep = "\t", col.names = T, row.names = F, quote = F)}
  
  enricher.df$logP = -log10(as.numeric(enricher.df$P.value))
  enricher.df$logP[is.na(enricher.df$logP)] = -1
  
  enricher.df$notes = ""
  
  enricher.df <- enricher.df %>% select(Database,Term,notes,logP,everything())
  enricher.df <- enricher.df %>% arrange(-logP)
  
  return(enricher.df)
  
}
