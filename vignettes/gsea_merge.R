#GSEA merge
#example script/temolate (not a function)


require(ggplot2)
require(data.table)

#### read in gsea files ####
BAP1all.gsea1 <- read.delim("/Users/tgraeber/Dropbox/glab/collab f/Kim Paraiso/glm/BAP1all my_analysis.GseaPreranked.1610680009928/gsea_report_for_na_pos_1610680009928.tsv")
BAP1all.gsea2 <- read.delim("/Users/tgraeber/Dropbox/glab/collab f/Kim Paraiso/glm/BAP1all my_analysis.GseaPreranked.1610680009928/gsea_report_for_na_neg_1610680009928.tsv")
BAP1all.gsea <- rbind(BAP1all.gsea1, BAP1all.gsea2)
BAP1all.gsea$NES <- as.numeric(BAP1all.gsea$NES)

BAP1pos.gsea1 <- read.delim("/Users/tgraeber/Dropbox/glab/collab f/Kim Paraiso/glm/BAP1positive my_analysis.GseaPreranked.1610680005101/gsea_report_for_na_pos_1610680005101.tsv")
BAP1pos.gsea2 <- read.delim("/Users/tgraeber/Dropbox/glab/collab f/Kim Paraiso/glm/BAP1positive my_analysis.GseaPreranked.1610680005101/gsea_report_for_na_neg_1610680005101.tsv")
BAP1pos.gsea <- rbind(BAP1pos.gsea1, BAP1pos.gsea2)
BAP1pos.gsea$NES <- as.numeric(BAP1pos.gsea$NES)

BAP1mix.gsea1 <- read.delim("/Users/tgraeber/Dropbox/glab/collab f/Kim Paraiso/glm/BAP1mixed my_analysis.GseaPreranked.1610679995815/gsea_report_for_na_pos_1610679995815.tsv")
BAP1mix.gsea2 <- read.delim("/Users/tgraeber/Dropbox/glab/collab f/Kim Paraiso/glm/BAP1mixed my_analysis.GseaPreranked.1610679995815/gsea_report_for_na_neg_1610679995815.tsv")
BAP1mix.gsea <- rbind(BAP1mix.gsea1, BAP1mix.gsea2)
BAP1mix.gsea$NES <- as.numeric(BAP1mix.gsea$NES)

BAP1neg.gsea1 <- read.delim("/Users/tgraeber/Dropbox/glab/collab f/Kim Paraiso/glm/BAP1negative my_analysis.GseaPreranked.1610680147807/gsea_report_for_na_pos_1610680147807.tsv")
BAP1neg.gsea2 <- read.delim("/Users/tgraeber/Dropbox/glab/collab f/Kim Paraiso/glm/BAP1negative my_analysis.GseaPreranked.1610680147807/gsea_report_for_na_neg_1610680147807.tsv")
BAP1neg.gsea <- rbind(BAP1neg.gsea1, BAP1neg.gsea2)
BAP1neg.gsea$NES <- as.numeric(BAP1neg.gsea$NES)

  

if (1) {
#### convert undefined NES scores to max value with correct sign ####
max_nes = max(abs(as.numeric(BAP1all.gsea$NES)), na.rm = TRUE)
BAP1all.gsea$NES = ifelse(is.na(BAP1all.gsea$NES), sign(BAP1all.gsea$ES)*max_nes, BAP1all.gsea$NES)

max_nes = max(abs(as.numeric(BAP1pos.gsea$NES)), na.rm = TRUE)
BAP1pos.gsea$NES = ifelse(is.na(BAP1pos.gsea$NES), sign(BAP1pos.gsea$ES)*max_nes, BAP1pos.gsea$NES)

max_nes = max(abs(as.numeric(BAP1mix.gsea$NES)), na.rm = TRUE)
BAP1mix.gsea$NES = ifelse(is.na(BAP1mix.gsea$NES), sign(BAP1mix.gsea$ES)*max_nes, BAP1mix.gsea$NES)

max_nes = max(abs(as.numeric(BAP1neg.gsea$NES)), na.rm = TRUE)
BAP1neg.gsea$NES = ifelse(is.na(BAP1neg.gsea$NES), sign(BAP1neg.gsea$ES)*max_nes, BAP1neg.gsea$NES)
}

#### end of basic read and modify files section ####

#### merge gsea NES files ####

merge1 <- merge(BAP1all.gsea, BAP1pos.gsea, by="NAME", suffixes = c("all", "pos"), all=F)
merge2 <- merge(BAP1mix.gsea, BAP1neg.gsea, by="NAME", suffixes = c("mix", "neg"), all=F)
merge <- merge(merge1, merge2, by="NAME", all=F)
merge$NESall <- as.numeric(merge$NESall)
merge$NESpos <- as.numeric(merge$NESpos)
merge$NESmix <- as.numeric(merge$NESmix)
merge$NESneg <- as.numeric(merge$NESneg)

ggplot(merge, aes(NESall, NESpos)) + geom_point()
ggplot(merge, aes(NESall, NESmix)) + geom_point()
ggplot(merge, aes(NESpos, NESmix)) + geom_point()

if (0) {
  merge[grepl("IFN|interferon", merge$NAME, ignore.case = T),grepl("NES|NAME",colnames(merge))]
  merge[grepl("TNF", merge$NAME),grepl("NES|NAME",colnames(merge))]
}

mergeNES <- merge[,grepl("NAME|ES|NAME|p.val",colnames(merge))]

write.table(mergeNES, file = "mergeNES.out", sep="\t")


#### look up geneset-based results ####

geneset = "BIOCARTA_IL6_PATHWAY"  #rsl3 sensitive - NESmix & NESpos

BAP1all.gsea.geneset <- read.delim(paste0("/Users/tgraeber/Dropbox/glab/collab f/Kim Paraiso/glm/BAP1all my_analysis.GseaPreranked.1610680009928/",geneset,".tsv"))
BAP1pos.gsea.geneset <- read.delim(paste0("/Users/tgraeber/Dropbox/glab/collab f/Kim Paraiso/glm/BAP1positive my_analysis.GseaPreranked.1610680005101/",geneset,".tsv"))
BAP1mix.gsea.geneset <- read.delim(paste0("/Users/tgraeber/Dropbox/glab/collab f/Kim Paraiso/glm/BAP1mixed my_analysis.GseaPreranked.1610679995815/",geneset,".tsv"))
BAP1neg.gsea.geneset <- read.delim(paste0("/Users/tgraeber/Dropbox/glab/collab f/Kim Paraiso/glm/BAP1negative my_analysis.GseaPreranked.1610680147807/",geneset,".tsv"))





