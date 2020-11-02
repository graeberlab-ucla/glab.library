### NOTES ####

#LB # (LB..)
#IOIS # (IOIS..)
#Short ID
#Line # (Line..)

#Each patient is given a unique IOIS ID (and an Medical Record Number (MRN))
#Each tumor/case/resection is given an unique LB#
#A patient with multiple tumor resections is given a unique LB# for each surgery
#Using LB# is more appropriate than using IOIS#, as each resection will be associated with unique tumor model
#for patient ids we are transitioning to the "short ID" column to help identify its matched models more easily 
#(e.g. PT025 and GS025). For grouping samples you can still use the LB number or the "Line #" column.


#GE/Gene Expression
#  Fusion
#  Splice

#WES BAM
#  Consensus MT
#  CNA matrix
#
#  Sequenza = purity/ploidy

#will add 
#'match normal'
#'germmline mutation call'
#'what type of somatic mutation called'

#Flags:

#PDX Yes	
#NS Yes

#BLK: Frozen Bulk Tissue, no PT was available for seq, needed tumor DNAseq so proceeded with bulk tissue
#CD45: CD45+ cell isolate (incl. Tumor infiltrating immune cells)
#DGC: Differentiated Gliomasphere Cultures, gliomaspheres grown in presence of FBS
#GS: Gliomasphere, direct from patient, serum-free
#NRM: Normal samples, specifics of sample type (PBMC, Blood, CD45) can be found in Sample ID
#PRE: PRE-Sorted Patient Tumor samples, tumor samples still containing CD45+ cells, similar to TCGA
#PT: Purified (CD45-) Patient Tumors 
#SDX: Sphere-Derived Xenograft (PT>GS>SDX), first grown as GS culture, then implanted orthotopically
#SQDX: Subcu(Q)taneous-Derived (Orthotopic) Xenograft (PT>SQX>SQDX), first grown in flank, then serially implanted orthotopically
#SQX: Subcu(Q)taneous Xenograft (PT>SQX), patient brain directly into mouse flank
#XDS: Xenograft-Derived Sphere (PT>PDOX>GS), first grown orthotopically in PDOX, then cultured as GS
#XG: Xenograft...should relabel (PT>PDOX), patient brain directly into mouse orthotopically



if(!require(dplyr)) install.packages("dplyr", repos = "http://cran.us.r-project.org")
require(plyr)

#smd.c <- read.csv("/Users/tgraeber/Box/Graeber collaboration/GliomaPDOX Project/Bioinformatics/Metadata/tg_copies/Sequencing Metadata.csv",
#                header = TRUE, sep = ",", quote = "\"",
#                dec = ".", fill = TRUE) #, comment.char = "")


smd <- read.delim("/Users/tgraeber/Box/Graeber collaboration/GliomaPDOX Project/Bioinformatics/Metadata/tg_copies/Sequencing Metadata.txt",
           header = TRUE, sep = "\t", quote = "\"",
           dec = ".", fill = TRUE, comment.char = "")

#remove leading and trailing whitespace
smd <- trimws(smd, which = "both", whitespace = "[ \t\r\n]")

#remove samples that are flagged
smd.noflag <-subset(smd, FLAG != "Y")


#focus on PT, XG, GS sample types
smd.noflag.ptxggs <-subset(smd.noflag, Sample.Type == "PT" | Sample.Type == "XG" | Sample.Type == "GS")



#focus on samples with gene expression data
smd.noflag.ptxggs.ge <-subset(smd.noflag.ptxggs, GE == "Y")


#### calculate which patient cohorts have triplets (PT-XG-GS), or only doublets (PT-XG)

smd.noflag.pt.ge = subset(smd.noflag.ptxggs.ge, Sample.Type == "PT")
smd.noflag.xg.ge = subset(smd.noflag.ptxggs.ge, Sample.Type == "XG")
smd.noflag.gs.ge = subset(smd.noflag.ptxggs.ge, Sample.Type == "GS")


if (0) {
# identifying cases of multiple tumnor models created from the same patient
x <- ddply(smd.noflag.ptxggs.ge, .(IOIS..), mutate, count = length(unique(LB..)))
#y <- ddply(smd, .(IOIS..), mutate, count = length(unique(LB..)))
#x$count
x <- subset(x, count > 1)
#move the new variable to the front for quicker manual browsing
x <- x %>%
  select(count, Line.., everything())
}


#count each category of PT, XG, GS pre lab id
smd.noflag.ptxggs.ge <- smd.noflag.ptxggs.ge %>%
  group_by(LB..) %>%
  mutate(nPT = sum(Sample.Type=="PT")) %>%
  ungroup()

smd.noflag.ptxggs.ge <- smd.noflag.ptxggs.ge %>%
  group_by(LB..) %>%
  mutate(nXG = sum(Sample.Type=="XG")) %>%
  ungroup()

smd.noflag.ptxggs.ge <- smd.noflag.ptxggs.ge %>%
  group_by(LB..) %>%
  mutate(nGS = sum(Sample.Type=="GS")) %>%
  ungroup()

#batch control GS sample 
#smd.noflag.ptxggs.ge.3132 <- subset(smd.noflag.ptxggs.ge, LB..=="3132")


smd.noflag.ptxggs.ge$triple = ifelse((!is.na(smd.noflag.ptxggs.ge$nPT) & smd.noflag.ptxggs.ge$nPT > 0) & 
                               (!is.na(smd.noflag.ptxggs.ge$nXG) & smd.noflag.ptxggs.ge$nXG > 0) & 
                               (!is.na(smd.noflag.ptxggs.ge$nXG) & smd.noflag.ptxggs.ge$nGS > 0),1,0)

smd.noflag.ptxggs.ge$double = ifelse((!is.na(smd.noflag.ptxggs.ge$nPT) & smd.noflag.ptxggs.ge$nPT > 0) & 
                               (!is.na(smd.noflag.ptxggs.ge$nXG) & smd.noflag.ptxggs.ge$nXG > 0),1,0)

smd.noflag.ptxggs.ge$only_double = ifelse(smd.noflag.ptxggs.ge$triple < 1 & smd.noflag.ptxggs.ge$double > 0,1,0)

#list the triples (PT-XG-GS)
unique(subset(smd.noflag.ptxggs.ge$LB.., smd.noflag.ptxggs.ge$triple == 1))

#list the doubles (PT-XG)
unique(subset(smd.noflag.ptxggs.ge$LB.., smd.noflag.ptxggs.ge$only_double == 1))

#list the triples (PT-XG-GS)
triples.ge <- unique(subset(smd.noflag.ptxggs.ge$LB.., smd.noflag.ptxggs.ge$triple == 1))
triples.ge

#list the doubles (PT-XG)
doubles.ge <- unique(subset(smd.noflag.ptxggs.ge$LB.., smd.noflag.ptxggs.ge$only_double == 1))
doubles.ge





#focus on samples with WES data
smd.noflag.ptxggs.wes <-subset(smd.noflag.ptxggs, WES.BAM == "Y")

#may need to switch from WES.BAM to Consensus.MT and/or CNA.matrix

#### calculate which patient cohorts have triplets (PT-XG-GS), or only doublets (PT-XG)

smd.noflag.pt.wes = subset(smd.noflag.ptxggs.wes, Sample.Type == "PT")
smd.noflag.xg.wes = subset(smd.noflag.ptxggs.wes, Sample.Type == "XG")
smd.noflag.gs.wes = subset(smd.noflag.ptxggs.wes, Sample.Type == "GS")



#count each category of PT, XG, GS pre lab id
smd.noflag.ptxggs.wes <- smd.noflag.ptxggs.wes %>%
  group_by(LB..) %>%
  mutate(nPT = sum(Sample.Type=="PT")) %>%
  ungroup()

smd.noflag.ptxggs.wes <- smd.noflag.ptxggs.wes %>%
  group_by(LB..) %>%
  mutate(nXG = sum(Sample.Type=="XG")) %>%
  ungroup()

smd.noflag.ptxggs.wes <- smd.noflag.ptxggs.wes %>%
  group_by(LB..) %>%
  mutate(nGS = sum(Sample.Type=="GS")) %>%
  ungroup()

#batch control GS sample 
#smd.noflag.ptxggs.wes.3132 <- subset(smd.noflag.ptxggs.wes, LB..=="3132")


smd.noflag.ptxggs.wes$triple = ifelse((!is.na(smd.noflag.ptxggs.wes$nPT) & smd.noflag.ptxggs.wes$nPT > 0) & 
                                       (!is.na(smd.noflag.ptxggs.wes$nXG) & smd.noflag.ptxggs.wes$nXG > 0) & 
                                       (!is.na(smd.noflag.ptxggs.wes$nXG) & smd.noflag.ptxggs.wes$nGS > 0),1,0)

smd.noflag.ptxggs.wes$double = ifelse((!is.na(smd.noflag.ptxggs.wes$nPT) & smd.noflag.ptxggs.wes$nPT > 0) & 
                                       (!is.na(smd.noflag.ptxggs.wes$nXG) & smd.noflag.ptxggs.wes$nXG > 0),1,0)

smd.noflag.ptxggs.wes$only_double = ifelse(smd.noflag.ptxggs.wes$triple < 1 & smd.noflag.ptxggs.wes$double > 0,1,0)



#list the triples (PT-XG-GS)
triples.wes <- unique(subset(smd.noflag.ptxggs.wes$LB.., smd.noflag.ptxggs.wes$triple == 1))
triples.wes
  
#list the doubles (PT-XG) [i.e., no GS]
doubles.wes <- unique(subset(smd.noflag.ptxggs.wes$LB.., smd.noflag.ptxggs.wes$only_double == 1))
doubles.wes




### compare ge and wes triples and doubles
#the output are "LB #" (LB..)

#list triples in ge that are not triples in wes
triples.ge[!(triples.ge %in% triples.wes)]

#list triples in ge that are not triples nor doubles in wes
triples.ge[!(triples.ge %in% triples.wes) & !(triples.ge %in% doubles.wes)]

##list doubles in ge that are not doubles in wes
#doubles.ge[!(doubles.ge %in% doubles.wes)]

#list doubles in ge that are not triples or doubles in wes
doubles.ge[!(doubles.ge %in% triples.wes) & !(doubles.ge %in% doubles.wes)]


#list triples in wes that are not triples in ge
triples.wes[!(triples.wes %in% triples.ge)]

#list triples in wes that are not triples nor doubles in ge
triples.wes[!(triples.wes %in% triples.ge) & !(triples.wes %in% doubles.ge)]

##list doubles in wes that are not doubles in ge
#doubles.wes[!(doubles.wes %in% doubles.ge)]

#list doubles in wes that are not triples or doubles in ge
doubles.wes[!(doubles.wes %in% triples.ge) & !(doubles.wes %in% doubles.ge)]








#smd.noflag.ptxggs.3195 <- subset(smd.noflag.ptxggs, LB..=="3195")
#smd.noflag.ptxggs.3195$GE
#smd.noflag.ptxggs.3195$WES.BAM


### example output -- the output are "LB #" (LB..)
# ### compare ge and wes triples and doubles
# #the output are "LB #" (LB..)
# 
# #list triples in ge that are not triples in wes
# triples.ge[!(triples.ge %in% triples.wes)]
#[1] 2857 3166 3386
# 
# #list triples in ge that are not triples nor doubles in wes
# triples.ge[!(triples.ge %in% triples.wes) & !(triples.ge %in% doubles.wes)]
#[1] 2857 3166 3386
# 
# ##list doubles in ge that are not doubles in wes
# #doubles.ge[!(doubles.ge %in% doubles.wes)]
# 
# #list doubles in ge that are not triples or doubles in wes
# doubles.ge[!(doubles.ge %in% triples.wes) & !(doubles.ge %in% doubles.wes)]
#[1] 2870 2873 2978 3063 3505 3811
# 
# 
# #list triples in wes that are not triples in ge
# triples.wes[!(triples.wes %in% triples.ge)]
#[1] 3133 3865 3872
# 
# #list triples in wes that are not triples nor doubles in ge
# triples.wes[!(triples.wes %in% triples.ge) & !(triples.wes %in% doubles.ge)]
#[1] 3133 3865 3872
# 
# ##list doubles in wes that are not doubles in ge
# #doubles.wes[!(doubles.wes %in% doubles.ge)]
# 
# #list doubles in wes that are not triples or doubles in ge
# doubles.wes[!(doubles.wes %in% triples.ge) & !(doubles.wes %in% doubles.ge)]
#[1] 3195





#length(unique(smd$LB..))
#length(unique(smd$IOIS..))

#length(unique(smd.noflag.ptxggs$LB..))
#length(unique(smd.noflag.ptxggs$IOIS..))






### not used ###


#for (i in 1:length(unique(smd1.ge.mod$LB..))) {
#  print(unique(smd1.ge.mod$LB..)[i])
#  #smd1.ge.mod$nPT[]
#}

#smd1.ge.mod <- smd1.ge %>%
#  group_by(LB..) %>%
#  mutate(n_cohort = n()) %>%
#  ungroup()

#smd1.ge.mod_xx = subset(smd1.ge.mod, n_cohort == 1)
#smd1.ge.mod_xxx = subset(smd1.ge.mod, n_cohort == 3)


