# devtools::install_github("danielbraas/MetabFUN")
#devtools::install_github("juyeki/MetabR")
library(MetabFUN)
library(MetabR)
library(tidyr)
library(readxl)
library(dplyr)
library(ggrepel)
#Sys.setenv(JAVA_HOME = "C:\\Program Files\\Java\\jdk-11.0.8\\")
library(xlsx)
library(writexl)
####
library(tidyverse)

data_dir <- "K:/Teitell lab/JH-07112020_7-8_Vanq"
output_dir <- "K:/Teitell lab/JH-07112020_7-8_Vanq/test3"
diff_computers <- FALSE
ifelse (grepl("ICS", data_dir), ICS <- TRUE, ICS <- FALSE)
#set tracer_type to "none", "partial", or "full"
tracer_type <- "full"
#set normalize to TRUE or FALSE 
normalize <- FALSE
#set correct_for to "CN" of "N"
correct_for <- "C"

#Abbreviation Data
library(googlesheets4)
#Abbrev_NEW
#Abbrev <- read_sheet("https://docs.google.com/spreadsheets/d/118M3rvfJAOQrOoYEHZVjEFg_k0CMwK9i8wq56u_ZXP4/edit?ts=5eaa1a9e#gid=220628921")
#Abbrev_NEW2
#Abbrev <- read_sheet("https://docs.google.com/spreadsheets/d/1He49QJYE1ld0VUgzNTNht-AuoznydEL9ysjCPRkgk1Y/edit#gid=220628921")
Abbrev <- read_excel("C:/Users/FTsang/Downloads/Abbrev_NEW2 (1).xlsx")
Abbrev <- Abbrev[-c(248,249),]

#Title
setwd(data_dir)
Title <- gsub('.xls[x]?','', list.files(pattern='.xls[x]?'))
Title <- gsub('_sample info sheet|_sample info|sample form|_sampleinfosheet|_Sampleinfosheet', '', Title)
Title <- gsub(' ', '', Title)
Title <- Title[1]

#info 
info <- read_excel(list.files()[grep('.xls[x]?',list.files())][1])
info$Condition <- gsub("^\\s+|\\s+$", "", info$Condition)     
info$Condition <- gsub('/', '-', info$Condition)
info$Condition <- gsub('_', '-', info$Condition)
info$Sample <- gsub('-', '.', info$Sample)
info <- info[!grepl("QC.blank", info$Sample), ]
info$Condition <- factor(info$Condition, levels = unique(info$Condition))   

setwd(data_dir)
data <- read_csv(list.files()[grep('filtered',tolower(list.files()))])
colnames(data)[1] <- "Used_ID"
data <- data %>% group_by(Used_ID) %>% do( data.frame(with(data=., .[order(Iso),] )) )

to_keep <- info$Sample
to_keep <- gsub("-", "\\.", to_keep)
data <- data[,colnames(data) %in% c("Used_ID", "Iso", to_keep)]

#data$Used_ID[ which(data$Used_ID == "5M-thioadenosine")] <- "5'-Methylthioadenosine"
for (i in 1:nrow(Abbrev))
  data$Used_ID[as.character(data$Used_ID) == Abbrev$Abb[i]] <- as.character(Abbrev$Used_ID[i])

#ISTD Plot
### Added code to ensure that Norvaline is not included in the ISTD plot
std <- data[which(data$Iso == "Std" & data$Used_ID != "Norvaline"), ]
### Keep Norvaline in data (?)
data <- data[-which(data$Iso == "Std" & data$Used_ID != "Norvaline"), ]
if(sum(grepl("blank", colnames(std))) > 0 )
  std <- std[,-which(grepl("blank", colnames(std)))]

#250K
std_250K <- std %>% dplyr::select(-Iso) %>%
  gather(Sample, Value, -Used_ID) %>%
  group_by(Used_ID) %>%
  mutate(Av = mean(Value, na.rm=T),
         Std = sd(Value, na.rm=T),
         CV = Std / Av,
         ID = paste0(gsub(' Std','',Used_ID), ' ',round(CV*100,0),'%')) %>%
  ungroup()

for(i in 1:nrow(info))
  std_250K$Sample[which(std_250K$Sample == info$Condition[i])] <- info$Sample[i]

setwd(output_dir)

pdf(paste0("QC-ISTDs with 250K-", Title, ".pdf"), width = 14, height = 10)
ggplot(std_250K, aes(Sample, log(Value,10)))+
  geom_point(size=3)+
  facet_wrap(~ID, scales='free')+
  theme_bw()+
  theme(text = element_text(face='bold',size=14),
        axis.text.x = element_text(angle=90, vjust=.3, hjust=1))+
  labs(x='',y='Relative response of ISTD (A.U., log10)') + scale_x_discrete(limits = info[order(info$Run.Order),]$Sample)
dev.off()

#only samples 
std <- std[,-which(grepl("250K", colnames(std)))]
Std <- std %>% dplyr::select(-Iso) %>%
  gather(Sample, Value, -Used_ID) %>%
  group_by(Used_ID) %>%
  mutate(Av = mean(Value, na.rm=T),
         Std = sd(Value, na.rm=T),
         CV = Std / Av,
         ID = paste0(gsub(' Std','',Used_ID), ' ',round(CV*100,0),'%')) %>%
  ungroup()

for(i in 1:nrow(info))
  Std$Sample[which(Std$Sample == info$Condition[i])] <- info$Sample[i]

info_std <- info
info_std <- info_std[-which(grepl("250K", info_std$Sample)),]
pdf(paste0("QC-ISTDs-", Title, ".pdf"), width = 14, height = 10)
ggplot(Std, aes(Sample, log(Value,10)))+
  geom_point(size=3)+
  facet_wrap(~ID, scales='free')+
  theme_bw()+
  theme(text = element_text(face='bold',size=14),
        axis.text.x = element_text(angle=90, vjust=.3, hjust=1))+
  labs(x='',y='Relative response of ISTD (A.U., log10)') + scale_x_discrete(limits = info_std[order(info_std$Run.Order),]$Sample)
dev.off()


#Normalize 
normalize <- FALSE
if(normalize)
{
  output_file <- paste0(Title, "_raw data table normalized.csv")
  for (i in 1: length(info$Cell.Number))
    for (j in 1: nrow(data)) 
      data[j,i+2] <- data[j,i+2]/info$Cell.Number[i]
} else
{
  output_file <- paste0(Title, "_raw data table unnormalized.csv")
}


data_output <- left_join(x = data, y = Abbrev, by = "Used_ID")
for(i in 1:nrow(data_output))
{
  if(is.na(data_output$Abb[i]))
  {
    if(data_output$Used_ID[i] %in% Abbrev$Synonyms)
    {
      j <- which(Abbrev$Synonyms == data_output$Used_ID[i])
      
      data_output[i, which(colnames(data_output) == 'Synonyms'): which(colnames(data_output) == 'Notes')] <-
        Abbrev[j, which(colnames(Abbrev) == 'Synonyms'): which(colnames(Abbrev) == 'Notes')]
    }
  }
}

for(i in 1:nrow(data_output))
{
  if(is.na(data_output$Abb[i]))
  {
    if(data_output$Used_ID[i] %in% Abbrev$Abb)
    {
      j <- which(Abbrev$Abb == data_output$Used_ID[i])
      
      data_output[i, which(colnames(data_output) == 'Synonyms'): which(colnames(data_output) == 'Notes')] <-
        Abbrev[j, which(colnames(Abbrev) == 'Synonyms'): which(colnames(Abbrev) == 'Notes')]
    }
  }
}

data_output <- data_output %>% rename(Name = Abb)
num <- ncol(data_output)
data_output <- data_output %>% 
  select(Name, Used_ID, KEGG.ID,Nr.C, Iso, everything())

data_output[,c('Synonyms', 'Formula', 'Nr.N', 'Polarity', 'Vanq.Method', 'ICS.Method', 'Nitrogens.Method', 'Notes')] <- NULL

data(anno)
for (i in 1:nrow(Anno))
  data_output$Used_ID[as.character(data_output$Name) == Anno$abbrev[i]] <- as.character(Anno$full_name[i])
num_samples <- sum(!grepl('QC',info$Sample))

data_output <- data_output[,c(1:(num_samples + 5), grep(pattern = "QC.blank", x = colnames(data_output)), 
                              grep(pattern = "QC.250K|QC.50K", x = colnames(data_output)), grep(pattern = "QC_|QC.0", x = colnames(data_output)))]
combined_sample_condition <- c()

cond <- c(rep(NA, 5), as.vector(info$Condition[1:num_samples]), rep(NA, (ncol(data_output) - 5 - num_samples)))
cond <- paste(colnames(data_output), cond, sep = " - ")
cond <- gsub(" - NA", "", cond)

colnames(data_output) <- cond
data_output[,6:ncol(data_output)] <- round(data_output[,6:ncol(data_output)])
setwd(output_dir)
write.csv(x = data_output, file = output_file, row.names = F)

data <- data[ , -which(grepl('QC.', colnames(data)))]

#getting sample conditions
info <- info[-which(grepl("QC.", info$Sample)), ]
Freq <- data.frame(table(info$Condition))
Freq <- Freq[Freq$Freq != 0, ]
Sample.Name <- vector()
for (i in 1:nrow(Freq)){
  for (j in 1:Freq$Freq[i]){
    Sample.Name <- append(Sample.Name, paste(Freq$Var1[i], j, sep='_'))
  }
}
info$Sample.Name <- Sample.Name
num_conditions <- length(unique(info$Condition))

new_colnames <- info$Sample
non_blanks_names <- new_colnames[!grepl(paste0("blank",collapse="|"),new_colnames)]
blanks_names <- new_colnames[grepl(paste0("blank",collapse="|"),new_colnames)]
col.order <- c("Used_ID","Iso",non_blanks_names,blanks_names)
#col.order <- gsub("\\.","-",col.order)
data <- data[,col.order]
for (i in 1:length(info$Sample.Name)) 
  colnames(data)[i+2] <- info$Sample.Name[i]
info_ordered <- info[order(info$Run.Order),]

#### Added code to produce QC Norvaline plot
#Creating Norvaline Plot
setwd(output_dir)
if(!ICS)
{
  Norv <- data %>%
    filter(Used_ID=='Norvaline') %>%
    gather(Sample, Norv,-Used_ID, -Iso) %>%
    .$Norv
  
  pdf(file = paste0("QC-Norvaline-", Title, ".pdf"), width=12, height=9, pointsize=12)
  norv_plot<-data %>%
    filter(Used_ID=='Norvaline') %>%
    gather(Sample, Value,-Used_ID, -Iso) %>%
    mutate(Sample=factor(Sample, levels=info_ordered$Sample.Name)) %>%  
    dplyr::select(Sample, Value) %>%
    ggplot(., aes(Sample, Value)) +
    geom_point(size=3) +
    geom_line(aes(as.integer(Sample), Value), color='blue') +
    labs(x='Sample Name',y='Response',title='Norvaline Response') +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + scale_x_discrete(limit = info_ordered$Sample.Name)
  print(norv_plot)  #print statement is required to plot within if statements, loops, functions
  dev.off()
}
####


norm <- data
norm <- norm[!norm$Used_ID=="Norvaline",] 


data2 <- norm %>%
  gather(Exp, Value, -Used_ID, -Iso) %>%
  separate(Exp, c('Condition', 'Exp'), sep='_') %>%
  left_join(., Abbrev, by='Used_ID') %>% 
  dplyr::rename(Name = Abb) %>% 
  dplyr::select(Name, Iso, Condition, Exp, Value, KEGG.ID, Nr.C, everything())
data2$Condition <- factor(data2$Condition, levels=as.character(unique(data2$Condition)))   
data2$Exp <- paste0('Exp', data2$Exp)
NA_Names <- data2 %>%
  group_by(Name) %>%
  mutate(NA_list=sum(is.na(Value)),
         NA_potential=n()) %>%
  filter(NA_list==NA_potential) %>%
  distinct()
data2 <- data2 %>%
  filter(!(Name %in% NA_Names$Name))
data2$Name <- as.character(data2$Name)
data2 <- data2[!is.na(data2$Name),]


select <- dplyr::select
#### RelAmounts

data4 <- data2 %>%
  select(Name, Condition, Iso, KEGG.ID, Exp, Value) %>%
  mutate(Exp = factor(Exp, levels = unique(Exp))) %>%
  group_by(Name, Condition, Exp) %>%
  mutate(Amount=sum(Value, na.rm=T)) %>%
  ungroup() 

data4$ID <- paste(data4$Name, paste(data4$Condition, data4$Exp, sep = "-"),sep = "-")
data4 <- data4[!duplicated(data4$ID), ]
data4$ID <- NULL
data4 <- data4 %>%
  select(Name, Condition, KEGG.ID, Exp, Amount) %>% 
  spread(Exp, Amount)

ATP_ADP=try(cbind(Name="ADP/ATP",data4[1:length(levels(data4$Condition)),2], 'KEGG.ID'=NA,
                  data4[data4$Name=="ADP",4:length(data4)]/data4[data4$Name=="ATP",4:length(data4)]),
            silent=T)
if (exists('ATP_ADP')==T & class(ATP_ADP) != 'try-error') data4 <- rbind(data4, ATP_ADP)
ATP_AMP=try(cbind(Name="AMP/ATP",data4[1:length(levels(data4$Condition)),2], 'KEGG.ID'=NA,
                  data4[data4$Name=="AMP",4:length(data4)]/data4[data4$Name=="ATP",4:length(data4)]),
            silent = T)
if (exists('ATP_AMP')==T & class(ATP_AMP) != 'try-error') data4 <- rbind(data4, ATP_AMP)
GSH_GSSG=try(cbind(Name="GSH/GSSG",data4[1:length(levels(data4$Condition)),2], 'KEGG.ID'=NA,
                   data4[data4$Name=="GSH",4:length(data4)]/data4[data4$Name=="GSSG",4:length(data4)]),
             silent = T)
if (exists('GSH_GSSG')==T & class(GSH_GSSG) != 'try-error') data4 <- rbind(data4, GSH_GSSG)
Creatine_PCreatine=try(cbind(Name="Creatine/P-Creatine",data4[1:length(levels(data4$Condition)),2], 'KEGG.ID'=NA,
                             data4[data4$Name=="Creatine",4:length(data4)]/data4[data4$Name=="P-Creatine",4:length(data4)]),
                       silent = T)
if (exists('Creatine_PCreatine')==T & class(Creatine_PCreatine) != 'try-error') data4 <- rbind(data4, Creatine_PCreatine)

data4 <- data4 %>%
  gather(Exp, Amount, -Name, -Condition,-KEGG.ID)
data4$Amount <- suppressMessages(mapvalues(data4$Amount, c('Inf','NaN'), c(NA,NA)))
data4$Amount[data4$Amount==0] <- NA
test1 <- split(data4, data4[c('Name','Condition')])
NA_function <- function(x) if (sum(is.na(x)) < length(x)) return(x=x) else return(x=rep(0, length(x)))
new.Value <- as.vector(sapply(test1, function(x) NA_function(x$Amount)))

data4 <- suppressWarnings(data4 %>%
                            arrange(Condition, Name) %>%
                            mutate(Amount=new.Value) %>%
                            group_by(Name, Condition) %>%
                            mutate(Av=mean(Amount, na.rm=T),
                                   Std=sd(Amount, na.rm=T),
                                   CV=sd(Amount, na.rm=T)/mean(Amount, na.rm=T)) %>%
                            ungroup())
data8=split(data4, data4[,1])
ANOVA=suppressWarnings(sapply(data8, function(x) anova(aov(x$Amount~x$Condition))$Pr[1]))
ANOVA=rep(ANOVA,1,each=length(levels(data4$Condition)))

data4 <- data4 %>%
  select(Name, KEGG.ID, Condition, Exp, Amount, Av, Std, CV) %>%
  spread(Exp, Amount) %>%
  arrange(Name, Condition) %>%
  cbind('ANOVA'=ANOVA, 'Sig'=NA) %>%
  group_by(Name) %>%
  mutate(Norm_Av=Av/Av[1],
         Norm_Std=Std/Av[1]) %>%
  ungroup()
for (i in 1:nrow(data4))
{
  if (is.na(data4$ANOVA[i])==T) data4$Sig[i]=""
  else if (data4$ANOVA[i] <= 0.001) data4$Sig[i]="***"
  else if (data4$ANOVA[i] <= 0.01) data4$Sig[i]="**"
  else if (data4$ANOVA[i] <= 0.05) data4$Sig[i]="*"
  else data4$Sig[i]=""
}

#### end of RelAmounts
#info$Cell.Number <- rep("",length(info$Cell.Number))
setwd(output_dir)
write.csv(data4, file = paste0(Title, "-Amounts.csv"), row.names = FALSE)

RelA <- make_matrix(data4)
colors <- c("turquoise","red","plum4","steelblue1","red4","springgreen2","slateblue2","sienna1","darkgreen","lightpink1","navy","olivedrab1",
            "orangered","darkslateblue","lightseagreen","magenta2","royalblue","yellowgreen","lightsalmon","cyan","maroon1","indianred3","mediumseagreen",
            "slateblue3","hotpink","lemonchiffon1","orangered4","lightcoral","tomato")


#make_heatmap(RelA, info)
ext <- 'Relative Amounts'
matrix <- RelA
if (exists('samples')==F) samples <- info
#ann <- select(samples, Condition, Cell.Number) %>%
  #as.data.frame()
#ann <- samples[1:6,] %>% select(Sample, Condition, Cell.Number) %>% as.data.frame()
ann <- samples[1:6,] %>% select(Condition, Cell.Number) %>% as.data.frame()
rownames(ann) <- colnames(matrix)
ann$Norvaline <- 1
ann$Cell.Number <- 1
ann_colors = list(
  Condition=colors[1:length(unique(gsub('_(.)*','',colnames(matrix))))],
  Norvaline = c("white", "purple"),
  Cell.Number = c("white", "green")
)
names(ann_colors[['Condition']]) <- unique(gsub('_(.)*','',colnames(matrix)))
matrix[is.na(matrix)] <- 0
heatmap.title=paste(Title, '-Heatmap-',ext,'.pdf', sep='')

pheatmap::pheatmap(matrix, cluster_row = T, cluster_col = T,
                   clustering_distance_rows='correlation',
                   clustering_distance_cols='correlation',
                   color = colorRampPalette(normal)(100),
                   border_color="black", scale="row",
                   cellwidth = 20, cellheight = 10,
                   annotation=ann, annotation_colors = ann_colors,
                   show_colnames = F, main=paste(Title,ext,sep='-'),
                   filename=heatmap.title, width = 9, height = 11)
dev.set(dev.next())

heatmap.title <- paste0(Title, "-Heatmap-", ext, "-Unclustered.pdf")
pheatmap::pheatmap(matrix, cluster_row = T, cluster_col = F,
                   clustering_distance_rows='correlation',
                   clustering_distance_cols='correlation',
                   color = colorRampPalette(normal)(100),
                   border_color="black", scale="row",
                   cellwidth = 20, cellheight = 10,
                   annotation=ann, annotation_colors = ann_colors,
                   show_colnames = F, main=paste(Title,ext,sep='-'),
                   filename=heatmap.title, width = 9, height = 11)
dev.set(dev.next())

#### end of make_heatmap

make_PCA2_update(RelA)
names(data4)[names(data4) == "Norm_Av"] <- "RelAmounts_Ave"
names(data4)[names(data4) == "Norm_Std"] <- "RelAmounts_Std"

## make_MID(data2)
if(tracer_type == "partial" | tracer_type == "full")
{
  data2 <- data2 %>% ungroup()
  percent_label <- data2 %>%
    select(Name, KEGG.ID, Condition, Iso, Nr.C, Exp, Value) 
  
  MID <- percent_label %>%
    group_by(Name, Condition, Exp) %>% 
    arrange(Name, Condition, Exp) %>%
    mutate(Sum = sum(Value, na.rm = T), Fraction = Value * 100 / Sum) %>% 
    ungroup()
  
  MID <- MID %>% 
    group_by(Name, Condition, Iso) %>%
    mutate(Norm_Av=mean(Fraction, na.rm=T),
           Norm_Std=sd(Fraction, na.rm=T),
           CV=sd(Fraction, na.rm=T)/mean(Fraction, na.rm=T),
           Av = Norm_Av) %>% 
    ungroup() %>% 
    select(Name, Condition, Iso, Exp, Fraction, Norm_Av, Norm_Std, CV, Av, Nr.C)
  MID$Exp <- gsub('Exp','MID', MID$Exp)
  MID$Fraction[is.na(MID$Fraction)] <- 0
  MID <- MID[!duplicated(MID[,1:4]), ]
  
  data8 <- split(MID, MID[,c(3,1)], drop=TRUE)
  
  
  ANOVA <- sapply(data8, function(x) {if(sum(x$Fraction, na.rm=T)==0) return(NA) else {anova(aov(x$Fraction~x$Condition))$Pr[1]}})
  
  
  data3 <- MID %>%
    spread( Exp, Fraction) %>% 
    left_join(., percent_label, by = c("Name", "Condition", "Iso")) %>% 
    arrange(Name, Iso)
  ### Rearranged columns so Exp is in the 4th column
  data3 <- data3 %>% relocate(Exp, .after = Iso)
  ### data3.dup <- data3[duplicated(data3[,c(1:4)]), ]
  data3 <- data3[!duplicated(data3[,c(1:4)]), ]
  
  ###data3$ANOVA <- rep(ANOVA, each = num_conditions)
  data3$ANOVA <- rep(ANOVA, each = num_samples)
  for (i in 1:nrow(data3)){
    if (is.na(data3$ANOVA[i]) | data3$ANOVA[i] == "NaN") data3$Sig[i]=""
    else if (data3$ANOVA[i] <= 0.001) data3$Sig[i]="***"
    else if (data3$ANOVA[i] <= 0.01) data3$Sig[i]="**"
    else if (data3$ANOVA[i] <= 0.05) data3$Sig[i]="*"
    else data3$Sig[i]=""
  }
  
  ### These lines were originally commented out
  setwd(output_dir)
  #write.csv(data3, file=paste0(Title,"-Isotopomer data.csv"), row.names=FALSE)
  write.csv(data3, file=paste0(Title,"-Isotopologue data uncorrected.csv"), row.names=FALSE)
  ##write.csv(data3, file=paste0(Title, '-uncorrected_MID.csv'), row.names=F)
  ## end of make_MID
}

library(stringr)
library(writexl)
library(data.table)
# molecule_C$Molecule[molecule_C$Molecule == "5M-adenosine"] <- "5M-thioadenosine"
##CORRECTION (using IsoCorrectoR)
if(tracer_type == "partial" | tracer_type == "full")
{
  library(IsoCorrectoR)
  temp <- correct_iso2(data3, correct_for, correct_for)
  ### C13-10 Iso got turned into C13-1C12 PARENT in the correct_iso function
  ## Need to fix correct_iso gsub
  temp$Iso[temp$Iso == "C13-1C12 PARENT"] <- "C13-10"
  
  temp$Value <- NULL
  colnames(temp)[which(colnames(temp) == "Corrected_Value")] <- "Value"
  temp <- temp %>% 
    select(Name, KEGG.ID, Condition, Iso, Nr.C.x, Exp, Value)
  colnames(temp)[which(colnames(temp) == "Nr.C.x")] <- "Nr.C"   #temp looks like percent labeled.
  
  MID_corrected <- temp %>%
    group_by(Name, Condition, Exp) %>% 
    arrange(Name, Condition, Exp) %>%
    mutate(Sum = sum(Value, na.rm = T), Fraction = Value * 100 / Sum) %>% 
    ungroup()
  MID_corrected <- MID_corrected %>% 
    group_by(Name, Condition, Iso) %>%
    mutate(Norm_Av=mean(Fraction, na.rm=T),
           Norm_Std=sd(Fraction, na.rm=T),
           CV=sd(Fraction, na.rm=T)/mean(Fraction, na.rm=T),
           Av = Norm_Av) %>% 
    ungroup() %>% 
    select(Name, Condition, Iso, Exp, Fraction, Norm_Av, Norm_Std, CV, Av, Nr.C)
  MID_corrected$Exp <- gsub('Exp','MID', MID_corrected$Exp)
  MID_corrected$Fraction[is.na(MID_corrected$Fraction)] <- 0
  
  data8 <- split(MID_corrected, MID_corrected[,c(3,1)], drop=TRUE)
  
  
  ANOVA <- sapply(data8, function(x) {if(sum(x$Fraction, na.rm=T)==0) return(NA) else {anova(aov(x$Fraction~x$Condition))$Pr[1]}})
  
  
  
  data3 <- MID_corrected %>%
    spread(Exp, Fraction) %>% 
    left_join(., temp, by = c("Name", "Condition", "Iso")) %>%     ###changed to temp 
    arrange(Name, Iso)
  
  #data3$ANOVA <- rep(ANOVA, each = num_conditions)
  data3$ANOVA <- rep(ANOVA, each = num_samples)
  for (i in 1:nrow(data3)){
    if (is.na(data3$ANOVA[i]) | data3$ANOVA[i] == "NaN") data3$Sig[i]=""
    else if (data3$ANOVA[i] <= 0.001) data3$Sig[i]="***"
    else if (data3$ANOVA[i] <= 0.01) data3$Sig[i]="**"
    else if (data3$ANOVA[i] <= 0.05) data3$Sig[i]="*"
    else data3$Sig[i]=""
  }
  
  data3$Condition <- factor(data3$Condition, levels = levels(data4$Condition))
  
  #### Added to output corrected MID, which comes from running IsoCorrectR
  #write.csv(data3, file=paste0(Title, '-corrected_MID.csv'), row.names=F)
  write.csv(data3, file=paste0(Title,"-Isotopologue data corrected.csv"), row.names=FALSE)
}

#MIDs heatmaps 
if(tracer_type == "partial" | tracer_type == "full")
{
  # Checking for duplicates
  MS_data <- data3
  anova <- 0.05
  MS_data$ID <- paste(MS_data$Name, MS_data$Condition, MS_data$Iso, sep = "_")
  MS_data <- MS_data[!duplicated(MS_data$ID),]
  MS_data$ID <- NULL
  
  data9 <- MS_data %>%
    #filter(ANOVA <= anova) %>%
    select(Name, Condition, Iso, contains('MID')) %>%
    gather(Exp, Value, -Name, -Condition, -Iso) %>%
    arrange(Name, Condition) %>%
    #mutate(Condition_Exp = Condition:Exp, Condition=NULL, Exp=NULL) %>%
    #mutate(Name_Iso = Name:Iso, Name=NULL, Iso=NULL) %>%
    unite(Condition_Exp, c(Condition, Exp), sep='_') %>%
    mutate(Condition_Exp=factor(Condition_Exp, levels=unique(Condition_Exp))) %>%
    unite(Name_Iso, c(Name, Iso), sep='_') %>%
    spread(Condition_Exp, Value)
  
  Name_Iso <- data9$Name_Iso
  data9 <- data9 %>%
    select(-Name_Iso)%>%
    mutate_if(is.character,as.numeric) 
  data9 <- cbind(Name_Iso, data9)
  
  data9 <- data9[!(rowSums(data9[2:length(data9)], na.rm=T)==0),]
  data9[is.na(data9)] <- 0
  
  data5=as.matrix(data9[2:length(data9)])
  rownames(data5) <- data9$Name
  data5 <- data5[!(rowSums(data5))==0,]
  data5 <- data5[,!(colSums(data5))==0]
  rownames(data5) <- data9$Name_Iso
  MID <- data5
  ## end of make_matrix
  ### Had to change function make_heatmap() into one called make_heatmap2() where
  ## rows with zero variance are taken out before pheatmap is called
  make_heatmap2(MID, info, width = 15, height = 60)
  dev.set(dev.next()) 
  make_heatmap2(MID, info, width = 15, height = 60, cluster_samples = F)
}

## make_data_labeled
if(tracer_type == "partial")
{
  data_labeled <- data3 %>%
    select(Name, KEGG.ID, Condition, Iso, starts_with('MID')) %>%
    filter(Iso=='C12 PARENT') %>%
    gather(Exp, Value, -Name, -KEGG.ID, -Condition, -Iso) 
  data_labeled$Value <- as.numeric(data_labeled$Value)
  data_labeled <- data_labeled %>% 
    mutate(Labeled=(1-Value/100)*100) %>%
    group_by(Name, Condition) %>%
    mutate(Norm_Av=mean(Labeled, na.rm=T),
           Norm_Std=sd(Labeled, na.rm=T),
           CV=Norm_Std/Norm_Av,
           Av = Norm_Av) %>%
    select(-Iso, -Value) %>%
    ungroup()
  data_labeled$Exp <- gsub('MID','Labeled', data_labeled$Exp)
  data_labeled <- unique(data_labeled)
  #####fixing rounding problem 
  data_labeled$Labeled <- round(data_labeled$Labeled, 13)
  
  test1 <- split(data_labeled, data_labeled[c('Name', 'Condition')])
  NA_function <- function(x) if (sum(is.na(x)) < length(x)) return(x=x) else return(x=rep(0, length(x)))
  new.Value <- as.vector(sapply(test1, function(x) NA_function(x$Labeled)))
  #a <- unlist(new.Value)
  data_labeled <- data_labeled %>%
    arrange(Condition, Name) %>%
    mutate(Labeled=new.Value)
  
  data8 <- split(data_labeled, data_labeled[,1])
  ANOVA <- suppressWarnings(sapply(data8, function(x) anova(aov(x$Labeled~x$Condition))$Pr[1]))
  ANOVA <- rep(ANOVA,1,each=length(unique(info$Condition)))
  
  data_labeled <- spread(data_labeled, Exp, Labeled) %>%
    arrange(Name)
  
  data_labeled <- cbind(data_labeled, ANOVA, 'Sig'=NA)
  for (i in 1:nrow(data_labeled)){
    if (data_labeled$ANOVA[i] == "NaN") data_labeled$Sig[i]=""
    else if (data_labeled$ANOVA[i] <= 0.001) data_labeled$Sig[i]="***"
    else if (data_labeled$ANOVA[i] <= 0.01) data_labeled$Sig[i]="**"
    else if (data_labeled$ANOVA[i] <= 0.05) data_labeled$Sig[i]="*"
    else data_labeled$Sig[i]=""
  }
  write.csv(data_labeled, file=paste0(Title,'-labeled data.csv'), row.names=F)
  
  data_labeled_mat <- make_matrix(data_labeled)
  make_heatmap(data_labeled_mat, info, width = 10, height = 10)
  dev.set(dev.next())
  make_heatmap(data_labeled_mat, info, cluster_samples = F, width = 10, height = 7)
  dev.set(dev.next())
}


###make_FC (might need fixing.)
if(tracer_type == "full")
{
  DF <- data3
  
  DF$num <- NA
  ### Here, there are 12 in Df$Iso that are C13-1C12 Parent
  for (i in 1:nrow(DF))
  {
    if(DF$Iso[i] == "C12 PARENT") #Parent
      DF$num[i] <- 0
    else if(grepl("C13-", DF$Iso[i])) #Carbons
    {
      DF$num[i] <- as.numeric(strsplit(DF$Iso[i], "-")[[1]][2])
    }
    else if(grepl("C13N15-", DF$Iso[i])) #Carbons/Nitrogens
    {
      DF$num[i] <- sum(as.numeric(strsplit(DF$Iso[i], "-")[[1]][2:3]))
    }
    else #rest should be Nitrogens
    {
      DF$num[i] <- as.numeric(strsplit(DF$Iso[i], "-")[[1]][2])
    }
  }
  DF$Exp <- NULL
  DF$ANOVA <- NULL
  DF$Sig <- NULL
  DF$Nr.C.y <- NULL
  colnames(DF)[which(colnames(DF) == "Nr.C.x")] <- "Nr.C"
  ### Make Value NULL and remove duplicates
  DF$Value <- NULL
  DF <- DF[!duplicated(DF[,]), ]
  
  FC <- suppressWarnings(DF %>%
                           #    inner_join(., Abbrev, by=c('Name'='Abb', 'KEGG.ID'='KEGG.ID')) %>%
                           select(Name, KEGG.ID, Condition, Iso, Nr.C, starts_with('MID'), num) %>%
                           gather(Exp, MID, -Name, -KEGG.ID, -Condition, -Iso,-Nr.C, -num) %>%
                           mutate(iSi = num*MID) %>%
                           group_by(Name, Condition, Exp) %>%
                           mutate(FC=sum(iSi, na.rm=T)/Nr.C) %>%
                           filter(Iso=='C12 PARENT') %>%
                           ungroup() %>%
                           mutate(FC=mapvalues(FC, 0, NA)) %>%                #this is important when a sample is missing or all MID values were 0
                           select(Name, KEGG.ID, Condition, Exp, FC) %>%
                           group_by(Name, Condition) %>%
                           mutate(Norm_Av=mean(FC, na.rm=T),
                                  Norm_Std=sd(FC, na.rm=T),
                                  CV=Norm_Std/Norm_Av,
                                  Av = Norm_Av) %>%
                           ungroup())
  FC$Exp <- gsub('MID','FC',FC$Exp)
  test1 <- split(FC, FC[c('Name', 'Condition')])
  NA_function <- function(x) if (sum(is.na(x)) < length(x)) return(x=x) else return(x=rep(0, length(x)))
  new.Value <- as.vector(sapply(test1, function(x) NA_function(x$FC)))
  
  FC <- FC %>%
    arrange(Condition, Name) %>%
    mutate(FC=new.Value)
  
  data8=split(FC, FC[,1])
  ANOVA=suppressWarnings(sapply(data8, function(x) anova(aov(x$FC~x$Condition))$Pr[1]))
  
  ANOVA=rep(ANOVA,1,each=length(unique(info$Condition)))
  
  ### This gives an error: Each row of output must be identified by a unique combination of keys.
  ### No longer needed? FC <- FC[!duplicated(FC[,]), ]
  #FC <- FC %>% distinct()
  FC <- spread(FC, Exp, FC) %>%
    arrange(Name) %>%
    mutate(Sig='NA')
  FC$ANOVA <- ANOVA
  
  for (i in 1:nrow(FC)){
    if (FC$ANOVA[i] == "NaN") FC$Sig[i]=""
    else if (FC$ANOVA[i] <= 0.001) FC$Sig[i]="***"
    else if (FC$ANOVA[i] <= 0.01) FC$Sig[i]="**"
    else if (FC$ANOVA[i] <= 0.05) FC$Sig[i]="*"
    else FC$Sig[i]=""
  }
  
  write.csv(FC, file=paste0(Title, '-fractional contribution.csv'), row.names=F)
  
  #make_matrix
  MS_data <- FC
  data9 <- MS_data %>%
    # filter(ANOVA <= anova) %>%
    select(Name, Condition, grep('Exp|FC|Labeled', names(MS_data))) %>%
    gather(Exp, Value, -Name, -Condition) %>%
    arrange(Name, Condition) %>%
    #mutate(Exp = as.numeric(gsub('Exp|FC|Labeled','',.$Exp))) %>%
    group_by(Name) %>%
    mutate(Nr.NA = sum(is.na(Value)),
           Nr.Samples = n()) %>%
    filter(Nr.NA < Nr.Samples - 1) %>%
    ungroup() %>%
    select(-Nr.NA, -Nr.Samples) %>%
    #arrange(Name, Condition, Exp) %>%
    unite(Condition_Exp, c(Condition, Exp), sep='_') %>%
    mutate(Condition_Exp=factor(Condition_Exp, levels=unique(Condition_Exp))) %>%
    spread(Condition_Exp, Value)
  data9[is.na(data9)] <- 0
  data5=as.matrix(data9[2:length(data9)])
  rownames(data5) <- data9$Name
  data5 <- data5[!(rowSums(data5))==0,]
  data5 <- data5[,!(colSums(data5))==0]
  
  FC_mat <- data5
  
  #FC heatmap
  ext = 'Fractional Contribution'
  matrix <- FC_mat
  if (exists('samples')==F) samples <- info
  ann <- select(samples, Condition, Cell.Number) %>%
    as.data.frame()
  rownames(ann) <- colnames(matrix)
  
  if (exists('Norv')==T) {
    if(all(is.na(Norv)))
      Norv <- rep(1, length(Norv))
    ann$Norvaline <- Norv
  } else ann$Norvaline <- 1
  
  ann_colors = list(
    Condition=colors[1:length(unique(gsub('_(.)*','',colnames(matrix))))],
    Norvaline = c("white", "purple"),
    Cell.Number = c("white", "green")
  )
  names(ann_colors[['Condition']]) <- unique(gsub('_(.)*','',colnames(matrix)))
  matrix[is.na(matrix)] <- 0
  heatmap.title=paste(Title, '-Heatmap-',ext, '-Unclustered', '.pdf', sep='')
  
  ann$Cell.Number <- rep(1, 6)
  # Might need to adjust width and height of heatmap here
  pheatmap::pheatmap(matrix, cluster_row = T, cluster_col = F,
                     clustering_distance_rows='correlation',
                     clustering_distance_cols='correlation',
                     color = colorRampPalette(normal)(100),
                     border_color="black", scale="row",
                     cellwidth = 20, cellheight = 10,
                     annotation=ann, annotation_colors = ann_colors,
                     show_colnames = F, main=paste(Title,ext,sep='-'),
                     filename=heatmap.title, width = 10, height = 20)
}


CoAs <- c(CoAs, "Isobutyryl-CoA") 
pathway <- c(glycolysis, TCA, PPP, Curr, Cys, AA, FA, Hex, Adenine, Cytosine, Guanine, Thymine, Uracil, Fru, CoAs)
all_data4 <- unique(data4$Name)
non_data4 <- setdiff(all_data4, pathway)
all_data3 <- unique(data3$Name)
non_data3 <- setdiff(all_data3, pathway)
all_data_labeled <- unique(data_labeled$Name)
non_data_labeled <- setdiff(all_data_labeled, pathway)
all_FC <- unique(FC$Name)
non_FC <- setdiff(all_FC, pathway)

non_data4 <- split(non_data4, ceiling(seq_along(non_data4)/20))
non_data3 <- split(non_data3, ceiling(seq_along(non_data3)/20))
non_data_labeled <- split(non_data_labeled, ceiling(seq_along(non_data_labeled)/20))
non_FC <- split(non_FC, ceiling(seq_along(non_FC)/20))

temp <- (data %>% group_by(Used_ID) %>% filter(n() == 1))
temp <- left_join(temp, Abbrev, by = 'Used_ID')
only_M0 <- temp$Abb

names(data4)[names(data4) == "Norm_Av"] <- "RelAmounts_Ave"
names(data4)[names(data4) == "Norm_Std"] <- "RelAmounts_Std"
names(data3)[names(data3) == "Norm_Av"] <- "RelAmounts_Ave"
names(data3)[names(data3) == "Norm_Std"] <- "RelAmounts_Std"
names(data_labeled)[names(data_labeled) == "Norm_Av"] <- "RelAmounts_Ave"
names(data_labeled)[names(data_labeled) == "Norm_Std"] <- "RelAmounts_Std"
names(FC)[names(FC) == "Norm_Av"] <- "RelAmounts_Ave"
names(FC)[names(FC) == "Norm_Std"] <- "RelAmounts_Std"

data3$ID <- paste(data3$Name, data3$Condition, data3$Iso, sep = "_")
data3 <- data3[!duplicated(data3$ID), ]
data3$ID <- NULL


data3 <- check_50(data3)
data4 <- check_50(data4)
data_labeled <- check_50(data_labeled)
FC <- check_50(FC)

for(i in 1:nrow(data3))
{
  if(data3$under_50_percent[i] == "X")
    data3$RelAmounts_Std[i] <- "0.00"
}

for(i in 1:nrow(data_labeled))
{
  if(data_labeled$under_50_percent[i] == "X")
    data_labeled$RelAmounts_Std[i] <- 0
}

data3$RelAmounts_Ave <- as.numeric(data3$RelAmounts_Ave)
data3$RelAmounts_Std <- as.numeric(data3$RelAmounts_Std)


metabs <- unique(data4$Name)
for(i in metabs)
{
  if(all(is.nan(data4[which(data4$Name == i), 'RelAmounts_Ave'])))
    data4 <- data4[-which(data4$Name == i), ]
}


#MIDs are uncorrected
setwd(output_dir)
plotname <- paste0(Title, "-Plots All.pdf")
if(tracer_type == "none")
  plotname <- paste0(Title, "-Plots RelAmounts.pdf")
pdf(file = plotname, width=14, height=10, pointsize=12)

bar_update_manual('glycolysis',data4, n = num_conditions, type = "tf")
Maven_MID_plot('glycolysis',data3, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot('glycolysis',data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot('glycolysis',FC, n = num_conditions, type = "tf", only_M0 = only_M0)

bar_update_manual("TCA",data4, n = num_conditions, type = "tf")
Maven_MID_plot("TCA",data3, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot("TCA",data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot("TCA",FC, n = num_conditions, type = "tf", only_M0 = only_M0)

bar_update_manual("PPP",data4, n = num_conditions, type = "tf")
Maven_MID_plot("PPP",data3, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot("PPP",data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot("PPP",FC, n = num_conditions, type = "tf", only_M0 = only_M0)

bar_update_manual("Curr",data4, n = num_conditions, type = "tf")
Maven_MID_plot("Curr",data3, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot("Curr",data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot("Curr",FC, n = num_conditions, type = "tf", only_M0 = only_M0)

bar_update_manual("Cys",data4, n = num_conditions, type = "tf")
Maven_MID_plot("Cys",data3, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot("Cys",data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot("Cys",FC, n = num_conditions, type = "tf", only_M0 = only_M0)

bar_update_manual("AA",data4, n = num_conditions, type = "tf")
Maven_MID_plot("AA",data3, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot("AA",data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot("AA",FC, n = num_conditions, type = "tf", only_M0 = only_M0)

bar_update_manual("FA",data4, n = num_conditions, type = "tf")
Maven_MID_plot("FA",data3, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot("FA",data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot("FA",FC, n = num_conditions, type = "tf", only_M0 = only_M0)

bar_update_manual("Hex",data4, n = num_conditions, type = "tf")
Maven_MID_plot("Hex",data3, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot("Hex",data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot("Hex",FC, n = num_conditions, type = "tf", only_M0 = only_M0)

bar_update_manual("Adenine",data4, n = num_conditions, type = "tf")
Maven_MID_plot("Adenine",data3, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot("Adenine",data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot("Adenine",FC, n = num_conditions, type = "tf", only_M0 = only_M0)

bar_update_manual("Cytosine",data4, n = num_conditions, type = "tf")
Maven_MID_plot("Cytosine",data3, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot("Cytosine",data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot("Cytosine",FC, n = num_conditions, type = "tf", only_M0 = only_M0)

bar_update_manual("Guanine",data4, n = num_conditions, type = "tf")
Maven_MID_plot("Guanine",data3, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot("GUanine",data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot("Guanine",FC, n = num_conditions, type = "tf", only_M0 = only_M0)

bar_update_manual("Thymine",data4, n = num_conditions, type = "tf")
Maven_MID_plot("Thymine",data3, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot("Thymine",data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot("Thymine",FC, n = num_conditions, type = "tf", only_M0 = only_M0)

bar_update_manual("Uracil",data4, n = num_conditions, type = "tf")
Maven_MID_plot("Uracil",data3, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot("Uracil",data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot("Uracil",FC, n = num_conditions, type = "tf", only_M0 = only_M0)

bar_update_manual('Fru',data4, n = num_conditions, type = "tf")
Maven_MID_plot('Fru',data3, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot("Fru",data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot("Fru",FC, n = num_conditions, type = "tf", only_M0 = only_M0)

bar_update_manual('CoAs',data4, n = num_conditions, type = "tf")
Maven_MID_plot('CoAs',data3, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot("CoAs",data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
Maven_MID_plot("CoAs",FC, n = num_conditions, type = "tf", only_M0 = only_M0)
#first set
bar_update_manual(unname(unlist(non_data4[1])),data4, n = num_conditions, type = "tf", title_type = "nonpathway1")
Maven_MID_plot(unname(unlist(non_data3[1])),data3, n = num_conditions, type = "tf", title_type = "nonpathway1", only_M0 = only_M0)
Maven_MID_plot(unname(unlist(non_data3[1])),data_labeled, n = num_conditions, type = "tf", title_type = "nonpathway1", only_M0 = only_M0)
Maven_MID_plot(unname(unlist(non_data3[1])),FC, n = num_conditions, type = "tf", title_type = "nonpathway1", only_M0 = only_M0)

#second set
bar_update_manual(unname(unlist(non_data4[2])),data4, n = num_conditions, type = "tf", title_type = "nonpathway2")
Maven_MID_plot(unname(unlist(non_data3[2])),data3, n = num_conditions, type = "tf", title_type = "nonpathway2", only_M0 = only_M0)
Maven_MID_plot(unname(unlist(non_data3[2])),data_labeled, n = num_conditions, type = "tf", title_type = "nonpathway2", only_M0 = only_M0)
Maven_MID_plot(unname(unlist(non_data3[2])),FC, n = num_conditions, type = "tf", title_type = "nonpathway2", only_M0 = only_M0)

#third set
bar_update_manual(unname(unlist(non_data4[3])),data4, n = num_conditions, type = "tf", title_type = "nonpathway3")
Maven_MID_plot(unname(unlist(non_data3[3])),data3, n = num_conditions, type = "tf", title_type = "nonpathway3", only_M0 = only_M0)
Maven_MID_plot(unname(unlist(non_data3[3])),data_labeled, n = num_conditions, type = "tf", title_type = "nonpathway3", only_M0 = only_M0)
Maven_MID_plot(unname(unlist(non_data3[3])),FC, n = num_conditions, type = "tf", title_type = "nonpathway3", only_M0 = only_M0)

dev.off()

