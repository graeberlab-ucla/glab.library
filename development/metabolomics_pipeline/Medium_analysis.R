#Sys.setenv(JAVA_HOME = "C:\\Program Files\\Java\\jdk-11.0.8\\")
library(MetabFUN)
library(tidyr)
library(readxl)
library(dplyr)
#devtools::install_github("juyeki/MetabR")
library(MetabR)
library(xlsx)
library(pheatmap)
####
library(ggplot2)

data_dir <- "N:/TraceFinderData/4.0/Projects/Medium/Goldstein Lab/JG-07222020_7-7_medium_Vanq"
output_dir <- "N:/TraceFinderData/4.0/Projects/Medium/Goldstein Lab/JG-07222020_7-7_medium_Vanq/test2"
#data_dir <- "/Users/tgraeber/Dropbox/glab/Metabolomics Center/Metabolomics Pipeline Scripts and Notes/footprinting troubleshooting/more than one unspent control set/2020.07.07_Jenna-Godstein_16D Enza_metabolic footprinting"
#output_dir <- "/Users/tgraeber/Dropbox/glab/Metabolomics Center/Metabolomics Pipeline Scripts and Notes/footprinting troubleshooting/more than one unspent control set/2020.07.07_Jenna-Godstein_16D Enza_metabolic footprinting/test2"
## This is the only place anno_dir is in the script
###anno_dir <- "C:/Users/djtan/Documents/metabolobics scripts"
unlabelled <- TRUE
tracer_type<-"none"

colors1<-c("turquoise","red","plum4","steelblue1","red4","springgreen2","slateblue2","sienna1","darkgreen","lightpink1","navy","olivedrab1",
           "orangered","darkslateblue","lightseagreen","magenta2","royalblue","yellowgreen","lightsalmon","cyan","maroon1","indianred3","mediumseagreen",
           "slateblue3","hotpink","lemonchiffon1","orangered4","lightcoral","tomato")
####data(abbrev_vanq)
#### Changed abbrev from one in package to the one in google drive
#### abbrev_vanq vs abbrev?
Abbrev <- read_excel("C:/Users/FTsang/Downloads/Abbrev_NEW2 (1).xlsx")
#Abbrev <- read_excel("/Users/tgraeber/Dropbox/glab/Metabolomics Center/Metabolomics Pipeline Scripts and Notes/workspace/Abbrev_NEW2.xlsx")
Abbrev <- Abbrev[-c(248,249),]

setwd(data_dir)
Title <- paste0('Footprint-',gsub('.xls[x]?','', list.files(pattern='.xls[x]?')))
#Title<-Title[2] #this approach should be improved
Title<-Title[1] #this approach should be improved


# enter cell names and cell numbers
#info <- read_excel(list.files()[grep('.xls[x]?',list.files())][2]) #this approach should be improved
info <- read_excel(list.files()[grep('.xls[x]?',list.files())][1]) #this approach should be improved
info$Samples <- 1:nrow(info)
info$Condition <- gsub("^\\s+|\\s+$", "", info$Condition)        #remove leading or trailing white space
info$Condition <- gsub('/', '-', info$Condition)
info$Condition <- gsub('_', '-', info$Condition)
info$Sample <- gsub('-', '.', info$Sample)
info$Condition <- factor(info$Condition, levels = unique(info$Condition)) 

samples <- info$Sample
samples <- samples[!grepl(paste0("QC", collapse = "|"), samples)]
num_samples<-length(samples)
################################################################################



setwd(data_dir)
folder.name <- gsub('(.)*Projects/','',data_dir) %>% gsub('/','-',.)
files <- "Medium/Medium.csv"


data <- lapply(files, read.csv, header=T, na.strings = c("N/F",'N/A'))
data <- do.call(rbind, data)
data$Experiment <- folder.name

#setwd("C:/Users/gjaviersanchez/Documents/Metabolomics/Temp2")
setwd(output_dir)
write.csv(data, 'all data raw.csv', row.names=F)

data <- data %>%
  dplyr::select(Compound, Filename, Area) %>%
  spread(Filename, Area)

keep <- gsub("\\.", "-", info$Sample)
keep <- c("Compound", keep)
data <- data[, which(colnames(data) %in%  keep)]
data_temp <- filter_data(data, info)

data$Compound = as.character(data$Compound)
for (i in 1:length(data$Compound)){                   #loop through IDs and fill up with Ms
  if (grepl(" M[0-9]+$", data$Compound[i])==FALSE) data$Compound[i]=paste(data$Compound[i],"M0")}
Used_ID <- gsub(" M[0-9]+$| Std", "", data$Compound)
name_split <- regexpr("M[0-9]+$|Std",data$Compound)
Iso <-regmatches(data$Compound, name_split)
data <- data.frame(Used_ID, Iso, data[,2:length(data)])
data <- data %>% select(-grep('QC_|QC.', colnames(data)), grep('QC_|QC.', colnames(data)))

#data$Iso <- factor(data$Iso, levels=as.character(unique(data$Iso)))

data_output <- left_join(x = data, y = Abbrev, by = "Used_ID")
data_output <- data_output %>% rename(Name = Abb)
#data_output <- data_output %>% select(Name, Used_ID, KEGG.ID, Pathway, Nr.C, Iso, info$Sample)
data_output <- data_output %>% select(Name, Used_ID, KEGG.ID, Vanq.Method, Nr.C, Iso, info$Sample)

if(tracer_type == "none")
  data_output <- subset(data_output, Iso == "M0" | Iso =="M1")
####################################
### Took out rows with NAs for Nr.C
data_output <- data_output[!is.na(data_output$Nr.C),]
isotopes <- rep(NA, max(data_output$Nr.C))
for (i in 1:(length(isotopes)+1))
  isotopes[i] <- paste0("M", as.character(i-1))
data_output <- data_output %>% 
  mutate(Iso = factor(Iso, levels = isotopes)) %>%
  arrange(Name, Iso)
num_samples <- sum(!grepl('QC',info$Sample))

data(anno)
for (i in 1:nrow(Anno))
  data_output$Used_ID[as.character(data_output$Name) == Anno$abbrev[i]] <- as.character(Anno$full_name[i])
levels(data_output$Iso) <- c(levels(data_output$Iso), "Std")
data_output$Iso[ which(is.na(data_output$Iso))] <- "Std"
data_output <- rbind(c(rep(NA, 6), as.vector(info$Condition[1:num_samples]), rep(NA, (ncol(data_output) - 6 - num_samples))), data_output)
data_output <- data_output[,c(1:(num_samples + 6), grep(pattern = "QC.blank", x = colnames(data_output)), 
                              grep(pattern = "QC.250K|QC.50K", x = colnames(data_output)), grep(pattern = "QC_|QC.0", x = colnames(data_output)))]

combined_sample_condition <- c()
for (i in 7:(num_samples+6))
{
  x <- paste(colnames(data_output)[i], data_output[1,i], sep=" - ")
  combined_sample_condition <- c(combined_sample_condition, x)
}

colnames(data_output)[7:(num_samples+6)]<-combined_sample_condition
data_output<-data_output[-1,]
setwd(output_dir)
write.xlsx(x = data_output, file = paste0(Title, "_raw data table.xlsx"), row.names = F, col.names = T, showNA = F)

####################################
info_with_blanks<-info
info <- info[!grepl("QC.blank", info$Sample), ] 
data <- data_temp
data_with_blanks<-data
data <- data[ , -which(grepl('QC.blank', colnames(data)))]
data$Compound = as.character(data$Compound)
data_with_blanks$Compound = as.character(data_with_blanks$Compound)
for (i in 1:length(data$Compound)){                   #loop through IDs and fill up with Ms
  if (grepl(" M[0-9]+$", data$Compound[i])==FALSE) data$Compound[i]=paste(data$Compound[i],"M0")}
Used_ID <- gsub(" M[0-9]+$| Std", "", data$Compound)
name_split <- regexpr("M[0-9]+$|Std",data$Compound)
Iso <-regmatches(data$Compound, name_split)
data <- data.frame(Used_ID, Iso, data[,2:length(data)])
data <- data %>% select(-grep('QC_|QC.', colnames(data)), grep('QC_|QC.', colnames(data)))


for (i in 1:length(data_with_blanks$Compound)){                   #loop through IDs and fill up with Ms
  if (grepl(" M[0-9]+$", data_with_blanks$Compound[i])==FALSE) data_with_blanks$Compound[i]=paste(data_with_blanks$Compound[i],"M0")}
Used_ID <- gsub(" M[0-9]+$| Std", "", data_with_blanks$Compound)
name_split <- regexpr("M[0-9]+$|Std",data_with_blanks$Compound)
Iso <-regmatches(data_with_blanks$Compound, name_split)
data_with_blanks <- data.frame(Used_ID, Iso, data_with_blanks[,2:length(data_with_blanks)])
data_with_blanks <- data_with_blanks %>% select(-grep('QC_|QC.', colnames(data_with_blanks)), grep('QC_|QC.', colnames(data_with_blanks)))


std <- make_istd(data, remove = FALSE)$Std
### Remove Norvaline from std plots
std <- std[std$Used_ID != "Norvaline",]
setwd(output_dir)
info_ordered <- info[order(info$Run.Order),]
info_ordered_with_blanks<-info_with_blanks[order(info_with_blanks$Run.Order),]

istd_title <- paste0("QC-ISTDs with 250ks-", Title, ".pdf")
plot_istd(Std = std, info = info_ordered, title = istd_title,pdf_width=20)

data <- data[ , !grepl('QC.', colnames(data))]
data_with_blanks <- data_with_blanks[ , !grepl('QC.250|QC.0', colnames(data_with_blanks))]
info <- info[!grepl("QC.", info$Sample), ]
info_with_blanks <- info_with_blanks[!grepl("QC.250|QC.0", info_with_blanks$Sample), ]
std <- make_istd(data, remove = TRUE)$Std
data <- make_istd(data, remove = TRUE)$data 
setwd(output_dir)
info_ordered <- info[order(info$Run.Order),]
plot_istd(Std = std, info = info_ordered, title = paste0("QC-ISTDs-", Title, ".pdf"),pdf_width = 20)

#removing standards
###data <- filter(data, !(grepl('Std', data$Iso)))
data <- filter(data, !(grepl('Std', data$Iso) & data$Used_ID != "Norvaline"))
data$Used_ID <- as.character(data$Used_ID)
data_with_blanks$Used_ID <- as.character(data_with_blanks$Used_ID)

data$Iso <- factor(data$Iso, levels=paste0('M',0:50))
data_with_blanks$Iso <- factor(data_with_blanks$Iso, levels=paste0('M',0:50))
data <- plyr::arrange(data, Used_ID, Iso)
data_with_blanks <- plyr::arrange(data_with_blanks, Used_ID, Iso)

# Sample.Name<-Sample.Name[!grepl("QC-250",Sample.Name)]
# ####
# col.order<-c("Name","Iso",info$Sample)
# col.order<-gsub("-",".",col.order)
# data<-data[,col.order]
# # check ordering
# ####
# for (i in 1:length(info$Sample.Name)) colnames(data)[i+2] <- info$Sample.Name[i]
# 
# 
# #substitute long names and add numbers of carbons in molecule
# for (i in nrow(Abbrev):1) {
#   data$Name <- sub(Abbrev$Used_ID[i], Abbrev$Abb[i], data$Name, fixed=T)
# }

info_with_blanks$Condition <- factor(info_with_blanks$Condition, levels = unique(info_with_blanks$Condition)) 
Freq <- data.frame(table(info$Condition))
Freq_with_blanks<- data.frame(table(info_with_blanks$Condition))
Freq <- Freq[Freq$Freq != 0, ]
Freq_with_blanks <- Freq_with_blanks[Freq_with_blanks$Freq != 0, ]

Sample.Name <- vector()
Sample.Name_with_blanks<-vector()
for (i in 1:nrow(Freq)){
  for (j in 1:Freq$Freq[i]){
    Sample.Name <- append(Sample.Name, paste(Freq$Var1[i], j, sep='_'))
  }
}

for (i in 1:nrow(Freq_with_blanks)){
  for (j in 1:Freq_with_blanks$Freq[i]){
    Sample.Name_with_blanks <- append(Sample.Name_with_blanks, paste(Freq_with_blanks$Var1[i], j, sep='_'))
  }
}

info$Sample.Name <- Sample.Name
info_with_blanks$Sample.Name <- Sample.Name_with_blanks
num_conditions <- length(unique(info$Condition))
num_conditions_with_blanks <- length(unique(info_with_blanks$Condition))

#adding conditions to data
new_colnames <- info$Sample
new_colnames_with_blanks <- info_with_blanks$Sample
non_blanks_names <- new_colnames[!grepl(paste0("blank",collapse="|"),new_colnames)]
non_blanks_names_with_blanks <- new_colnames_with_blanks[!grepl(paste0("blank",collapse="|"),new_colnames_with_blanks)]
blanks_names <- new_colnames[grepl(paste0("blank",collapse="|"),new_colnames)]
blanks_names_with_blanks <- new_colnames_with_blanks[grepl(paste0("blank",collapse="|"),new_colnames_with_blanks)]
col.order <- c("Used_ID","Iso",non_blanks_names,blanks_names)
col.order_with_blanks <-c("Used_ID","Iso",non_blanks_names_with_blanks,blanks_names_with_blanks)
col.order <- gsub("-",".",col.order)
col.order_with_blanks <- gsub("-",".",col.order_with_blanks)
data <- data[,col.order]
data_with_blanks<- data_with_blanks[,col.order_with_blanks]
for (i in 1:length(info$Sample.Name)) 
  colnames(data)[i+2] <- info$Sample.Name[i]
info_ordered <- info[order(info$Run.Order),]

for (i in 1:length(info_with_blanks$Sample.Name)) 
  colnames(data_with_blanks)[i+2] <- info_with_blanks$Sample.Name[i]
info_ordered_with_blanks <- info_with_blanks[order(info_with_blanks$Run.Order),]

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
dev.off()

#Normalizing Norvaline 
if (length(Norv)==0){
  Norv <- rep(1, length(samples))
  #Norv <- rep(1,nrow(samples))
  print("No norvaline")
}

norm <- data
norm_with_blanks<-data_with_blanks

if (length(Norv)!=0){
  ###norm <- norm[!norm$Name=="Norvaline",]     #remove Norvaline from list
  norm <- norm[!norm$Used_ID=="Norvaline",]
}

if (length(Norv)!=0){
  ###norm_with_blanks <- norm_with_blanks[!norm_with_blanks$Name=="Norvaline",]     #remove Norvaline from list
  norm_with_blanks <- norm_with_blanks[!norm_with_blanks$Used_ID=="Norvaline",]
}

if(length(unique(info$Cell.Number))==1)
{info$Cell.Number=1}

if(length(unique(info_with_blanks$Cell.Number[1:length(new_colnames)]))==1)
{info_with_blanks$Cell.Number=1}

to_normalize <- FALSE
info$Cell.Number <- 1
###to_normalize <- !(all(is.na(info$Cell.Number)))
if(to_normalize)
  if(length(unique(info$Cell.Number)) == 1)
    to_normalize <- FALSE
if (to_normalize & is.numeric(info$Cell.Number)==F) print('Problem: You need to convert the Cell.Number column')
if(to_normalize)
{
  for (i in 1: length(info$Cell.Number)){
    for (j in 1: nrow(norm)) norm[j,i+2]=norm[j,i+2]/info$Cell.Number[i]
  }
}
if(all(is.na(info$Cell.Number)))
  info$Cell.Number <- 1
#####
###to_normalize_with_blanks <- !(all(is.na(info_with_blanks$Cell.Number)))
to_normalize_with_blanks <- FALSE
info_with_blanks$Cell.Number <- 1
if(to_normalize_with_blanks)
  if(length(unique(info_with_blanks$Cell.Number)) == 1)
    to_normalize_with_blanks <- FALSE
if (to_normalize_with_blanks & is.numeric(info_with_blanks$Cell.Number)==F) print('Problem: You need to convert the Cell.Number column')
if(to_normalize_with_blanks)
{
  for (i in 1: length(info_with_blanks$Cell.Number)){
    for (j in 1: nrow(norm_with_blanks)) norm_with_blanks[j,i+2]=norm_with_blanks[j,i+2]/info_with_blanks$Cell.Number[i]
  }
}
if(all(is.na(info_ordered_with_blanks$Cell.Number)))
  info_with_blanks$Cell.Number <- 1

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

data2_with_blanks <- norm_with_blanks %>%
  gather(Exp, Value, -Used_ID, -Iso) %>%
  separate(Exp, c('Condition', 'Exp'), sep='_') %>%
  left_join(., Abbrev, by='Used_ID') %>% 
  dplyr::rename(Name = Abb) %>% 
  dplyr::select(Name, Iso, Condition, Exp, Value, KEGG.ID, Nr.C, everything())
data2_with_blanks$Condition <- factor(data2_with_blanks$Condition, levels=as.character(unique(data2_with_blanks$Condition)))   
data2_with_blanks$Exp <- paste0('Exp', data2_with_blanks$Exp)
NA_Names <- data2_with_blanks %>%
  group_by(Name) %>%
  mutate(NA_list=sum(is.na(Value)),
         NA_potential=n()) %>%
  filter(NA_list==NA_potential) %>%
  distinct()
data2_with_blanks <- data2_with_blanks %>%
  filter(!(Name %in% NA_Names$Name))
data2_with_blanks$Name <- as.character(data2_with_blanks$Name)
data2_with_blanks <- data2_with_blanks[!is.na(data2_with_blanks$Name),]


if(tracer_type == "none")
  data2 <- subset(data2, Iso == "M0" | Iso =="M1")

if(tracer_type == "none")
  data2_with_blanks <- subset(data2_with_blanks, Iso == "M0" | Iso =="M1")

if(to_normalize){
  output_file <- paste0(Title, "_raw data table-Filtered and Normalized.xlsx")
}else
  output_file <- paste0(Title, "_raw data table-Filtered.xlsx")
data2_output <- data2
data2_output$rep <- gsub(pattern = 'Exp', replacement = '', x = data2_output$Exp)
data2_output$Sample <- paste(data2_output$Condition, data2_output$rep, sep = "_")
data2_output$Condition <- NULL
data2_output$Exp <- NULL
data2_output$rep <- NULL
data2_output <- data2_output %>% spread(Sample, Value)
data2_output$Rt <- NULL
###data2_output <- data2_output[,c('Name', 'Used_ID', 'KEGG.ID', 'Pathway', 'Nr.C', 'Iso',info$Sample.Name)]
#data2_output <- data2_output[,c('Name', 'Used_ID', 'KEGG.ID', 'Vanq.Method', 'Nr.C', 'Iso',info$Sample.Name)]
data2_output <- data2_output[,c('Name', 'Used_ID', 'KEGG.ID', 'Pathway', 'Nr.C', 'Iso',info$Sample.Name)]



colnames(data2_output) <- c(colnames(data2_output)[1:6], info$Sample)
data2_output <- data2_output[order(data2_output$Name),]
num_samples <- sum(!grepl('QC',info$Sample))

for (i in 1:nrow(Anno))
  data2_output$Used_ID[as.character(data2_output$Name) == Anno$abbrev[i]] <- as.character(Anno$full_name[i])

data2_output <- rbind(c(rep(NA, 6), as.vector(info$Condition[1:num_samples]), rep(NA, (ncol(data_output) - 6 - num_samples))), data2_output)
combined_sample_condition<-c()
for (i in 7:(num_samples+6)){
  x<-paste(colnames(data2_output)[i],data2_output[1,i],sep=" - ")
  combined_sample_condition<-c(combined_sample_condition,x)
}

colnames(data2_output)[7:(num_samples+6)]<-combined_sample_condition
data2_output <- data2_output[-1,]
for(i in 7:ncol(data2_output))
{
  data2_output[,i] <- as.numeric(data2_output[,i])
  data2_output[,i] <- round(data2_output[,i])
}
setwd(output_dir)
write.xlsx(x = data2_output, file = output_file, row.names = F, col.names = T, showNA = F)








amounts <- data2 %>%
  dplyr::select(Name, KEGG.ID, Condition, Iso, Exp, Value) %>%
  group_by(Name, Condition, Exp) %>%
  mutate(Amount=sum(Value, na.rm=T)) %>%
  ungroup() %>% 
  filter(Iso=='M0')
amounts$Amount[amounts$Amount==0] <- NA     #need to convert 0 to NA for Av calculation

amounts_with_blanks <- data2_with_blanks %>%
  dplyr::select(Name, KEGG.ID, Condition, Iso, Exp, Value) %>%
  group_by(Name, Condition, Exp) %>%
  mutate(Amount=sum(Value, na.rm=T)) %>%
  ungroup() %>% 
  filter(Iso=='M0')
amounts_with_blanks$Amount[amounts_with_blanks$Amount==0] <- NA     #need to convert 0 to NA for Av calculation

amounts <- amounts %>%
  group_by(Name, Condition) %>%
  mutate(Av=mean(Amount, na.rm=T),
         Std=sd(Amount, na.rm=T),
         CV=sd(Amount, na.rm=T)/mean(Amount, na.rm=T)) %>%
  dplyr::select(Name, KEGG.ID, Condition, Exp, Amount, Av, Std, CV) %>%
  ungroup() %>%
  arrange(Condition, Name)

amounts_with_blanks <- amounts_with_blanks %>%
  group_by(Name, Condition) %>%
  mutate(Av=mean(Amount, na.rm=T),
         Std=sd(Amount, na.rm=T),
         CV=sd(Amount, na.rm=T)/mean(Amount, na.rm=T)) %>%
  dplyr::select(Name, KEGG.ID, Condition, Exp, Amount, Av, Std, CV) %>%
  ungroup() %>%
  arrange(Condition, Name)

test1 <- split(amounts, amounts[c('Name','Condition')])
NA_function <- function(x) if (sum(is.na(x)) < length(x)) return(x=x) else return(x=rep(0, length(x)))
new.Value <- as.vector(sapply(test1, function(x) NA_function(x$Amount)))
if (class(new.Value)=='list') new.Value <- unlist(new.Value)

test1_with_blanks <- split(amounts_with_blanks, amounts_with_blanks[c('Name','Condition')])
NA_function <- function(x) if (sum(is.na(x)) < length(x)) return(x=x) else return(x=rep(0, length(x)))
new.Value_with_blanks <- as.vector(sapply(test1_with_blanks, function(x) NA_function(x$Amount)))
if (class(new.Value_with_blanks)=='list') new.Value_with_blanks <- unlist(new.Value_with_blanks)

amounts <- amounts %>%
  mutate(Amount=new.Value)



amounts_with_blanks <- amounts_with_blanks %>%
  mutate(Amount=new.Value_with_blanks)

filtered_metabolites<-unique(amounts$Name)

amounts_with_blanks<-amounts_with_blanks[amounts_with_blanks$Name %in% filtered_metabolites,]


data8 <- split(amounts, amounts[,1])

data8_infinite <- vector("list", length = length(data8))
data8_finite <- vector("list", length = length(data8))
### Extract rows with NA or Inf Amount in loop
for (i in 1:length(data8)) {
  #not_finite <- data8[[1]][!is.finite(data8[[1]][,5]),]
  non_finite <- data8[[i]][!is.finite(rowSums(data8[[i]][,5])),]
  finite <- data8[[i]][is.finite(rowSums(data8[[i]][,5])),]
  #df <- df[is.finite(rowSums(df)),]
  data8_infinite[[i]] <- non_finite
  data8_finite[[i]] <- finite
}

### Run ANOvA on the rest
ANOVA <- suppressWarnings(sapply(data8_finite, function(x) anova(aov(x$Amount~x$Condition))$Pr[1]))
###ANOVA <- suppressWarnings(sapply(data8, function(x) anova(aov(x$Amount~x$Condition))$Pr[1]))
### Add the subsetted rows back
ANOVA <- rep(ANOVA, each=length(unique(info$Condition)))

data8_with_blanks <- split(amounts_with_blanks, amounts_with_blanks[,1])

data8_blanks_infinite <- vector("list", length = length(data8_with_blanks))
data8_blanks_finite <- vector("list", length = length(data8_with_blanks))
### Extract rows with NA or Inf Amount in loop
for (i in 1:length(data8_with_blanks)) {
  #not_finite <- data8[[1]][!is.finite(data8[[1]][,5]),]
  non_finite <- data8_with_blanks[[i]][!is.finite(rowSums(data8_with_blanks[[i]][,5])),]
  finite <- data8_with_blanks[[i]][is.finite(rowSums(data8_with_blanks[[i]][,5])),]
  #df <- df[is.finite(rowSums(df)),]
  data8_blanks_infinite[[i]] <- non_finite
  data8_blanks_finite[[i]] <- finite
}

###ANOVA_with_blanks <- suppressWarnings(sapply(data8_with_blanks, function(x) anova(aov(x$Amount~x$Condition))$Pr[1]))
ANOVA_with_blanks <- suppressWarnings(sapply(data8_blanks_finite, function(x) anova(aov(x$Amount~x$Condition))$Pr[1]))
ANOVA_with_blanks <- rep(ANOVA_with_blanks, each=length(unique(info_with_blanks$Condition)))

amounts <- amounts %>%
  arrange(Name) %>%
  spread(Exp, Amount) %>%
  cbind('ANOVA'=ANOVA, 'Sig'=NA)

amounts_with_blanks <- amounts_with_blanks %>%
  arrange(Name) %>%
  spread(Exp, Amount) %>%
  cbind('ANOVA'=ANOVA_with_blanks, 'Sig'=NA)


for (i in 1:nrow(amounts)){
  if (amounts$ANOVA[i] == "NaN") amounts$Sig[i]=""
  else if (amounts$ANOVA[i] <= 0.001) amounts$Sig[i]="***"
  else if (amounts$ANOVA[i] <= 0.01) amounts$Sig[i]="**"
  else if (amounts$ANOVA[i] <= 0.05) amounts$Sig[i]="*"
  else amounts$Sig[i]=""
}
write.csv(amounts, paste0(Title,"-Amounts",'.csv'), row.names=FALSE)
save(amounts, file='amounts.rdata')
amounts<-check_50(amounts)
amounts_with_blanks<-check_50(amounts_with_blanks)

plotname=paste(Title,"-plots-plus fresh media plus blanks1.pdf", sep='')
pdf(file = plotname, width=14, height=10, pointsize=12)

amounts %>%
  mutate(Name = paste(Name, sep = ' ')) %>%
  ggplot(., aes(Condition, Av, group=Condition, fill=Condition)) + 
  geom_bar(aes(linetype=under_50_percent,color = under_50_percent, size=under_50_percent),position="dodge", stat="identity", width=0.9) +
  guides(linetype=FALSE)+ scale_linetype_manual(values=c("solid","58")) + scale_size_manual(values=c(0.3,0.8), guide = F) + scale_colour_manual(values = c("black","gray29"), guide = F) +
  facet_wrap( ~ Name, scales="free", nrow=floor(sqrt(length(unique(amounts$Name))))) +
  theme_bw() +
  labs(x="", y="Relative amounts", title=Title, fill=element_blank()) +
  theme(
    plot.title=element_text(size=20, face="bold", vjust=2),         #sets title properties
    axis.title=element_text(size=16, lineheight=20, face="bold"),   #sets theme for axis font
    axis.text=element_text(size=11, face="bold"),
    axis.text.x=element_blank(),
    legend.title=element_text(face="bold", size=12),
    legend.text=element_text(face="bold",size=12),                  #sets legend text
    strip.text=element_text(face="bold", size=15),           #sets theme for title in facets
    panel.grid.major=element_blank()) +
  scale_fill_manual("Conditions", values =colors1)  +
  geom_errorbar(aes(ymin=Av, ymax=Av+Std), position=position_dodge(0.9), width=.2)
dev.off()



num_non_blanks<-length(which(!grepl("blank",unique(amounts_with_blanks$Condition))))

colors1[num_non_blanks+1:length(unique(amounts_with_blanks$Condition))]<-"grey45"

plotname=paste(Title,"-plots-plus fresh media plus blanks2.pdf", sep='')
pdf(file = plotname, width=14, height=10, pointsize=12)

amounts_with_blanks %>%
  mutate(Name = paste(Name, sep = ' ')) %>%
  ggplot(., aes(Condition, Av, group=Condition, fill=Condition)) + 
  geom_bar(aes(linetype=under_50_percent,color = under_50_percent, size=under_50_percent),position="dodge", stat="identity", width=0.9) +
  guides(linetype=FALSE)+ scale_linetype_manual(values=c("solid","58")) + scale_size_manual(values=c(0.3,0.8), guide = F) + scale_colour_manual(values = c("black","gray29"), guide = F) +
  facet_wrap( ~ Name, scales="free", nrow=floor(sqrt(length(unique(amounts$Name))))) +
  theme_bw() +
  labs(x="", y="Relative amounts", title=Title, fill=element_blank()) +
  theme(
    plot.title=element_text(size=20, face="bold", vjust=2),         #sets title properties
    axis.title=element_text(size=16, lineheight=20, face="bold"),   #sets theme for axis font
    axis.text=element_text(size=11, face="bold"),
    axis.text.x=element_blank(),
    legend.title=element_text(face="bold", size=12),
    legend.text=element_text(face="bold",size=12),                  #sets legend text
    strip.text=element_text(face="bold", size=15),           #sets theme for title in facets
    panel.grid.major=element_blank()) +
  scale_fill_manual("Conditions", values =colors1)  +
  geom_errorbar(aes(ymin=Av, ymax=Av+Std), position=position_dodge(0.9), width=.2)
dev.off()


##May need to adjust the colors used inorder to match the number of conditions, sometimes the condition number is too great
##Consider changing the pallete by changin 'scale_fill_manual' to 'scale_fill_brewer("Conditions", palette= "Accent")'
amounts2 <- amounts %>%
  dplyr::select(Name, KEGG.ID, Condition, contains('Exp')) %>%
  gather(Exp, Amount, -Name, -KEGG.ID, -Condition) %>%
  arrange(Name, Condition) %>%
  unite(Condition_Exp, c(Condition,Exp), sep='_') %>%
  mutate(Condition_Exp=factor(Condition_Exp, levels=unique(Condition_Exp))) %>%
  #  mutate(Condition_Exp = Condition:Exp, Condition=NULL, Exp=NULL) %>%
  spread(Condition_Exp, Amount) %>%
  dplyr::select_if(~sum(!is.na(.))>0) %>%
  ungroup()


fresh_regex = '[Uu]nspent|[Bb]lank|[Ff]resh|[Cc]ontrol'
#fresh_regex = '[Uu]nspent|[Bb]lank|[Ff]resh|[Cc]ontrol|[Mm]edium|[Mm]edia'
if (0) {
  #special code for data set '2020.07.07_Jenna-Godstein_16D Enza_metabolic footprinting'
  fresh_regex = '[Uu]nspent|[Bb]lank|[Ff]resh|[Cc]ontrol|Media complete|Media -SG'
}

fresh <- info[grep(fresh_regex, info$Condition),]$Samples  #find blank medium samples
fresh2 <- info[grep(fresh_regex, info$Condition),]

if (length(which(sapply(amounts2[,3:length(amounts2)], sum, na.rm=T)==0)) > 0){
  amounts2 <- amounts2[-(which(sapply(amounts2[,3:length(amounts2)], sum, na.rm=T)==0)+2)]   #remove columns with only NAs
}  

#info$Medium[1:9] <- "A"

# the structure of the experimnental design and of the sample info file is expected to be such that 
# the unspent media controls that should be subtracted in the footprinting plots come in sets of 3 replicates, 
# and that the first one will be used for the samples with A in the Medium column, second B, etc.
# samples that do not need a subtraction, can be assigned an additional letter, and put into the unspent dataframne
# (if other designs are used, this code will need adjustment)
for (i in seq_len(length(unique(na.omit(info$Medium))))) {                                 #create variable with average amount in blank medium
  #print(i)
  assign(LETTERS[i],
         apply(amounts2[,fresh[c(i*3-2, i*3-1, i*3)]+2], 1, mean, na.rm=T))
  #confirm that the Medium coilumn in the sample info file is structured correctly (see above)
  test <- fresh2[grep(LETTERS[i], fresh2$Medium),]$Medium
  if (length(test) !=3) {
    stop("need to evolve the code to deal with unspent media that is not in replicates of 3")
  }
  if (length(unique(test)) != 1 | test[1] != LETTERS[i]) {
    stop("the Medium column in the sample info file does not have the expected structure")
  }
}

#change NAs in 'fresh' variables to 0
####samples <- info[-grep('unspent|[Bb]lank|[Ff]resh|[Cc]ontrol|[Mm]edium|[Mm]edia', info$Condition),]
samples <- info[-grep(fresh_regex, info$Condition),]
amounts3 <- amounts2



for (i in seq_len(length(unique(samples$Medium)))){
  print(i)
  spent <- samples[grep(LETTERS[i], samples$Medium),]$Samples + 2
  amounts3[,spent] <- amounts3[,spent]-get(LETTERS[i])
}
amounts3 <- amounts3[,-(fresh+2)]
amounts3[is.na(amounts3)] <- 0


condition = as.character(samples$Sample.Name)
### If Cell.Number is blank, run the line below
samples$Cell.Number <- 1
numbers=samples$Cell.Number
Norv=Norv[-fresh]


numbers<-as.numeric(numbers)
for (i in 1: length(numbers)){
  for (j in 1: nrow(amounts3)) amounts3[j,i+2]=amounts3[j,i+2]/numbers[i]
}
#if above gives error, make sure you have numerics in Cell.numer column
# write.csv(amounts3,paste0(Title,'-Normalized Amounts.csv'), row.names=F)
# save(amounts3, file='Normalized Amounts.rdata')

#making a heatmap
amounts3 <- amounts3[,!(colSums(amounts3[3:length(amounts3)], na.rm=T)==0)]
amounts3 <- amounts3[!(rowSums(amounts3[3:length(amounts3)], na.rm=T)==0),]
data5=as.matrix(amounts3[3:length(amounts3)])
rownames(data5)=amounts3[,1]
colnames(data5)=names(amounts3[3:length(amounts3)])

#fill with arbitrary values if numbers is N/A
for (i in 1:length(numbers))
{
  if(numbers[i] == "N/A")
    numbers[i] <- 1
}
numbers <- as.numeric(numbers)

Norv<-rep(1,length(numbers))
ann <- data.frame("Condition"=as.character(samples$Condition),'Cell Number'=numbers, 'Norvaline'=Norv)
Sig.color=c("blue", "red")
names(Sig.color)=unique(samples$Type.1)
Norv.color = c("white", "blue")
Cell.Number.color= c("white", "green")
names(colors1)=unique(samples$Condition)
my_cell_type <- c("#7FC97F", "#BEAED4", "#FDC086" ) 
###^^ make this a factor and make sure to add it into the ann_colors list below
ann_colors = list(Condition=colors1[1:length(unique(samples$Condition))], Norvaline = Norv.color, Cell.Number = Cell.Number.color, Cell.Type = my_cell_type )
#ann_colors = list(Condition=colors[1:length(unique(samples$Condition))], Norvaline = Norv.color, Cell.Number = Cell.Number.color)
#up for discussion whether or not to include all the annotation legends
rownames(ann)=colnames(data5)

write.csv(data5, paste0(Title,'-heatmap data.csv'), row.names=T)
data5[is.na(data5)] <- 0     #replace NAs
heatmap.title=paste(Title, "-All-Clustered Heatmap.pdf", sep='')
#dev.off()

pheatmap(data5, cluster_row=T, cluster_col=T, clustering_distance_rows='correlation', clustering_distance_cols='correlation',
         color = colorRampPalette(c("blue", "white", "red"))(100), border_color="black", scale="row",
         cellwidth = 30, cellheight = 18, annotation=ann, annotation_colors = ann_colors, show_colnames = F,
         main=Title, filename=heatmap.title)
dev.off()
dev.set(dev.next())
#make unclustered heatmap too if needed

heatmap.title=paste(Title, "-Unclustered Heatmap.pdf", sep='')
#dev.off()

pheatmap(data5, cluster_row=T, cluster_col=F, clustering_distance_rows='correlation', clustering_distance_cols='correlation',
         color = colorRampPalette(c("blue", "white", "red"))(100), border_color="black", scale="row",
         cellwidth = 30, cellheight = 18, annotation=ann, annotation_colors = ann_colors, show_colnames = F,
         main=Title, filename=heatmap.title)
dev.off()
dev.set(dev.next())

#If annotations legend becomes longer than the metabolite list, change size of the font i.e. fontsize=6
amounts3[,3:length(amounts3)][amounts3[,3:length(amounts3)]==0] <- NA



##### Considers the NA value that could be present in one of the replicates per condition and fixes 'NA' error in Average Calcualtion

amounts3 <- amounts3 %>%
  gather(Condition_Exp, Amount, -Name, -KEGG.ID) %>%
  separate(Condition_Exp, c( 'Condition', 'Exp'), sep='_') %>%
  mutate(Condition=factor(Condition, levels=unique(.$Condition)))
amounts3$Amount <- suppressMessages(mapvalues(amounts3$Amount, c('Inf','NaN'), c(NA,NA)))
amounts3$Amount[amounts3$Amount==0] <- NA
test1 <- split(amounts3, amounts3[c('Name','Condition')])
NA_function <- function(x) if (sum(is.na(x)) < length(x)) return(x=x) else return(x=rep(0, length(x)))
new.Value <- as.vector(sapply(test1, function(x) NA_function(x$Amount)))

#only in case new.Value is not a vector of the same length as amounts3 and is instead a list of shorter length, use this function
#####check to make sure amounts3 row numbers and new.Value are of the same length
new.Value <- unlist(new.Value, use.names = F)
###########################################################################


amounts3 <- suppressWarnings(amounts3 %>%
                               arrange(Condition, Name) %>%
                               mutate(Amount=new.Value) %>%
                               group_by(Name, Condition) %>%
                               mutate(Av=mean(Amount, na.rm=T),
                                      Std=sd(Amount, na.rm=T),
                                      CV=sd(Amount, na.rm=T)/mean(Amount, na.rm=T)) %>%
                               ungroup())

#####

#amounts3 <- amounts3 %>%
#  gather(Condition_Exp, Amount, -Name, -KEGG.ID) %>%
#  separate(Condition_Exp, c('Condition', 'Exp'), sep='_') %>%
#  mutate(Condition=factor(Condition, levels=unique(.$Condition))) %>%
#  group_by(Name, Condition) %>%
#  mutate(Av=mean(Amount),
#         Std=sd(Amount),
#         CV=Std/Av) %>%
#  ungroup()

amounts3$Amount <- mapvalues(amounts3$Amount, NA, 0)

data8=split(amounts3, amounts3[,1])
ANOVA=suppressWarnings(sapply(data8, function(x) anova(aov(x$Amount~x$Condition))$Pr[1]))
ANOVA=rep(ANOVA,1,each=length(unique(amounts3$Condition)))

amounts3 <- amounts3 %>%
  spread(Exp, Amount) %>%
  arrange(Name) %>%
  mutate(ANOVA = ANOVA,
         Sig = '')

for (i in 1:nrow(amounts3)){
  if (amounts3$ANOVA[i] == "NaN") amounts3$Sig[i]=""
  else if (amounts3$ANOVA[i] <= 0.001) amounts3$Sig[i]="***"
  else if (amounts3$ANOVA[i] <= 0.01) amounts3$Sig[i]="**"
  else if (amounts3$ANOVA[i] <= 0.05) amounts3$Sig[i]="*"
  else amounts3$Sig[i]=""
}
write.csv(amounts3, paste0(Title,'-Amounts normalized to unspent.csv'), row.names=FALSE)
save(amounts3, file='amounts3.rdata')
amounts3[is.na(amounts3)] <- 0

for (i in 1:nrow(amounts3)) {
  if (amounts3$Av[i] < 0) amounts3$Std[i] = amounts3$Std[i] * -1
}
amounts3<-check_50(amounts3)
plotname=paste(Title,"-plots-relative-to-unspent1.pdf", sep='')
pdf(file = plotname, width=14, height=10, pointsize=12)

amounts3 %>%
  mutate(Name = paste(Name, Sig, sep=' ')) %>%
  ggplot(., aes(Condition, Av, group=Condition, fill=Condition)) +
  geom_bar(aes(linetype=under_50_percent,color = under_50_percent, size=under_50_percent),position="dodge", stat="identity", width=0.9) +
  guides(linetype=FALSE)+ scale_linetype_manual(values=c("solid","58")) + scale_size_manual(values=c(0.3,0.8), guide = F) + scale_colour_manual(values = c("black","gray29"), guide = F) +
  facet_wrap( ~ Name, scales="free", nrow=floor(sqrt(length(unique(amounts3$Name))))) +
  theme_bw() +
  labs(x="", y="Relative amounts", title=Title, fill=element_blank()) +
  theme(
    plot.title=element_text(size=20, face="bold", vjust=2),         #sets title properties
    axis.title=element_text(size=16, lineheight=20, face="bold"),   #sets theme for axis font
    axis.text=element_text(size=11, face="bold"),
    axis.text.x=element_blank(),
    legend.title=element_text(face="bold", size=12),
    legend.text=element_text(face="bold",size=12),                  #sets legend text
    strip.text=element_text(face="bold", size=15),           #sets theme for title in facets
    panel.grid.major=element_blank()) +
  scale_fill_manual("Conditions", values = colors1)  +
  geom_errorbar(aes(ymin=Av, ymax=Av+Std), position=position_dodge(0.9), width=.2)
dev.off()
save.image(paste0(Title,'.rdata'))

#calculating isotopomer distribution
samples=info
samples$Sample=1:nrow(samples)
condition=samples$Sample.Name
condition = as.character(samples$Sample.Name)
samples$Condition=factor(samples$Condition, levels=unique(samples$Condition))
numbers=samples$Cell.Number
### Changed from colors to colors1
rm(colors1)


select <- dplyr::select
##Maybe incorporate this function at the beginning of this script to avoid 'dplyr::select' notation
data3 <- make_MID(data2)
names(data3)[names(data3) == "Norm_Av"] <- "RelAmounts_Ave"
names(data3)[names(data3) == "Norm_Std"] <- "RelAmounts_Std"
if (unlabelled == FALSE){
  pdf(paste0(Title,'-MID plots.pdf'), width=14, height=10)
  bar_update_manual(unique(data3$Name), data3,n=num_conditions,type="medium")
  dev.off()
  save.image(paste0(Title,'.rdata'))}

