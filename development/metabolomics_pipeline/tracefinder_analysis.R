# devtools::install_github("danielbraas/MetabFUN")
# devtools::install_github("juyeki/MetabR")
library(MetabFUN)
library(MetabR)
library(tidyr)
library(readxl)
library(dplyr)
library(ggrepel)
#Sys.setenv(JAVA_HOME = "C:\\Program Files\\Java\\jdk-11.0.8\\")
library(xlsx)

###### parameters ######

#data_dir <- "M:/TraceFinderData/4.0/Projects/Cells/Nakamura Lab/JN-08202020_7-12_Gln_Vanq" #Klara
#output_dir <- "M:/TraceFinderData/4.0/Projects/Cells/Nakamura Lab/JN-08202020_7-12_Gln_Vanq/test2" #Klara

#mac paths
data_dir <- "/Users/tgraeber/Dropbox/glab/Metabolomics Center/Metabolomics Pipeline Scripts and Notes/workspace/Projects/Cells/Nakamura JN-08202020_7-12_Gln_Vanq"
output_dir <- "/Users/tgraeber/Dropbox/glab/Metabolomics Center/Metabolomics Pipeline Scripts and Notes/workspace/Projects/Cells/Nakamura JN-08202020_7-12_Gln_Vanq/test3"

#set tracer_type to: "full", "partial", or "none"
tracer_type <- "full"
diff_computers <- FALSE
ifelse (grepl("ICS", data_dir), ICS <- TRUE, ICS <- FALSE)

#set natural isotope abundance parameters
if (tracer_type == "full" || tracer_type == "partial") {
  #turn on one of these options
  if (1) {
    isotope_correction_flag = "default" #for using the 1996 Fernandez et al. JofMS method (first method used by the UCLA Metabolomics Center)
  } else if (0) {
    isotope_correction_flag = "IsoCorrectoR" #for using the IsoCorrectoR method, via our correct_iso function
    label <- "C"
    correct_for <- "C"
  } else if (0) {
    isotope_correction_flag = "force_uncorrected"
  } else {
    stop("need to set isotope_correction_flag")
  }
}


#Abbreviation Data
library(googlesheets4)
#Abbrev_NEW
#Abbrev <- read_sheet("https://docs.google.com/spreadsheets/d/118M3rvfJAOQrOoYEHZVjEFg_k0CMwK9i8wq56u_ZXP4/edit?ts=5eaa1a9e#gid=220628921")
#Abbrev_NEW2
#Abbrev <- read_sheet("https://docs.google.com/spreadsheets/d/1He49QJYE1ld0VUgzNTNht-AuoznydEL9ysjCPRkgk1Y/edit#gid=220628921")

### Read in the most recent version of Abbrev
#Abbrev <-  read_excel("C:/Users/FTsang/Downloads/Abbrev_NEW2 (1).xlsx")
Abbrev <-  read_excel("/Users/tgraeber/Dropbox/glab/Metabolomics Center/Metabolomics Pipeline Scripts and Notes/workspace/Abbrev_NEW2.xlsx")
Abbrev <- Abbrev[-c(248,249),]

#Title
setwd(data_dir)
#Title <- gsub('.xls[x]?','', list.files(pattern='.xls[x]?'))
Title <- gsub('.xls[x]?','', list.files(pattern='[Ss]ample.*.xls[x]?'))
Title <- Title
Title <- gsub('_sample info sheet|_sample info|sample form|_sampleinfosheet|_Sampleinfosheet', '', Title)
Title <- gsub(' ', '', Title)

#info 
#info <- read_excel(list.files()[grep('.xls[x]?',list.files())])
info <- read_excel(list.files()[grep('[Ss]ample.*.xls[x]?',list.files())])
#info$Condition <- gsub("^\\s+|\\s+$", "", info$Condition)     
info$Condition <- gsub('/', '-', info$Condition)
info$Condition <- gsub('_', '-', info$Condition)
info$Sample <- gsub('-', '.', info$Sample)
info$Condition <- factor(info$Condition, levels = unique(info$Condition))   

samples <- info$Sample
samples <- samples[!grepl(paste0("QC", collapse = "|"), samples)]
num_samples<-length(samples)

#data
setwd(data_dir)
folder.name <- gsub('(.)*Projects/', '', data_dir) %>% gsub('/', '-', .)
if (ICS)
{
  files <- "ICS/ICS.csv"
} else
{
  files <- c("AA/AA.csv", "CoAs/CoAs.csv", "Currency/Currency.csv", "dNTPs/dNTPs.csv", "FA/FA.csv", "Glycolysis/Glycolysis.csv", "Hexosamine/Hexosamine.csv", "NTPs/NTPs.csv", "PPP/PPP.csv", "TCA/TCA.csv","Urea/Urea.csv")
  # files <- c("AA/AA-2.csv", "CoAs/CoAs-2.csv", "Currency/Currency-2.csv", "dNTPs/dNTPs-2.csv", "FA/FA-2.csv", "Glycolysis/Glycolysis-2.csv", "Hexosamine/Hexosamine-2.csv",
  #            "NTPs/NTPs-2.csv", "PPP/PPP-2.csv", "TCA/TCA-2.csv","Urea/Urea-2.csv")
}
if (sum(grepl('ICS',files)) > 0) Title <- paste0('ICS-', Title)
data <- lapply(files, read.csv, header=T, na.strings = c("N/F",'N/A'))

if (!diff_computers)
{
  data <- do.call(rbind, data)
} else
{
  library(plyr)
  data <- do.call(rbind.fill, data)
  detach(package:plyr, unload=TRUE)
}
data$Experiment <- folder.name
data[,1] <- NULL
setwd(output_dir)
write.csv(data, 'all data raw.csv', row.names=F)

data <- data %>%
  dplyr::select(Compound, Filename, Area) %>%
  spread(Filename, Area)

#info$Sample <- colnames(data)[2:17]
#info$Sample <- colnames(data)[-1]
keep <- gsub("\\.", "-", info$Sample)
keep <- c("Compound", keep)
data <- data[, which(colnames(data) %in%  keep)]
### Move QC columns behind actual samples
### data <- data[,c(1, which(grepl("SL-", colnames(data))), which(grepl("QC-", colnames(data))))]
data_temp <- filter_data(data, info)

data$Compound = as.character(data$Compound)
for (i in 1:length(data$Compound)){                   #loop through IDs and fill up with Ms
  if (grepl(" M[0-9]+$", data$Compound[i])==FALSE) data$Compound[i]=paste(data$Compound[i],"M0")}
Used_ID <- gsub(" M[0-9]+$| Std", "", data$Compound)
name_split <- regexpr("M[0-9]+$|Std",data$Compound)
Iso <-regmatches(data$Compound, name_split)
data <- data.frame(Used_ID, Iso, data[,2:length(data)])
data <- data %>% select(-grep('QC_|QC.', colnames(data)), grep('QC_|QC.', colnames(data)))

data_output <- left_join(x = data, y = Abbrev, by = "Used_ID")
data_output <- data_output %>% dplyr::rename(Name = Abb)

#data_output <- data_output %>% select(Name, Used_ID, KEGG.ID, Pathway, Nr.C, Iso, info$Sample)
info$Sample <- gsub("-",".",info$Sample)
data_output <- data_output %>% select(Name, Used_ID, KEGG.ID, Vanq.Method, Nr.C, Iso, info$Sample)

if(tracer_type == "none")
  data_output <- subset(data_output, Iso == "M0" | Iso =="M1")

data_output <- data_output[!is.na(data_output$Nr.C),]
isotopes <- rep(NA, max(data_output$Nr.C))
for (i in 1:(length(isotopes)+1))
  isotopes[i] <- paste0("M", as.character(i-1))
data_output <- data_output %>% 
  mutate(Iso = factor(Iso, levels = isotopes)) %>%
  arrange(Name, Iso)
num_samples <- sum(!grepl('QC',info$Sample))

data(anno)
#data(anno, package="MetabR")
# data(package = "MetabR") #lists what is available
#load(file = "/Users/tgraeber/Dropbox/glab/Metabolomics Center/Metabolomics Pipeline Scripts and Notes/pipeline scripts/githup sets - metabolomics R pipeline/MetabR-master/data/anno.rda")

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

info <- info[!grepl("QC.blank", info$Sample), ] 
data <- data_temp
data <- data[ , -which(grepl('QC.blank', colnames(data)))]
data$Compound = as.character(data$Compound)
for (i in 1:length(data$Compound)){                   #loop through IDs and fill up with Ms
  if (grepl(" M[0-9]+$", data$Compound[i])==FALSE) data$Compound[i]=paste(data$Compound[i],"M0")}
Used_ID <- gsub(" M[0-9]+$| Std", "", data$Compound)
name_split <- regexpr("M[0-9]+$|Std",data$Compound)
Iso <-regmatches(data$Compound, name_split)
data <- data.frame(Used_ID, Iso, data[,2:length(data)])
data <- data %>% select(-grep('QC_|QC.', colnames(data)), grep('QC_|QC.', colnames(data)))

### Don't plot Norvaline in the ISTD plots
std <- make_istd(data, remove = FALSE)$Std
std <- std[std$Used_ID != "Norvaline",]
setwd(output_dir)
info_ordered <- info[order(info$Run.Order),]
if(ICS)
{
  istd_title <- paste0("QC-ISTDs with 50ks-", Title, ".pdf")
} else
{
  istd_title <- paste0("QC-ISTDs with 250ks-", Title, ".pdf")
}
plot_istd(Std = std, info = info_ordered, title = istd_title)

#ISTD plot with only samples
data <- data[ , !grepl('QC.', colnames(data))]
info <- info[!grepl("QC.", info$Sample), ]
std <- make_istd(data, remove = TRUE)$Std
std <- std[std$Used_ID != "Norvaline",]
data <- make_istd(data, remove = TRUE)$data 
setwd(output_dir)
info_ordered <- info[order(info$Run.Order),]
plot_istd(Std = std, info = info_ordered, title = paste0("QC-ISTDs-", Title, ".pdf"))

#removing standards
data <- filter(data, !(grepl('Std', data$Iso) & data$Used_ID != "Norvaline"))
#data <- filter(data, !(grepl('Std', data$Iso)))
data$Used_ID <- as.character(data$Used_ID)
data$Iso <- factor(data$Iso, levels=paste0('M',0:50))
data <- plyr::arrange(data, Used_ID, Iso)

#getting sample conditions
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

#adding conditions to data
new_colnames <- info$Sample
non_blanks_names <- new_colnames[!grepl(paste0("blank",collapse="|"),new_colnames)]
blanks_names <- new_colnames[grepl(paste0("blank",collapse="|"),new_colnames)]
col.order <- c("Used_ID","Iso",non_blanks_names,blanks_names)
col.order <- gsub("-",".",col.order)
data <- data[,col.order]
for (i in 1:length(info$Sample.Name)) 
  colnames(data)[i+2] <- info$Sample.Name[i]
info_ordered <- info[order(info$Run.Order),]

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

#Normalizing Norvaline 
if (length(Norv)==0){
  Norv <- rep(1, length(samples))
  #Norv <- rep(1,nrow(samples))
  print("No norvaline")
}
norm <- data
norm <- norm[!norm$Used_ID=="Norvaline",]    

####
to_normalize <- !(all(is.na(info$Cell.Number)))
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
  dplyr::mutate(NA_list=sum(is.na(Value)),
         NA_potential=n()) %>%
  filter(NA_list==NA_potential) %>%
  distinct()
data2 <- data2 %>%
  filter(!(Name %in% NA_Names$Name))
data2$Name <- as.character(data2$Name)
data2 <- data2[!is.na(data2$Name),]

if(tracer_type == "none")
  data2 <- subset(data2, Iso == "M0" | Iso =="M1")

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
data2_output <- data2_output[,c('Name', 'Used_ID', 'KEGG.ID', 'Vanq.Method', 'Nr.C', 'Iso',info$Sample.Name)]
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

temp <- (data %>% group_by(Used_ID) %>% filter(n() == 1))
temp <- left_join(temp, Abbrev, by = 'Used_ID')
only_M0 <- temp$Abb

select <- dplyr::select
data4 <- RelAmounts(data2, anova = TRUE, output = T)
### Try changing NaN and Inf values in data4 columns CV, Norm_Av, and Norm_Std to NA
#data4$CV[!is.finite(data4$CV)] <- NA
#data4$Norm_Av[!is.finite(data4$Norm_Av)] <- NA
#data4$Norm_Std[!is.finite(data4$Norm_Std)] <- NA
# Just close R and rerun if error appears regarding "n()" - seems to have to do with the order packages are loaded into the environment
RelA <- make_matrix(data4)

colors <- c("turquoise","red","plum4","steelblue1","red4","springgreen2","slateblue2","sienna1","darkgreen","lightpink1","navy","olivedrab1",
            "orangered","darkslateblue","lightseagreen","magenta2","royalblue","yellowgreen","lightsalmon","cyan","maroon1","indianred3","mediumseagreen",
            "slateblue3","hotpink","lemonchiffon1","orangered4","lightcoral","tomato")

###### make MID, FC and data_labeled dataframes and execute natureal abundance correction ######
if(tracer_type == "full" || tracer_type == "partial")
{
  
  print("making MID files, and performing natural isotope correction (",isotope_correction_flag,")")
  if (isotope_correction_flag == "default") {
  data3 <- make_MID(data2) #make_Mid also executes the natural isotope abundance corrections using an interative algorithm
  #data3 <- make_MID_tg(data2) #make_Mid also executes the natural isotope abundance corrections using an interative algorithm
  
  #tfdata2 <- data2
  #tfdata3 <- data3
  #save(tfdata2, file = "tracefinder_data2.rda")
  #save(tfdata3, file = "tracefinder_data3.rda")
  # load("tracefinder_data2.rda")
  # load("tracefinder_data3.rda")
  
  } else if (isotope_correction_flag == "IsoCorrectoR" || isotope_correction_flag == "force_uncorrected") {
    #uncorrected & correction using correct_iso / IsoCorrectoR
    
    if (isotope_correction_flag == "force_uncorrected") {
      dir.name = paste0(output_dir,"/uncorrected")
      if (!dir.exists(dir.name)) {dir.create(dir.name)}
      setwd(dir.name)
      # setwd(output_dir)
    }
    
    data3_uncorr <- make_MID_uncorrected(data2)
    #save(data3_uncorr, file = "tracefinder_data3_uncorr.rda")
    
    data3_corr <- data3
    data3 <- data3_uncorr
    
    if (isotope_correction_flag == "IsoCorrectoR") {
      
      dir.name = paste0(output_dir,"/corrected_IsoCorrectoR")
      if (!dir.exists(dir.name)) {dir.create(dir.name)}
      setwd(dir.name)
      # setwd(output_dir)

      data2_isocorrector <- correct_iso2_tg(data3_uncorr, label, correct_for)
      colnames(data2_isocorrector)[colnames(data2_isocorrector) == "Corrected_Value"] <- "Value"
      data3_isocorrector <- make_MID_uncorrected(data2_isocorrector)
      
      #save(data2_isocorrector, file = "tracefinder_data2_isocorrector.rda")
      #save(data3_isocorrector, file = "tracefinder_data3_isocorrector.rda")

      data3 <- data3_isocorrector
      
    }
  }
  
  MID <- make_matrix(data3)
  if(tracer_type == "full")
  {
    FC <- make_FC(data3)
    FC_mat <- make_matrix((FC))
  }
  if(tracer_type == "partial")
  {
    data_labeled <- make_data_labeled(data3)
    data_labeled_mat <- make_matrix(data_labeled)
  }
}




###### output plots ######

print("making output plots")

RelA <- RelA[which(rownames(RelA) != "ADPfromATP"),]
make_heatmap(RelA, info)
dev.set(dev.next()) 
make_heatmap(RelA, info, cluster_samples = F)
dev.set(dev.next())

MID <- MID[which(rownames(MID) != "ADPfromATP"), ]
data_labeled_mat <- data_labeled_mat[which(rownames(data_labeled_mat) != "ADPfromATP"), ]
FC_mat <- FC_mat[which(rownames(FC_mat) != "ADPfromATP"), ]

if(tracer_type == "full" || tracer_type == "partial")
{ 
  dev.set(dev.next()) 
  make_heatmap(MID, info)
  if(tracer_type == "full")
    make_heatmap(FC_mat, info)
  if(tracer_type == "partial")
    make_heatmap(data_labeled_mat, info)
  dev.set(dev.next())
}
dev.set(dev.next()) 

if(tracer_type == "full")
{
  make_heatmap(FC_mat, info, cluster_samples = F)
  dev.set(dev.next()) 
}

if(tracer_type == "partial")
{
  make_heatmap(data_labeled_mat, info, cluster_samples = F)
  dev.set(dev.next())
}
dev.set(dev.next())

samples <- info

make_PCA2_update(RelA)
if(tracer_type == "full" || tracer_type == "partial")
{
  make_PCA2_update(MID)
  if(tracer_type == "full")
    make_PCA2_update(FC_mat)
  if(tracer_type == "partial")
    make_PCA2_update(data_labeled_mat)
}

names(data4)[names(data4) == "Norm_Av"] <- "RelAmounts_Ave"
names(data4)[names(data4) == "Norm_Std"] <- "RelAmounts_Std"
names(data3)[names(data3) == "Norm_Av"] <- "RelAmounts_Ave"
names(data3)[names(data3) == "Norm_Std"] <- "RelAmounts_Std"
names(FC)[names(FC) == "Norm_Av"] <- "RelAmounts_Ave"
names(FC)[names(FC) == "Norm_Std"] <- "RelAmounts_Std"
names(data_labeled)[names(data_labeled) == "Norm_Av"] <- "RelAmounts_Ave"
names(data_labeled)[names(data_labeled) == "Norm_Std"] <- "RelAmounts_Std"

#nonpathway metabolites
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

data4 <- check_50(data4)
data3 <- check_50(data3)
data_labeled <- check_50(data_labeled)
FC <- check_50(FC)
for (metab in unique(data4$Name))
{
  if(all(is.na(data4[which(data4$Name == metab), 'RelAmounts_Ave'])))
    data4 <- data4[-which(data4$Name == metab), ]
}

##making graphs
plotname <- paste0(Title, "-Plots All.pdf")
if(tracer_type == "none")
  plotname <- paste0(Title, "-Plots RelAmounts.pdf")
pdf(file = plotname, width=14, height=10, pointsize=12)

bar_update_manual('glycolysis',data4, n = num_conditions, type = "tf")
bar_update_manual('glycolysis',data3, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual('glycolysis',data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual('glycolysis',FC, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("TCA",data4, n = num_conditions, type = "tf")
bar_update_manual("TCA",data3, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("TCA",data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("TCA",FC, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("PPP",data4, n = num_conditions, type = "tf")
bar_update_manual("PPP",data3, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("PPP",data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("PPP",FC, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("Curr",data4, n = num_conditions, type = "tf")
bar_update_manual("Curr",data3, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("Curr",data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("Curr",FC, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("Cys",data4, n = num_conditions, type = "tf")
bar_update_manual("Cys",data3, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("Cys",data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("Cys",FC, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("AA",data4, n = num_conditions, type = "tf")
bar_update_manual("AA",data3, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("AA",data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("AA",FC, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("FA",data4, n = num_conditions, type = "tf")
bar_update_manual("FA",data3, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("FA",data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("FA",FC, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("Hex",data4, n = num_conditions, type = "tf")
bar_update_manual("Hex",data3, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("Hex",data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("Hex",FC, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("Adenine",data4, n = num_conditions, type = "tf")
bar_update_manual("Adenine",data3, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("Adenine",data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("Adenine",FC, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("Cytosine",data4, n = num_conditions, type = "tf")
bar_update_manual("Cytosine",data3, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("Cytosine",data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("Cytosine",FC, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("Guanine",data4, n = num_conditions, type = "tf")
bar_update_manual("Guanine",data3, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("Guanine",data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("Guanine",FC, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("Thymine",data4, n = num_conditions, type = "tf")
bar_update_manual("Thymine",data3, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("Thymine",data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("Thymine",FC, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("Uracil",data4, n = num_conditions, type = "tf")
bar_update_manual("Uracil",data3, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("Uracil",data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("Uracil",FC, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual('Fru',data4, n = num_conditions, type = "tf")
bar_update_manual('Fru',data3, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("Fru",data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual('Fru',FC, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual('CoAs',data4, n = num_conditions, type = "tf")
bar_update_manual('CoAs',data3, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual("CoAs",data_labeled, n = num_conditions, type = "tf", only_M0 = only_M0)
bar_update_manual('CoAs',FC, n = num_conditions, type = "tf", only_M0 = only_M0)
#first set
bar_update_manual(unname(unlist(non_data4[1])),data4, n = num_conditions, type = "tf", title_type = "nonpathway1")
bar_update_manual(unname(unlist(non_data3[1])),data3, n = num_conditions, type = "tf", title_type = "nonpathway1", only_M0 = only_M0)
bar_update_manual(unname(unlist(non_data_labeled[1])),data_labeled, n = num_conditions, type = "tf", title_type = "nonpathway1", only_M0 = only_M0)
bar_update_manual(unname(unlist(non_FC[1])),FC, n = num_conditions, type = "tf", title_type = "nonpathway1", only_M0 = only_M0)
#second set
bar_update_manual(unname(unlist(non_data4[2])),data4, n = num_conditions, type = "tf", title_type = "nonpathway2")
bar_update_manual(unname(unlist(non_data3[2])),data3, n = num_conditions, type = "tf", title_type = "nonpathway2", only_M0 = only_M0)
bar_update_manual(unname(unlist(non_data_labeled[2])),data_labeled, n = num_conditions, type = "tf", title_type = "nonpathway2", only_M0 = only_M0)
bar_update_manual(unname(unlist(non_FC[2])),FC, n = num_conditions, type = "tf", title_type = "nonpathway2", only_M0 = only_M0)
#third set
bar_update_manual(unname(unlist(non_data4[3])),data4, n = num_conditions, type = "tf", title_type = "nonpathway3")
bar_update_manual(unname(unlist(non_data3[3])),data3, n = num_conditions, type = "tf", title_type = "nonpathway3", only_M0 = only_M0)
bar_update_manual(unname(unlist(non_data_labeled[3])),data_labeled, n = num_conditions, type = "tf", title_type = "nonpathway3", only_M0 = only_M0)
bar_update_manual(unname(unlist(non_FC[3])),FC, n = num_conditions, type = "tf", title_type = "nonpathway3", only_M0 = only_M0)

dev.off()



