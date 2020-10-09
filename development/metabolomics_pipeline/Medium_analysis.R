library(mixOmics)
library(MetabFUN)
library(readxl)

if (dir.exists('C:/Users/dbraas/Dropbox/projects')) dropbox <- 'C:/Users/dbraas/Dropbox/projects/'
if (dir.exists('C:/Users/Danie/Dropbox/projects')==T) dropbox <- 'C:/Users/Danie/Dropbox/projects/'
if (dir.exists('C:/Users/Daniel/Dropbox/projects')==T) dropbox <- 'C:/Users/Daniel/Dropbox/projects/'
if (dir.exists('E:/Dropbox_Metabolomics/Dropbox/projects')) dropbox <- 'E:/Dropbox_Metabolomics/Dropbox/projects/'
if (dir.exists('//10.47.221.23/d/Dropbox/projects')==T) dropbox <- '//10.47.221.23/d/Dropbox/projects/'
Abbrev <- read.csv(paste0(dropbox,'Metabolomics/Abbreviations.csv'), header = T)

Title <- paste0('Footprint-',gsub('.xls[x]?','', list.files(pattern='.xls[x]?')))

# enter cell names and cell numbers
info <- read_excel(list.files()[grep('.xls[x]?',list.files())])
info$Sample <- 1:nrow(info)
info$Condition <- gsub("^\\s+|\\s+$", "", info$Condition)        #remove leading or trailing white space
info$Condition <- gsub('/', '-', info$Condition)
info$Condition <- gsub('_', '-', info$Condition)
info$Condition <- factor(info$Condition, levels = unique(info$Condition))    #important when rearranging the order of samples
if (is.null(info$Medium)) info$Medium <- 'A'
Freq <- data.frame(table(info$Condition))
Sample.Name <- vector()
for (i in 1:length(levels(info$Condition))){
  for (j in 1:Freq$Freq[i]){
    Sample.Name <- append(Sample.Name, paste(Freq$Var1[i], j, sep='_'))
  }
}
if (is.numeric(info$Cell.Number)==F) print('Problem: You need to convert the Cell.Number column')
info$Sample.Name <- Sample.Name
samples <- info
numbers <- samples$Cell.Number
save(info, file='sample info.rdata')

wd = getwd()
folder.name <- gsub('(.)*Projects/','',wd) %>% gsub('/','-',.)
files=list.files(pattern='.csv', recursive=T)
files <- files[!grepl('Suggested|Retention time', files)]
files <- files[grep('(.)*/.(.)*.csv', files)]

data=lapply(files, read.csv, header=T, na.strings = c("N/F",'N/A'))
data=do.call(rbind, data)
data$Experiment <- folder.name

write.csv(data, 'all data raw.csv', row.names=F)

#write.csv(info, paste0(dropbox,'all_medium_raw/',
#                       folder.name, '_', 'sample info.csv'), row.names = F)
#write.csv(data, paste0(dropbox, 'all_medium_raw/',
#                       folder.name, 'all data raw-', '_',
#                       gsub(':','-',file.info(list.files(pattern='.raw', recursive=T)[2])$mtime),
#                       '.csv'),
#          row.names=F)

data <- data %>%
  select(Compound, Filename, Area) %>%
  spread(Filename, Area)

data$Compound = as.character(data$Compound)
for (i in 1:length(data$Compound)){                   #loop through IDs and fill up with Ms
  if (grepl(" M[0-9]+$", data$Compound[i])==FALSE) data$Compound[i]=paste(data$Compound[i],"M0")}
Name <- gsub(" M[0-9]+$", "", data$Compound)
name_split <- regexpr("M[0-9]+$",data$Compound)
Iso <-regmatches(data$Compound, name_split)
data <- data.frame(Name, Iso, data[,2:length(data)])
data$Iso <- factor(data$Iso, levels=as.character(unique(data$Iso)))
save(data, file='data.rdata')

#background subtraction
if (sum(grepl('[Bb]lank|water|[Bb]uffer', names(data))) > 0) {
  blank <- data[,grep('[Bb]lank|water|[Bb]uffer', names(data))]
  data <- data[,-grep('[Bb]lank|water|[Bb]uffer', names(data))]
  blank[is.na(blank)] <- 0
}

if (sum(grepl('500k_cells|50[Kk]|100[Kk]|X250k|200[kK]', names(data))) > 0) {                          ##collect all standard samples and remove from regular samples
  Standard <- data[,c(1, grep('500k_cells|50[Kk]|100[Kk]|X250k|200[kK]', names(data)))]
  data <- data[,-grep('500k_cells|50[Kk]|100[Kk]|X250k|200[kK]', names(data))]
}

for (i in 1:length(Sample.Name)) colnames(data)[i+2]=Sample.Name[i]

# analyzing the metabolite standards --------------------------------------

Std <- filter(data, grepl('Std', data$Name)) %>% 
  select(-Iso) %>% 
  gather(Sample, Value, -Name) %>% 
  filter(!(grepl('500k_cells|50[Kk]|100[Kk]|250[Kk]|X500[Kk]', Sample))) %>% 
  group_by(Name) %>% 
  mutate(Av = mean(Value, na.rm=T),
         Std = sd(Value, na.rm=T),
         CV = Std / Av,
         ID = paste0(gsub(' Std','',Name), ' ',round(CV*100,0),'%')) %>% 
  ungroup()  

pdf('Relative response of ISTDs with CV.pdf', width=14, height=10)
ggplot(Std, aes(Sample, log(Value,10)))+
  geom_point(size=3)+
  facet_wrap(~ID, scales='free')+
  theme_bw()+
  theme(text = element_text(face='bold',size=14),
        axis.text.x = element_text(angle=90, vjust=.3, hjust=1))+
  labs(x='',y='Relative response of ISTD (A.U., log10)')
dev.off()

data <- filter(data, !(grepl('Std', data$Name)))

#substitute long names and add numbers of carbons in molecule
for (i in nrow(Abbrev):1) {
  data$Name <- sub(Abbrev$Used_ID[i], Abbrev$Abb[i], data$Name, fixed=T)
}

Norv <- data %>%
  filter(Name=='Norvaline') %>%
  gather(Sample, Norv,-Name, -Iso) %>%
  .$Norv
Norv.title=paste("Norvaline Response-", Title, ".pdf",sep="")
pdf(file = Norv.title, width=12, height=9, pointsize=12)
data %>%
  filter(Name=='Norvaline') %>%
  gather(Sample, Value,-Name, -Iso) %>%
  mutate(Sample=factor(Sample, levels=unique(Sample))) %>%
  select(Sample, Value) %>%
  ggplot(., aes(Sample, Value)) +
  geom_point(size=3) +
  geom_line(aes(as.integer(Sample), Value), color='blue') +
  labs(list(x='Sample Name',y='Response',title='Norvaline Response')) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
dev.off()

if (length(Norv)==0) Norv <- rep(1,nrow(samples))

norm <- data
#for (i in 3:length(norm)){
#  t <- norm[norm$Name=="Norvaline",i]
#  norm[i] = norm[i] / t
#}
norm=norm[!norm$Name=="Norvaline",]     #remove Norvaline from list

#for (i in 1: length(numbers)){
#  for (j in 1: nrow(norm)) norm[j,i+2]=norm[j,i+2]/numbers[i]
#}

data2 <- norm %>%
  gather(Exp, Value, -Name, -Iso) %>%
  separate(Exp, c('Condition', 'Exp'), sep='_') %>%
  left_join(., Abbrev, by=c('Name'='Abb'))
data2$Condition <- factor(data2$Condition, levels=levels(info$Condition))
data2$Exp <- paste0('Exp', data2$Exp)
NA_Names <- data2 %>%
  group_by(Name) %>%
  mutate(NA_list=sum(is.na(Value)),
         NA_potential=n()) %>%
  filter(NA_list==NA_potential) %>%
  distinct()
data2 <- data2 %>%
  filter(!(Name %in% NA_Names$Name))
write.csv(data2, file=paste0(Title,"-all data normalized.csv"), row.names=FALSE)
save(data2, file='data2.rdata')

amounts <- data2 %>%
  select(Name, KEGG.ID, Condition, Iso, Exp, Value) %>%
  group_by(Name, Condition, Exp) %>%
  mutate(Amount=sum(Value, na.rm=T)) %>%
  ungroup() %>%
  filter(Iso=='M0')
amounts$Amount[amounts$Amount==0] <- NA     #need to convert 0 to NA for Av calculation

amounts <- amounts %>%
  group_by(Name, Condition) %>%
  mutate(Av=mean(Amount, na.rm=T),
         Std=sd(Amount, na.rm=T),
         CV=sd(Amount, na.rm=T)/mean(Amount, na.rm=T)) %>%
  select(Name, KEGG.ID, Condition, Exp, Amount, Av, Std, CV) %>%
  ungroup() %>%
  arrange(Condition, Name)

test1 <- split(amounts, amounts[c('Name','Condition')])
NA_function <- function(x) if (sum(is.na(x)) < length(x)) return(x=x) else return(x=rep(0, length(x)))
new.Value <- as.vector(sapply(test1, function(x) NA_function(x$Amount)))
if (class(new.Value)=='list') new.Value <- unlist(new.Value)

amounts <- amounts %>%
  mutate(Amount=new.Value)

data8=split(amounts, amounts[,1])
ANOVA=suppressWarnings(sapply(data8, function(x) anova(aov(x$Amount~x$Condition))$Pr[1]))
ANOVA=rep(ANOVA, each=length(unique(info$Condition)))

amounts <- amounts %>%
  arrange(Name) %>%
  spread(Exp, Amount) %>%
  cbind('ANOVA'=ANOVA, 'Sig'=NA)

for (i in 1:nrow(amounts)){
  if (amounts$ANOVA[i] == "NaN") amounts$Sig[i]=""
  else if (amounts$ANOVA[i] <= 0.001) amounts$Sig[i]="***"
  else if (amounts$ANOVA[i] <= 0.01) amounts$Sig[i]="**"
  else if (amounts$ANOVA[i] <= 0.05) amounts$Sig[i]="*"
  else amounts$Sig[i]=""
}
write.csv(amounts, paste0(Title,"-Amounts",'.csv'), row.names=FALSE)
save(amounts, file='amounts.rdata')

plotname=paste(Title,"-plots-plus blank.pdf", sep='')
pdf(file = plotname, width=14, height=10, pointsize=12)

amounts %>%
  mutate(Name = paste(Name, Sig, sep = ' ')) %>%
  ggplot(., aes(Condition, Av, group=Condition, fill=Condition)) +
    geom_bar(position="dodge", stat="identity", colour="black", width=0.9) +
    facet_wrap( ~ Name, scales="free", nrow=floor(sqrt(length(unique(amounts$Name))))) +
    theme_bw() +
    labs(list(x="", y="Relative amounts", title=Title, fill=element_blank())) +
    theme(
      plot.title=element_text(size=20, face="bold", vjust=2),         #sets title properties
      axis.title=element_text(size=16, lineheight=20, face="bold"),   #sets theme for axis font
      axis.text=element_text(size=11, face="bold"),
      axis.text.x=element_blank(),
      legend.title=element_text(face="bold", size=12),
      legend.text=element_text(face="bold",size=12),                  #sets legend text
      strip.text=element_text(face="bold", size=15),           #sets theme for title in facets
      panel.grid.major=element_blank()) +
      scale_fill_manual("Conditions", values = colors)  +
    geom_errorbar(aes(ymin=Av, ymax=Av+Std), position=position_dodge(0.9), width=.2)
dev.off()

amounts2 <- amounts %>%
  select(Name, KEGG.ID, Condition, contains('Exp')) %>%
  gather(Exp, Amount, -Name, -KEGG.ID, -Condition) %>%
  arrange(Name, Condition) %>%
  unite(Condition_Exp, c(Condition,Exp), sep='_') %>%
  mutate(Condition_Exp=factor(Condition_Exp, levels=unique(Condition_Exp))) %>%
  #  mutate(Condition_Exp = Condition:Exp, Condition=NULL, Exp=NULL) %>%
  spread(Condition_Exp, Amount) %>%
  ungroup()

fresh <- info[grep('unspent|[Bb]lank|[Ff]resh|[Cc]ontrol|[Mm]edium|[Mm]edia', info$Condition),]$Sample  #find blank medium samples

if (length(which(sapply(amounts2[,3:length(amounts2)], sum, na.rm=T)==0)) > 0){
  amounts2 <- amounts2[-(which(sapply(amounts2[,3:length(amounts2)], sum, na.rm=T)==0)+2)]   #remove columns with only NAs
}  

for (i in seq_len(length(unique(info$Medium)))) {                                 #create variable with average amount in blank medium
  assign(LETTERS[i],
         apply(amounts2[,fresh[c(i*3-2, i*3-1, i*3)]+2], 1, mean, na.rm=T))
}
#change NAs in 'fresh' variables to 0
samples <- info[-grep('unspent|[Bb]lank|[Ff]resh|[Cc]ontrol|[Mm]edium|[Mm]edia', info$Condition),]
amounts3 <- amounts2

for (i in seq_len(length(unique(samples$Medium)))){
  spent <- samples[grep(LETTERS[i], samples$Medium),]$Sample + 2
  amounts3[,spent] <- amounts3[,spent]-get(LETTERS[i])
}
amounts3 <- amounts3[,-(fresh+2)]
amounts3[is.na(amounts3)] <- 0

condition = as.character(samples$Sample.Name)
numbers=samples$Cell.Number

for (i in 1: length(numbers)){
  for (j in 1: nrow(amounts3)) amounts3[j,i+2]=amounts3[j,i+2]/numbers[i]
}
write.csv(amounts3,paste0(Title,'-Normalized Amounts.csv'), row.names=F)
save(amounts3, file='Normalized Amounts.rdata')

#making a heatmap
amounts3 <- amounts3[,!(colSums(amounts3[3:length(amounts3)], na.rm=T)==0)]
amounts3 <- amounts3[!(rowSums(amounts3[3:length(amounts3)], na.rm=T)==0),]
data5=as.matrix(amounts3[3:length(amounts3)])
rownames(data5)=amounts3[,1]
colnames(data5)=names(amounts3[3:length(amounts3)])

ann=data.frame("Condition"=as.character(samples$Condition),'Cell Number'=numbers)
Sig.color=c("blue", "red")
names(Sig.color)=unique(samples$Type.1)
Norv.color = c("white", "blue")
Cell.Number.color= c("white", "green")
names(colors)=unique(samples$Condition)
ann_colors = list(Condition=colors[1:length(unique(samples$Condition))], Norvaline = Norv.color, Cell.Number = Cell.Number.color)
rownames(ann)=colnames(data5)

write.csv(data5, paste0(Title,'-heatmap data.csv'), row.names=T)
data5[is.na(data5)] <- 0     #replace NAs
heatmap.title=paste(Title, "-All-Clustered Heatmap.pdf", sep='')
dev.off()

pheatmap::pheatmap(data5, cluster_row=T, cluster_col=T, clustering_distance_rows='correlation', clustering_distance_cols='correlation',
         color = colorRampPalette(c("blue", "white", "red"))(100), border_color="black", scale="row",
         cellwidth = 20, cellheight = 10, annotation=ann, annotation_colors = ann_colors, show_colnames = F,
         main=Title, filename=heatmap.title)
dev.off()

amounts3[,3:length(amounts3)][amounts3[,3:length(amounts3)]==0] <- NA

amounts3 <- amounts3 %>%
  gather(Condition_Exp, Amount, -Name, -KEGG.ID) %>%
  separate(Condition_Exp, c('Condition', 'Exp'), sep='_') %>%
  mutate(Condition=factor(Condition, levels=unique(.$Condition))) %>%
  group_by(Name, Condition) %>%
  mutate(Av=mean(Amount),
         Std=sd(Amount),
         CV=Std/Av) %>%
  ungroup()

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

plotname=paste(Title,"-plots.pdf", sep='')
pdf(file = plotname, width=14, height=10, pointsize=12)

amounts3 %>%
  mutate(Name = paste(Name, Sig, sep=' ')) %>%
  ggplot(., aes(Condition, Av, group=Condition, fill=Condition)) +
    geom_bar(position="dodge", stat="identity", colour="black", width=0.9) +
    facet_wrap( ~ Name, scales="free", nrow=floor(sqrt(length(unique(amounts3$Name))))) +
    theme_bw() +
    labs(list(x="", y="Relative amounts", title=Title, fill=element_blank())) +
    theme(
      plot.title=element_text(size=20, face="bold", vjust=2),         #sets title properties
      axis.title=element_text(size=16, lineheight=20, face="bold"),   #sets theme for axis font
      axis.text=element_text(size=11, face="bold"),
      axis.text.x=element_blank(),
      legend.title=element_text(face="bold", size=12),
      legend.text=element_text(face="bold",size=12),                  #sets legend text
      strip.text=element_text(face="bold", size=15),           #sets theme for title in facets
      panel.grid.major=element_blank()) +
    scale_fill_manual("Conditions", values = colors)  +
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
rm(colors)

data3 <- make_MID(data2)

metab <- filter(data3, Iso=='M1') %>% 
  .$Name %>% 
  as.character()

pdf(paste0(Title,'-MID plots.pdf'), width=14, height=10)
bar(unique(data3$Name), filter(data3, Name %in% metab))
dev.off()
save.image(paste0(Title,'.rdata'))
