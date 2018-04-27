#This program is meant to be an example of how to process data using ComBat
#In this example, we will combine data from two different file types, and correct batch effect using ComBat


#Load packages from BioConductor
bioPackages = c('Biobase','sva','PECA','ReactomePA','pcaMethods')
for (pack in bioPackages){
  if (!require(pack, character.only = T)) {
    source('https://bioconductor.org/biocLite.R')
    biocLite(pack)
    }
  library(pack, character.only = T)
}

#Load packages from CRAN
packages = c('pamr','limma','reshape2','stringr','matrixStats',
             'RColorBrewer','gplots','pheatmap','tidyverse','readr',
             'magrittr','scatterplot3d')

for (pack in packages){
  if (!require(pack, character.only = T)) {install.packages(pack)}
  library(pack, character.only = T)
}


#Selects columns of dataframe that have numeric data (usually expression levels)
dataRan <- function(pepFile = pepData){
  pepFile %>%
    sapply(., is.numeric) %>%
    which(. == TRUE) %>%
    return
}

mn = function(x) min(x, na.rm = T)
mx = function(x) max(x, na.rm = T)

#Set working directory
setwd('/Users/jozhang/Documents/Templates')

#Load Protein Names, Descriptions
#MasterData is table of Accessions for human proteins and their corresponding Gene/Description
MasterData = read_tsv('./datasets/uniprot-all.tab')
colnames(MasterData) = sapply(colnames(MasterData), function(x) sub(' ', '.', x))
colnames(MasterData)[c(1,3,4)] = c('Accession', 'Description', 'Gene')
MasterData$Gene = sapply(MasterData$Gene, function(x) unlist(strsplit(x, split = ' '))[1])
MasterData = MasterData[,c(1,3,4)]


#Load datafiles =====================================================
#There are two types of datafiles that are usually used:

#The first one is the raw spectrum data; if you have the spectrum data and matching protein table, use this function.
#     The columns used in processData may need to be adjusted to match the columns you want from your table
#     The desired columns from PSM are: Annotated.Sequence, Contaminant, Number.Protein.Groups, Master.Protein.Accessions,
#         adn the data columns X126, X127N, X127C, ... X131

PSM=read.table("./datasets/HLGSC/HLGSC_PSMs.txt", header=TRUE, sep='\t')
Proteins=read.table("./datasets/HLGSC/HLGSC_Proteins.txt", header=TRUE, sep='\t')

#List of sample names that match the order of the data columns in PSM
sNames1 = c('lgsc1','lgsc2','lgsc3','lgsc4','lgsc5','hgsc1' ,'hgsc2' ,'hgsc3' ,'hgsc4' ,'hgsc5')


###	FILTER AND AGGREGATE DATA TO PEPTIDE LEVEL 
processData <- function(PSM, Proteins, sNames, ... ){
  tablePSM=as.data.frame(PSM[,c(5,7:8,10,37:46)]) # pick the useful PSM columns
  print('Filtering Data')
  tablePSM=as.data.frame(tablePSM[tablePSM$Contaminant=="False",]) # remove peptides identified as 'contaminants'
  tablePSM=as.data.frame(tablePSM[tablePSM$Number.of.Protein.Groups=="1",]) # remove peptides that are not unique to one protein
  tableProt=as.data.frame(Proteins[,c(4,5)]) # pick the useful Proteins columns (accession and description)
  print('Merging PSM and Proteins Data')
  merge_table=merge(x=tablePSM,y=tableProt,by.x="Master.Protein.Accessions",by.y="Accession") # merge the PSM data to the Proteins description info using the accession
  print('Organizing Data Frame')
  peptide_table=as.data.frame(merge_table[,c(1,2,5:15)]) # remove Std and Contaminant columns
  peptide_table$Gene=sub(".*GN=(.*?)( .*$)","\\1",peptide_table$Description) # extract the gene info from the description column
  peptide_table$Description=sub(" OS=Homo sapiens GN=(\\w+(-\\w+)*) PE=(\\d+) SV=(\\d+)","",peptide_table$Description) # remove unneccessary description info
  peptide_table$Annotated.Sequence=sub("^.*[(.*).](.*?)([.(.*).].*$)","\\1",peptide_table$Annotated.Sequence) # extract just the peptide sequence from the annotated sequence column
  peptide_table$Annotated.Sequence=toupper(peptide_table$Annotated.Sequence) # convert all peptide sequences to upper case
  peptide_table$Gene=toupper(peptide_table$Gene) # convert all gene names to upper case
  peptide_table=peptide_table[,c(1,14,13,2,3:12)] # reorder the columns
  print('Aggregating PSM Data to Peptide Data')
  peptide_table=aggregate(cbind(X126,X127N,X127C,X128N,X128C,X129N,X129C,X130N,X130C,X131)~Master.Protein.Accessions+Gene+Description+Annotated.Sequence,peptide_table,median,na.action=na.pass,na.rm=TRUE) # aggregate the PSM data to the peptide-level taking the median PSM value as representative
  print('Aggregated')
  names(peptide_table)[5:14]<-sNames # rename the sample columns according to the names defined above
  names(peptide_table)[c(1,4)] = c('Accession', 'Sequence')
  print(head(peptide_table))
  print(dim(peptide_table))
  return(peptide_table)
}

#Create the data table for the peptide expression data
pepData = processData(PSM = PSM, Proteins = Proteins, sNames = sNames1) 

#Replace the Gene/Description data from the PSM file with the MasterData file to ensure consistency across datasets
pepData = pepData[,-which(colnames(pepData) %in% c('Gene','Description'))]

#The second type of data file is peptide expression tables (basically the ones that are created by processData in .txt format)
#     These are much easier to use.


#Set of datafiles to load
files = c('./datasets/tcc/ch_apr2017_TCC_Set1_peptideSet.txt',
          './datasets/tcc/ch_apr2017_TCC_Set2_peptideSet.txt',
          './datasets/tcc/ch_apr2017_TCC_Set3_peptideSet.txt')


#This program creates a dataframe (df) with the sample names, and the batch (i) they came from
batchDf = function(df, i=1){
  Sample = names(df)[dataRan(df)]
  output = data.frame(Sample)
  batch = rep(i, length(dataRan(df)))
  output$batch = batch
  return(output)
}

#pheno is the dataframe for information that can be used to describe each sample, i.e. batch, type, histology
pheno = batchDf(pepData, i = 1) #Create the information for the first dataset

#read through the other datasets (files) and add them to the existing pepData.
#   At the same time, use batchDf to append pheno, with the added samples and the batch they came from

for (i in 1:length(files)) {
  df = read_tsv(files[i])
  df = df[,-length(df)] #Remove the pools
  pepData = merge(pepData, df, all = TRUE, by = c('Accession', 'Sequence'))
  add_batch = batchDf(df = df, i = i+1)
  pheno = rbind(pheno, add_batch)
}
pepData[pepData == 0] <- NA #expression values should not be exactly 0; replaced with NA


#Remove HGSC4, it was a bad sample.
pepData = pepData[,-which(colnames(pepData) == 'hgsc4')]
pheno = pheno[-which(pheno$Sample == 'hgsc4'),]


#Add in Gene names and descriptions for each protein
pepData = merge(MasterData, pepData, by = 'Accession', all.y = T)
colnames(pepData)[dataRan(pepData)] = unlist(lapply(pheno$Sample,toString))

#batches is the list of all batches
batches = as.numeric(unique(pheno$batch))
#num_b is the number of batches
num_b = length(batches)

#batchNA tells us if a peptide is mostly missing from one batch
#ComBat will raise an error if the data for a peptide is completely missing from one batch
# Errors might also be raised if some, but not all, intensity levels are missing from a batch
# The current threshold for a peptide to be missing data from a batch is if it has more than 3 missing values,
#     but feel free to play around with the cutoff.
batchNA = function(probe){
  missing = F
  for (batch in levels(factor(pheno$batch))){
    if (!missing){
      sel = sapply(pheno[which(pheno$batch == batch),]$Sample, toString)
      missing = or(missing, sum(is.na(probe[sel])) > 3)
    }
  }
  return(missing)
}
#cleaned_data will be the peptide expression data for peptides that were present across all batches
cleaned_data = pepData[!apply(pepData, 1, batchNA),]


#Normalize all samples so that each one has the same total intensity; 
#We assume that samples should have roughly the same protein content

normPEP <- function(pepFile){
  pep=pepFile
  ran = dataRan(pep)
  norm_table=as.data.frame(colSums(pep[,ran],na.rm=TRUE,dims=1))
  max=max(norm_table)
  norm_table[c(1:length(ran)),]=norm_table[c(1:length(ran)),]/max
  pep_mat = as.matrix(pep[,ran])
  pep[,ran] = t(t(pep_mat)/norm_table[c(1:length(ran)),])
  print(head(pep))
  print(dim(pep))
  return(pep)
  
}

#Normalize all peptides to the same total log-intensity;
#normalizing peptides does not change the differential expression, and 
#     removes over-correlation between samples caused by large shifts in total expression between peptides
normPeptides = function(pepFile){
  df = pepFile
  dat = as.matrix(df[,dataRan(df)])
  dat = log2(dat)
  dat = t(scale(t(dat), scale = F))
  dat = 2**dat
  df[,dataRan(df)] = dat
  return(df)

}

#pepData = normPeptides(pepData)
cleaned_data = normPEP(cleaned_data)
cleaned_data = normPeptides(cleaned_data)


#Assign other information to each sample, such as which hospital the sample came from or what type of tumour it is

#The tcc dataset had high grade and tcc samples; the table matches samples to their type
#Load the tracking table to match samples to phenotypes
tracker = read_tsv('./datasets/tcc/tcc_tracker.txt')


#track_search matches tcc samples with their type
track_search = function(sam,name = 'Classifier'){
  sampl = toString(sam)
  if (grepl('pool', sam)){
    ans = "pool"
  } else {if (nrow(tracker[tracker$Sample == sampl,]) == 1) {
    ans = tracker[tracker$Sample == sampl,][,toString(name)]
  }else{
    ans = tracker[,toString(name)][[1]]
  }
  }
  return(ans)
}

pheno$type = sapply(sapply(pheno$Sample, function(x) track_search(x, "Type")), toString)

#the HLGSC dataset did not have a tracking file, but the names correspond to if it was a high grade or low grade tumour
#Therefore we edit the types in pheno to reflect this
pheno[!grepl('tcc', pheno$Sample),]$type = sub('(\\d+)', '', pheno[!grepl('tcc', pheno$Sample),]$Sample)



#Create a table of protein values ============================================================================================


med2 <- function(vec) {
  return (median(vec, na.rm = TRUE))
}

#pepToPro takes the peptide expression table and converts the peptide expression data to protein expression data
#Assumes that protein expression is the median value of the peptide expression
pepToPro <- function(pepTable = pepData){
  molten.pep = melt(pepTable, id=c('Accession', 'Sequence', 'Gene', 'Description'))
  proData = dcast(molten.pep, Accession+Gene+Description~variable, med2)
  pepNum = as.data.frame(table(unlist(pepTable[,'Accession'])))
  names(pepNum) = c('Accession', 'pepNum')
  proData = merge(pepNum, proData, by='Accession')
  print(head(proData))
  print(dim(proData))
  return(proData)
}

#Takes protein expression data and returns a matrix with only the expression values, with Accessions as the row names
prot_clean = function(proFile = proData){
  mat = as.matrix(proFile[,-c(1:4)])
  rownames(mat) = proFile$Accession
  return(mat)
}

#Takes peptide expression data and returns a matrix with only the expression values, with Accession-Sequence as the row names
pep_clean = function(pepFile){
  mat = as.matrix(pepFile[,dataRan(pepFile)])
  rownames(mat) = paste(pepFile$Accession, pepFile$Sequence, sep = '-')
  return(mat)
}
pre_data = pep_clean(cleaned_data)

#Run ComBat on logged data to prevent negative values for peptide expression
pre_data = log2(pre_data)

#scale makes the total logged expression value for each sample 0, and the variance for each sample 1
pre_data = scale(pre_data)

head(pre_data)




#Run ComBat on data ===============================================================
#model is the basic experimental design that we wish to correct with
#~1 assumes that no factors affect expression data
# If we add factors such as ~as.factor(type), ComBat will correct the data under the assumption
#     that tumour type affects expression data
# This is useful if batches do not have an even distribution of tumour types, at the cost of biasing the data correction
modcombat = model.matrix(~as.factor(type), data = pheno)
adjusted_data = ComBat(dat = pre_data, batch = pheno$batch, mod = modcombat)

#turn expression values to unlogged data
adjusted_data = 2**(adjusted_data)
pre_data = 2**(pre_data)

#add back columns of Gene/Accession/Sequence/Description
pre_df = cleaned_data

adjusted_df = data.frame(matrix(unlist(strsplit(rownames(adjusted_data),'-')), ncol = 2, byrow = T), adjusted_data)
rownames(adjusted_df) = c()
colnames(adjusted_df)[c(1,2)] = colnames(pepData)[c(1,4)]
colnames(adjusted_df) = sapply(colnames(adjusted_df), function(x) sub('\\.', '-', x))

adjusted_df = merge(MasterData, adjusted_df, by = 'Accession', all.y = T)

#set directory to where you want to store graphs/processed data
setwd('./Graphs/tcc')


#heatmap of data before batch correction ===================================================================
###	format data for correlation heat map
data = pre_data # set peptide intensity data frame as a matrix
data=log2(data) # log2-transform the data table
data=scale(data)
head(data)
dim(data)
###	pick a correlation method for your data
plot=cor(data,method="spearman", use = 'complete.obs') # use spearman correlation function to make the correlation table

plot = round(plot,2)

#Create coloured labels for batch and tumour type
batches = sapply(pheno$batch, toString)
annot_col = data.frame(row.names = pheno$Sample, Batch = batches, Type = pheno$type)
batchcols = brewer.pal(length(unique(batches)), 'Dark2')
names(batchcols) = unique(batches)
typecols = brewer.pal(9, 'Set1')[1:length(unique(pheno$type))]
names(typecols) = unique(pheno$type)

mycolors = list(Batch = batchcols, Type = typecols)

#Create of pdf of the correlation heatmap
#Clustering method lets you pick different ways to cluster the data
#   #Includes ward.D2, complete, average, and any other method that hclust() can use
pdf('TCC_PrC.pdf', onefile = F)
pheatmap(plot,
         col=rev(brewer.pal(11, 'RdBu')),
         scale = 'none',
         clustering_method = 'ward.D2',
         legend = T,
         annotation_col = annot_col,
         annotation_colors = mycolors,
         #clustering rows makes a cleaner-looking graph, at the cost of computation time
         cluster_rows = T,
         treeheight_row = 0,
         main = 'Spearman Correlation Heatmap',
         show_rownames = F,
         border_color = NA)

dev.off()

#heatmap of data after batch correction========================================================================

data2 = adjusted_data
data2=log2(data2) # log2-transform the data table
data2=scale(data2)
head(data2)
dim(data2) 
plot2=cor(data2,method="spearman", use = 'complete.obs') # use spearman correlation function to make the correlation table

plot2 = round(plot2,2)


pdf('TCC_PoC.pdf', onefile = F)
pheatmap(plot2,
         col=rev(brewer.pal(11, 'RdBu')),
         scale = 'none',
         clustering_method = 'ward.D2',
         legend = T,
         annotation_col = annot_col,
         annotation_colors = mycolors,
         cluster_rows = T,
         treeheight_row = 0,
         main = 'Spearman Correlation Heatmap',
         show_rownames = F,
         border_color = NA)

dev.off()


#Remove least variable peptides =============================================================
#Looking at peptides that have a higher variance tends to improve clustering by type
#var_cutoff is the minimum variance, adjust probs to change your percentile cut-off

variance = apply(adjusted_data, MARGIN = 1, function(x) var(x, na.rm = T))
adjusted_df$Variance = variance
var_cutoff = quantile(variance, probs = c(0.4), names = FALSE)
top_pepData= adjusted_df[adjusted_df$Variance >= var_cutoff,][,-length(adjusted_df)]
rownames(top_pepData) = c()



#look at heatmap for only the most variable peptides===============================

data3 = top_pepData # set ComBat adjusted protein intensity
data3 = as.matrix(data3[, dataRan(data3)])
rownames(data3) = paste(top_pepData$Accession,top_pepData$Sequence,sep = '-')
data3=log2(data3) # log2-transform the data table
data3=scale(data3)

head(data3)
dim(data3)


pdf('TCC_Large_model_top30.pdf', onefile = F)
pheatmap(data3,
         col=rev(brewer.pal(11, 'RdBu'))[c(3:6, 8:10)],
         breaks=  seq(-3, 3, length.out = 8),
         scale = 'none',
         clustering_method = 'ward.D2',
         legend = T,
         annotation_col = annot_col,
         annotation_colors = mycolors,
         cluster_rows = T,
         treeheight_row = 0,
         main = 'Peptide Intensity Heatmap',
         show_rownames = F,
         border_color = NA)

dev.off()

#Create subgroups based on clustering, compare correlation to known subgroups============

#Find new subgroups based on clustering
data4 = top_pepData[complete.cases(top_pepData),] # set ComBat adjusted protein intensity
data4 = as.matrix(data4[, dataRan(data4)])
data4=log2(data4) # log2-transform the data table
data4=scale(data4)

plot4 = cor(data4, method = 'spearman', use = 'complete.obs')

#Decide how many clusters to group data into
lim = 10
x_dat = c(1:lim)
y_dat = rep(0,lim)
for (k in x_dat){
  kclust = kmeans(t(data4), centers = k, nstart = 100)
  print(k)
  y_dat[k] = kclust$tot.withinss
}

#Plot of data roughly shows how much of the variance in the data is not accounted for by the clusters
#The more clusters you make, the more variance that is accounted for
#However, there are diminishing returns with more clusters; pick the number of clusters k to be roughly where diminishing returns set in
plot(x_dat, y_dat, main = 'Distortion vs. k', xlab = 'Number of clusters', ylab = "Distortion")
lines(x_dat,y_dat)

curvature = function(x,y){
  sl1 = (y[2] - y[1])/(x[2]-x[1])
  sl2 = (y[3] - y[2])/(x[3]-x[2])
  return (sl1-sl2)/(x[3] - x[1])
}

x2 = c(2:(lim-1))
y2 = sapply(x2, function(x) (curvature(x_dat[(x-1):(x+1)], y_dat[(x-1):(x+1)])))

#Plot of the estimate of the second derivative of the 'Distortion vs. k' plot. 
#Basically another graph to help pick k
#Look for exceptionally low points; those are points where diminishing returns start to set in
plot(x2,y2, col = 'red', xlab = 'k', ylab = '2nd derivative of distortion')

#Cluster the data on your chosen number of clusters
kclust = kmeans(t(data4), centers = 4, nstart = 500)

#Record clustering data in pheno
pheno$kclust = paste('k', kclust$cluster, sep = '')

#Test two groups against each other by their average correlation
X = levels(factor(pheno$kclust))
Y = unique(pheno$batch)

#Adjust ave_cor to use the columns of pheno that X and Y are from
ave_cor <- function(grp1, grp2){
  x = unlist(which(pheno$kclust == grp1))
  y = unlist(which(pheno$batch == grp2))
  return(mean(plot4[x,y]))
}

#Compare correlations between clusters
comp = sapply(Y, function(y) sapply(X, function(x) ave_cor(x,y)))
comp = round(comp, 3)

#Test how well the clusters and tumour typing line up
conting_k = as.data.frame.matrix(table(pheno$kclust, pheno$type))

chi_k1 = chisq.test(conting_k)
chi_k2 = chisq.test(conting_k/3)


heatmap.2(as.matrix(conting_k),
          dendrogram = 'none',
          cellnote = as.matrix(conting_k),
          notecex = 1.5,
          trace = 'none',
          Rowv = F,
          notecol = 'black',
          key.xlab = 'Count',
          main = 'Contingency table',
          col = brewer.pal(11, 'OrRd'))

heatmap=heatmap.2(comp,
                  Colv = F,
                  Rowv = F,
                  #cellnote=round(comp,2),
                  notecex = 2,
                  notecol='black',
                  main = "Cluster correlation",
                  key.xlab= "Averaged Spearman Correlation",
                  col=rev(brewer.pal(11,"RdBu")),
                  dendrogram = 'none',
                  trace='none',
                  margins =c(8,4)
)

#Re-run clustering with hclust() ==============================================
data4 = top_pepData
data4 = as.matrix(top_pepData[,dataRan(top_pepData)])
samp_dist = dist(t(data4))
samp_clust = hclust(samp_dist, method = 'complete')


num_c = 3

#give the results of an hclust, follow through the merge steps in hclust to give the desired number of clusters (num_c)
clusters = (rep(list(list()), dim(data4)[1]))

#hclust has a 'merge' section that lists which two clusters are combined for that step.
#   a positive number refers to a pre-made cluster that will be merged,
#   while a negative number refers to a single sample to be added
follow = function(x){
  if(x < 0){
    ret = abs(x)
  }
  if(x > 0){
    ret = clusters[[x]]
    #clusters[[x]] = list()
  }
  return(ret)
}

#follow enough merge steps until num_c clusters are left
for (i in 1:(dim(data4)[2] - num_c)){
  instruct = samp_clust$merge[i,]
  elem = append(follow(instruct[1]), follow(instruct[2]))
  if(instruct[1] > 0){
    clusters[[instruct[1]]] = list()
  }
  if(instruct[2] > 0){
    clusters[[instruct[2]]] = list()
  }
  clusters[[i]] = elem
}

non_zero = function(lst){
  return(sum(unlist(lst))> 0)
}

#if any single samples were not clustered, put them into their own clusters
clusters = Filter(non_zero, clusters)
singles = Filter(function(x) !(x %in% unlist(clusters)), 1:dim(data4)[2])
clusters = append(clusters, singles)

cluster = rep(0, dim(data4)[2])
for (i in 1:length(clusters)){
  cluster[unlist(clusters[i])] = i
}

#Record hclust results in pheno
pheno$hclust = paste('h', cluster, sep = '')

conting_h = as.data.frame.matrix(table(pheno$hclust, pheno$type))
chi_h1 = chisq.test(conting_h)
chi_h2 = chisq.test(conting_h/3)

#check agreement between k-means clustering and hierarchical clustering
compar = as.data.frame.matrix(table(pheno$hclust, pheno$kclust))


#Run PCA to check separation of clusters ================================
kclust_levels = levels(factor(pheno$kclust))
kclust_num = sapply(pheno$kclust, function(x) which(x==kclust_levels)[1])
kclust_type = brewer.pal(length(kclust_levels), 'Set1')[kclust_num]

hclust_levels = levels(factor(pheno$hclust))
hclust_num = sapply(pheno$hclust, function(x) which(x==hclust_levels)[1])
hclust_type = brewer.pal(length(hclust_levels), 'Set1')[hclust_num]

type_levels = unique(pheno$type)
type_num = sapply(pheno$type, function(x) which(x==type_levels)[1])
typeCol = brewer.pal(length(type_levels), 'Set1')[type_num]

top_data = as.matrix(top_pepData[,dataRan(top_pepData)])
pca_data = pca(t(adjusted_data), nPcs = 10)


pca1 = scores(pca_data)[,1]
pca2 = scores(pca_data)[,2]
pca3 = scores(pca_data)[,3]


plot(pca1, pca2, main = "PCA", 
     xlab = 'PCA1', ylab = 'PCA3', col = kclust_type, pch = 16)
legend('bottomleft',
       #legend = paste("Batch", levels(factor(pheno$batch)), sep = ' '),
       legend= kclust_levels,
       col = brewer.pal(8, 'Set1')[1:length(kclust_levels)],
       lty = 1,
       lwd = 10)

scatterplot3d(x = pca1, y = pca2, z = pca3, color = typeCol_type, pch = 16, type = 'h')
legend('bottomright',
       #legend = paste("Batch", levels(factor(pheno$batch)), sep = ' '),
       legend= type_levels,
       col = brewer.pal(4, 'Set1')[1:length(type_levels)],
       lty = 1,
       lwd = 10)


#Compare two subgroups with PECA ==========================

exprs = adjusted_df#set data frame
exprs[,dataRan(exprs)] = exprs[,dataRan(exprs)]-1
group1<- pheno[pheno$type == 'tcc',]$Sample %>%
  sapply(., toString)
group2<- pheno[pheno$type == 'hgsc',]$Sample %>%
  sapply(., toString)
xnorm = PECA_df(df = exprs,id = 'Accession',group1,group2,normalize=FALSE,test='modt')	#PECA analysis of peptide-level data
xnorm$Accession = rownames(xnorm)	# set row names using Gene column
rownames(xnorm) = NULL	# remove header from row names column
exprs$pepNum = 1	# add peptide counter
anno = aggregate(pepNum~Accession+Gene+Description,data=exprs,sum,na.action=na.pass,na.rm=TRUE)	# aggregate peptides to proteins and sum peptide counts
PECA = merge(anno,xnorm,by='Accession')	# merge PECA data to protein info 


PECA$extent = abs(PECA$slr - mean(PECA$slr))/sd(PECA$slr) #Measure how many standard devations away a fold change is from the mean
PECA[PECA$p.fdr > 0.05 | PECA$pepNum < 2,]$extent = 0 #Ignore proteins that aren't significant or are only measured by one peptide


write.table(PECA[PECA$pepNum>1, ],'TCC_vs_hgsc_PECA.txt',quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)

# PECA Table Headers
# slr 		(signal log ratio) = log2(group1/group2). The protein-level fold-change value
# Gene		Gene identifier
# Accession		Uniprot accession number
# Descriptions	Description based on gene identifier
# pepNum		Number of unique peptides quantified for that protein
# t			median t-statistic of all peptides assigned to that protein 
# p			protein p-value based on sampling from a beta distribution
# p.fdr		FDR (false discovery rate) adjusted p-value (q-value)


### VOLCANO PLOT FROM PECA DATA
peca=PECA
lowest = min(peca[peca$p.fdr > 0, ]$p.fdr)
peca[peca$p.fdr == 0, ]$p.fdr = lowest/10
peca$adj.p.fdr=-log10(peca$p.fdr) # calculate the -log10(p.fdr) values (y-axis)
peca$abs.score = abs(peca$slr) * peca$adj.p.fdr
cutoff = mean(peca$abs.score) + 1*sd(peca$abs.score)
#peca=as.data.frame(subset(peca,p.fdr<0.05)) # you can remove proteins with p.fdr > 0.05 if you would like

yMax=max(peca$adj.p.fdr) # the largest value in the adjusted p.fdr column

# set the colour parameters for the volcano plot. Red = significantly up-regulated in high risk, 
#Blue = significantly down-regulated in high risk, Grey = not significantly changed
colour_parameters <- ifelse((peca$abs.score < cutoff | peca$p.fdr > 0.05), 'gray22', 
                            ifelse(peca$slr > 0,'red4',
                                   ifelse(peca$slr < 0,'blue','gray22')))
colour_opacity <- adjustcolor(colour_parameters, alpha.f = 0.5) # adjust the opacity of the points to 50%

#make curved lines showing where the cutoff is
xlim = c(min(peca$slr),max(peca$slr))
ylim = c(-log10(0.05),yMax*1.1)

y_points = seq(min(ylim), max(ylim)+20, length.out = 100)
x_points = cutoff/y_points
smoothingLine = smooth.spline(x_points, y_points, spar = 0.1)
smoothingLine2 = smooth.spline(-x_points, y_points, spar = 0.1)


pdf('VolcanoPlot_TCC_vs_HGSC.pdf')
plot(x = peca$slr, # set x-axis values to be the slr values
     y = peca$adj.p.fdr, # set y-axis values to be the adjusted p.fdr values
     col = colour_opacity, 
     main="Volcano Plot", 
     pch = 16,
     xlim = xlim, # x-axis dimensions, adjust to fit your data
     ylim = ylim, # x-axis dimensions, adjust to fit your data
     abline(h=-log10(0.05),lty=2), #-log10 value of a p-value of 0.05
     xlab="slr",
     ylab="-log10(p.fdr)"
)
lines(smoothingLine, lty = 2)
lines(smoothingLine2, lty =2 )
dev.off()


#Pathway analysis: look at 500 most up-regulated genes in each type======================
setwd('/Users/jozhang/Documents/Templates')

#Use a different set of protein data that contains the Gene ID for reactome to work
MasterData2 = read_tsv('./datasets/uniprot-all2.tab')
colnames(MasterData2)[c(1:4)] = c('Accession', 'Description', 'Gene', 'Gene.ID')
MasterData2$Gene = sapply(MasterData2$Gene, function(x) unlist(strsplit(x, ';'))[1])
MasterData2$Gene.ID = sapply(MasterData2$Gene.ID, function(x) as.integer(unlist(strsplit(x, ';'))[1]))
MasterData2 = MasterData2[,c(1:4)]

sort_peca = peca[order(-peca$extent),]
sort_peca = sort_peca[sort_peca$extent > 0,]
sort_peca = merge(MasterData2, sort_peca, by = c('Accession','Description','Gene'))

pos_peca = sort_peca[sort_peca$slr > sd(peca$slr) & sort_peca$p.fdr < 0.05,]
if (dim(pos_peca)[1] > 500){
  pos_peca = pos_peca[c(1:500),]
}

neg_peca = sort_peca[sort_peca$slr < sd(peca$slr) & sort_peca$p.fdr < 0.05,]
if (dim(neg_peca)[1] > 500){
  neg_peca = neg_peca[c(1:500),]
}


setwd('./Graphs/tcc')
#Run Reactome pathway analysis
pos_path = enrichPathway(gene = pos_peca$Gene.ID, pvalueCutoff = 0.05, readable = T)
pdf('TCC_path_dotplot.pdf', width = 10)
dotplot(pos_path, showCategory = 10)
dev.off()

neg_path = enrichPathway(gene = neg_peca$Gene.ID, pvalueCutoff = 0.05, readable = T)
pdf('HGSC_path_dotplot.pdf', width = 10)
dotplot(neg_path, showCategory = 10)
dev.off()

#Write genes for Enrichr analysis (can be done online)
write_tsv(as.data.frame(pos_peca$Gene), 'TCC_genes.txt', col_names = F)
write_tsv(as.data.frame(neg_peca$Gene), 'HGSC_genes.txt', col_names = F)

