#This program is a modified way of using ComBat to analyse the entire set of peptide data,
#     instead of just peptides that are present across all batches.
#See powerpoint for additional theory/motivation

#Load packages from BioConductor
bioPackages = c('Biobase','sva','PECA', 'EmpiricalBrownsMethod')
for (pack in bioPackages){
  if (!require(pack, character.only = T)) {
    source('https://bioconductor.org/biocLite.R')
    biocLite(pack)
  }
  library(pack, character.only = T)
}

#Load packages from CRAN
packages = c('pamr','limma','reshape2','stringr',
             'matrixStats','tidyverse','readr','magrittr','RColorBrewer',
             'gplots','plyr','PECA','pheatmap','survcomp','gtools')

for (pack in packages){
  if (!require(pack, character.only = T)) {install.packages(pack)}
  library(pack, character.only = T)
}

#Set up some useful functions
dataRan <- function(pepFile = pepData){
  pepFile %>%
    sapply(., is.numeric) %>%
    which(. == TRUE) %>%
    return
}

mn = function(x) min(x, na.rm = T)
mx = function(x) max(x, na.rm = T)

#set working directory
setwd('/Users/jozhang/Documents/Templates')

#Load Protein Names, Descriptions
MasterData = read_tsv('./datasets/uniprot-all.tab')
colnames(MasterData) = sapply(colnames(MasterData), function(x) sub(' ', '.', x))
colnames(MasterData)[c(1,3,4)] = c('Accession', 'Description', 'Gene')
MasterData$Gene = sapply(MasterData$Gene, function(x) unlist(strsplit(x, split = ' '))[1])
MasterData = MasterData[,c(1,3,4)]

#Load all the files
files = c('./datasets/hgsc/ch_apr2017_HGSC_Set1_peptideSet.txt',
          './datasets/hgsc/ch_apr2017_HGSC_Set2_peptideSet.txt',
          './datasets/hgsc/ch_apr2017_HGSC_Set3_peptideSet.txt',
          './datasets/hgsc/ch_apr2017_HGSC_Set4_peptideSet.txt',
          './datasets/hgsc/ch_apr2017_HGSC_Set5_peptideSet.txt',
          './datasets/hgsc/ch_apr2017_HGSC_Set6_peptideSet.txt')


#Read through first dataset to set up data frame
pepData = readr::read_tsv(files[1])
pepData = pepData[,-length(pepData)] #remove pool

#for the HGSC dataset, batchDf is modified because the names of samples are repeated in the datasets.
#This is because this experiment ran each sample three times, so the names have -x, -y, or -z added to differentiate them in pheno.
batchDf = function(df, i=1, multi = FALSE){
  if(multi){
    bar = '-x'
    if (i%%3 == 2){
      bar = '-y'
    }
    if (i%%3 == 0){
      bar = '-z'
    }
    Sample = paste(names(df)[dataRan(df)], bar, sep = '')
  }else{
    Sample = names(df)[dataRan(df)]
  }
  output = data.frame(Sample)
  batch = rep(i, length(dataRan(df)))
  output$batch = batch
  return(output)
}
pheno = batchDf(pepData, i = 1, multi = T)

#read through the other datasets
for (i in 2:length(files)) {
  df = read_tsv(files[i])
  df = df[,-length(df)] #Remove the pools
  pepData = merge(pepData, df, all = TRUE, by = c('Accession','Sequence'))
  add_batch = batchDf(df = df, i = i, multi = T)
  pheno = rbind(pheno, add_batch)
}

colnames(pepData)[dataRan(pepData)] = unlist(lapply(pheno$Sample,toString))

#Check if a sequence was listed twice (i.e. for some reason was assigned to two or more proteins)
seqCheck = as.data.frame(table(unlist(pepData$Sequence)))
seqRem = seqCheck[seqCheck$Freq >1,]$Var1
#Remove all non-unique peptide sequences
pepData = pepData[!(pepData$Sequence %in% seqRem),]


pepData = merge(MasterData, pepData, by = 'Accession', all.y = T)

pepData[pepData == 0] <- NA
batches = as.numeric(levels(factor(pheno$batch)))
num_b = length(levels(factor(pheno$batch)))

#Create a 'hash code' to uniquely identity each peptide by what batches it was measured in
#If a peptide is in batches 1, 2, and 4, the code will be 2^1 + 2^2 + 2^4 = 22
batch_hash = function(probe){
  hash = 0
  for (batch in batches){
    sel = sapply(pheno[which(pheno$batch == batch),]$Sample, toString)
    pick = probe[sel]
    empty = length(pick[complete.cases(pick)]) < 4
    hash = hash + (2*!empty)^batch
  }
  return(hash)
}

pepData$hash = apply(pepData,1,batch_hash)
#In the edge case where a peptide is only present in one batch and a few samples,
#   the peptide gets removed to avoid issues with ComBat correction
pepData = pepData[pepData$hash != 0, ]

hashes = unique(pepData$hash)
hashes = hashes[order(hashes)]

#Converts the hashcode back to the T/F vector for the selected batches
#This allows us to recover which batches the peptide is present in without processing it again
hashBatch = function(hash){
  vec = as.integer(intToBits(hash))[2:(num_b+1)]
  vec = which(vec == 1)
  return(vec)
}

#Normalize sample intensity
normPEP <- function(pepFile){
  pep=pepFile
  ran = dataRan(pep)
  norm_table=as.data.frame(colSums(pep[,ran],na.rm=TRUE,dims=1))
  max=max(norm_table, na.rm = T)
  norm_table[c(1:length(ran)),]=norm_table[c(1:length(ran)),]/max
  pep_mat = as.matrix(pep[,ran])
  pep[,ran] = t(t(pep_mat)/norm_table[c(1:length(ran)),])
  return(pep)
  
}

#Normalize Probe intensity
normPeptides = function(pepFile){
  df = pepFile
  dat = as.matrix(df[,dataRan(df)]) %>%
    log2(.) %>%
    t %>%
    scale(., scale = F) %>%
    t %>%
    2**.
  df[,dataRan(df)] = dat
  return (df)
}

#Run processing and ComBat on each separate set of peptides that share the same batches
#If a peptide is only present in one batch, ComBat is not run on it, but the peptide is still normalized
for (hash in hashes){
  sub_data = pepData[pepData$hash == hash,]
  sel_batch = hashBatch(hash)
  sub_pheno = pheno[pheno$batch %in% sel_batch,]
  samp = sapply(sub_pheno$Sample, toString)
  sub_mat = sub_data[,samp] %>%
    normPEP %>%
    normPeptides %>%
    as.matrix %>%
    log2 %>%
    scale
  rownames(sub_mat) = paste(sub_data$Accession, sub_data$Sequence, sep = '_')
  #Apply ComBat if there is more than one batch
  if (length(sel_batch) != 1){
    modcombat = model.matrix(~1, data = sub_pheno)
    adj_data = ComBat(dat = sub_mat, batch = sub_pheno$batch, mod = modcombat) %>%
      2**.
  } else {
    adj_data = 2**sub_mat
  }
  sub_mat = 2**sub_mat
  raw_df = data.frame(matrix(unlist(strsplit(rownames(sub_mat),'_')), ncol = 2, byrow = T), sub_mat)
  rownames(raw_df) = c()
  colnames(raw_df)[c(1,2)] = colnames(pepData)[c(1,4)]
  colnames(raw_df) = sapply(colnames(raw_df), function(x) sub('\\.', '-', x)) 
  
  adj_df = data.frame(matrix(unlist(strsplit(rownames(adj_data),'_')), ncol = 2, byrow = T), adj_data)
  rownames(adj_df) = c()
  colnames(adj_df)[c(1,2)] = colnames(pepData)[c(1,4)]
  colnames(adj_df) = sapply(colnames(adj_df), function(x) sub('\\.', '-', x))
  if (hash == hashes[1]){
    out_df = adj_df
    pre_df = raw_df
  } else{
    out_df = rbind.fill(out_df, adj_df)
    pre_df = rbind.fill(pre_df, raw_df)
  }
  print(which(hashes == hash))
}

out_df = merge(MasterData, out_df, by = 'Accession', all.y = T)
pre_df = merge(MasterData, pre_df, by = 'Accession', all.y = T)
hash_out = merge(pepData[,c('Sequence','hash')],out_df, by = "Sequence")

#Create protein level data of the entire peptide expression data if desired

med2 <- function(vec) {
  return (median(vec, na.rm = TRUE))
}
#pepToPro takes the peptide expression table and converts the peptide expression data to protein expression data
#Assumes that protein expression is the median value of the peptide expression
pepToPro <- function(pepTable = pepData){
  molten.pep = melt(pepTable, id=c('Accession', 'Sequence', 'Gene', 'Description', 'hash'))
  proData = dcast(molten.pep, Accession+Gene+Description+hash~variable, med2)
  pepNum = as.data.frame(table(unlist(pepTable[,'Accession'])))
  names(pepNum) = c('Accession', 'pepNum')
  proData = merge(pepNum, proData, by='Accession')
  print(head(proData))
  print(dim(proData))
  return(proData)
}
#hash_prot = pepToPro(hash_out)


#Add tumour type to the pheno data
pheno$type = sub('\\d-\\w', '', pheno$Sample)

pre_data = pre_df[,dataRan(pre_df)]
out_data = out_df[,dataRan(out_df)]

setwd('./Graphs/hgsc')
#Correlation heatmap of data before batch correction ===================================================================
###	format data for correlation heat map
data = pre_data %>% # set peptide intensity data frame as a matrix
  log2 %>% # log2-transform the data table
  scale
head(data)
dim(data)
###	pick a correlation method for your data
plot=cor(data,method="spearman", use= 'complete.obs') # use pearson correlation function to make the correlation table

plot = round(plot,2)

batches = sapply(pheno$batch, toString)
annot_col = data.frame(row.names = pheno$Sample, Batch = batches, Type = pheno$type)
batchcols = brewer.pal(length(unique(batches)), 'Dark2')
names(batchcols) = unique(batches)
typecols = brewer.pal(length(unique(pheno$type)), 'Set1')
names(typecols) = unique(pheno$type)

mycolors = list(Batch = batchcols, Type = typecols)

pdf('HGSC_full_PrC.pdf', onefile = F)
pheatmap(plot,
         col=rev(brewer.pal(11,"RdBu")),
         scale = 'none',
         clustering_method = 'complete',
         legend = T,
         annotation_col = annot_col,
         annotation_colors = mycolors,
         cluster_rows = T,
         treeheight_row = 0,
         main = 'Spearman Correlation Heatmap',
         show_rownames = F,
         border_color = NA)

dev.off()

#heatmap of data after batch correction========================================================================

data2 = out_data %>%
  log2 %>%
  scale
head(data2)
dim(data2)
###	pick a correlation method for your data

plot2=cor(data2,method="spearman", use = 'complete.obs') # use spearman correlation function to make the correlation table


plot2 = round(plot2,2)


pdf('HGSC_full_PoC.pdf', onefile = F)
pheatmap(plot2,
         col=rev(brewer.pal(11,"RdBu")),
         scale = 'none',
         clustering_method = 'complete',
         legend = T,
         annotation_col = annot_col,
         annotation_colors = mycolors,
         cluster_rows = T,
         treeheight_row = 0,
         main = 'Spearman Correlation Heatmap',
         show_rownames = F,
         border_color = NA)

dev.off()

#Compare clustering between 'studies' ====================================================

#given a dataset and the desired number of clusters, clustering returns the grouping of the samples
#uses hierarchical clustering to determine clusters
#data can include non-numeric columns like accession or description, but all numeric columns must be sample data!
clustering = function(data, num){
  data = as.matrix(data[, dataRan(data)]) %>%
    log2 %>%
    scale
  distance = dist(t(data))
  clust = hclust(distance, method = 'complete')
  
  clusters = (rep(list(list()), dim(data)[1]))
  
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
  
  for (i in 1:(dim(data)[2] - num)){
    instruct = clust$merge[i,]
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
  
  clusters = Filter(non_zero, clusters)
  singles = Filter(function(x) !(x %in% unlist(clusters)), 1:dim(data)[2])
  clusters = append(clusters, singles)
  
  cluster = rep(0, dim(data)[2])
  for (i in 1:length(clusters)){
    cluster[unlist(clusters[i])] = i
  }
  return(cluster)
}
num_c = 3
#Record clustering over the entire set of peptides
pheno$hclust = paste('h', clustering(out_df, num_c), sep = '')

conting_h = as.data.frame.matrix(table(pheno$type, pheno$hclust))

#pick what to compare clustering to
contrast = 'type'
rows = unique(pheno[,contrast])
#Choose a variance cutoff using one of the two following options: best clustering versus a chosen contrast (usually tumour type),
#   or by most stable clustering (i.e. distance metric from consensus clustering is minimized)
#Pick variance cutoff by finding best clustering against contrast ===============
#Cycles through cutoffs from all peptides to top 10% of peptides by variance, and then records the variance that gave
#     clusters with the best agreement to tumour type
x = seq(0,0.9,0.05)
y = rep(0, length(x))
for (foo in c(1:length(x))){
  conts = list()
  options = permutations(n = num_c, r = num_c, v = 1:num_c)
  pb = txtProgressBar(min = 0, max = length(hashes), initial = 0, style = 3)
  for (hash in hashes){
    #select samples and peptides that match the hashcode
    exprs = hash_out[hash_out$hash == hash,][-c(1:5)]
    subpheno = pheno[pheno$batch %in% hashBatch(hash),]
    sel = sapply(subpheno$Sample, toString)
    exprs = exprs[,sel]
    variance = apply(exprs, 1, function(x) var(x, na.rm = T))
    var_cutoff = quantile(variance, probs = c(x[foo]), names = FALSE)
    rem = (variance >= var_cutoff)
    exprs = exprs[rem,]
    
    #cluster on selected data
    #number of clusters should equal number of types it is compared to
    num_s = length(unique(subpheno[,contrast]))
    #record clustering results in subpheno
    subpheno$clust = paste('c', clustering(exprs, num_s), sep = '')
    
    #creating a contingency table helps determine which variance cutoff gives the best agreement between clustering and tumour type
    cont = as.data.frame.matrix(table(subpheno[,contrast], subpheno$clust))
    #make a square matrix to see how well the clusters agree with the tumour types
    if (num_s < num_c){
      new_cont = matrix(0, ncol = num_c, nrow = num_c)
      rownames(new_cont) = rows
      for (r in rows){
        if (r %in% rownames(cont)){
          new_cont[r, ] = c(unlist(cont[r,]), rep(0, num_c - num_s))
        }
      }
      cont = as.data.frame.matrix(new_cont)
      colnames(cont) = paste('c', c(1:num_c), sep = '')
      rownames(cont) = rows
    }
    #aggregate contingency tables
    top = sum(cont)
    keep = cont
    #reorder columns of the contingency table (i.e. rename clusters) so that cluster 1 always is closest to the first contrast type,
    #   cluster 2 is always closest to the second contrast type, and so on
    for (i in c(1:dim(options)[1])){ 
      opt = cont[,c(options[i,])]
      for (j in c(1:num_c)){
        opt[j,j] = 0
      }
      opt[is.na(opt)] = 0
      if (sum(opt) < top) {
        keep = cont[,c(options[i,])]
        top = sum(opt)
      }
    }
    #record agreement between clustering
    n = which(hashes == hash)
    conts[[n]] = keep
    setTxtProgressBar(pb, n)
  }
  ave_cont = matrix(0, nrow = num_c, ncol = num_c)
  rownames(ave_cont) = rownames(conts[[1]])
  colnames(ave_cont) = paste('c', c(1:num_c), sep = '')
  for (i in c(1:num_c)){
    for (j in c(1:num_c)){
      for (k in c(1:length(hashes))){
        ave_cont[i,j] = ave_cont[i,j] + conts[k][[1]][i,j]
      }
    }
  }
  #Record how well the overall clustering agrees with tumour type using chi-squared as a metric
  y[foo] = chisq.test(ave_cont)$statistic
  cat('\n')
  print(x[foo])
}
plot(x,y, xlab = 'Variance Cutoff', ylab = 'Chi-Squared')
quant = x[which(y == max(y))]
#if more than one cutoff has the highest chi-squared, pick the one that keeps the most peptides
if (length(quant) > 1 ){
  quant = quant[1]
}

#Pick variance cutoff by finding most stable clustering====================

#The jaccard index measures how much the clustering of two samples agrees with each other
#Only compares studies where both samples were present
jaccard = function(X,Y){
  comp = cbind(X,Y)
  comp = comp[complete.cases(comp),]
  if (is.null(dim(comp))){
    inters = (comp[1] == comp[2])
    index = sum(inters)
  }else{
    inters = (comp[,1] == comp[,2])
    index = sum(inters)/(2*dim(comp)[1] - sum(inters))
  }
  
  return (1-index)
}
#The weighted jaccard index is modified to account for the fact that different peptide groups have different numbers of peptides
# weight can be any numeric vector of the same length as X and Y
jaccard_weight = function(X,Y, weight){
  comp = cbind(X,Y, weight)
  comp = comp[complete.cases(comp),]
  if (is.null(dim(comp))){
    inters = (comp[1] == comp[2])*comp[3]
    index = sum(inters)/(comp[3])
  } else {
    inters = (comp[,1] == comp[,2])*comp[,3]
    index = sum(inters)/(2*sum(comp[,3]) - sum(inters))
  }
  return (1-index)
}

weights = sqrt(table(pepData$hash))
x = seq(0,0.9,0.05)
y = rep(0, length(x))
for (foo in c(1:length(x))){
  conts = list()
  options = permutations(n = num_c, r = num_c, v = 1:num_c)
  pb = txtProgressBar(min = 0, max = length(hashes), initial = 0, style = 3)
  for (hash in hashes){
    #select samples and peptides that match the hashcode
    exprs = hash_out[hash_out$hash == hash,][-c(1:5)]
    subpheno = pheno[pheno$batch %in% hashBatch(hash),]
    sel = sapply(subpheno$Sample, toString)
    exprs = exprs[,sel]
    variance = apply(exprs, 1, function(x) var(x, na.rm = T))
    var_cutoff = quantile(variance, probs = c(x[foo]), names = FALSE)
    rem = (variance >= var_cutoff)
    exprs = exprs[rem,]
    
    #cluster on selected data
    #number of clusters should equal number of types it is compared to
    num_s = length(unique(subpheno[,contrast]))
    #record clustering results in subpheno
    clust = clustering(exprs, num_s)
    subpheno$clust = paste('c', clust, sep = '')
    #record agreement between clustering
    n = which(hashes == hash)
    if (n == 1){
      #cc is a table recording which clusters the samples were put into during each 'study'
      cc = data.frame(Sample = subpheno$Sample, clust = clust)
      colnames(cc)[2] = paste('clust', n, sep = '.')
    }else{
      add_cc = data.frame(Sample = subpheno$Sample, clust = clust)
      colnames(add_cc)[2] = paste('clust', n, sep = '.')
      cc = merge(cc, add_cc, all = T, by = 'Sample')
    }
    setTxtProgressBar(pb, n)
  }
  cc = cc[match(pheno$Sample, cc$Sample),]
  cc_data = t(as.matrix(cc[,dataRan(cc)]))
  dst = matrix(0, nrow = dim(cc_data)[2], ncol = dim(cc_data)[2])
  for (i in c(1:dim(cc_data)[2])){
    for (j in c(1:i)){
      dst[i,j] = jaccard_weight(cc_data[,i], cc_data[,j], weights)
    }
  }
  #Record how well the overall clustering agrees with tumour type using chi-squared as a metric
  y[foo] = sqrt(sum((0.5-dst)**2))
  cat('\n')
  print(x[foo])
}
plot(x,y, xlab = 'cutoff', ylab = 'Stability', main = 'Cutoff vs. Stability')
quant = x[which(y == max(y))]
#if more than one cutoff has the highest chi-squared, pick the one that keeps the most peptides
if (length(quant) > 1 ){
  quant = quant[1]
}

#Rerun clustering with the best cutoff (quant) for peptide variance ======================

#Re-run clustering with the best cut-off
options = permutations(n = num_c, r = num_c, v = 1:num_c)
conts = list()
metrics = rep(0, length(hashes)-1)
for (hash in hashes){
  exprs = hash_out[hash_out$hash == hash,][-c(1:5)]
  subpheno = pheno[pheno$batch %in% hashBatch(hash),]
  sel = sapply(subpheno$Sample, toString)
  exprs = exprs[,sel]
  variance = apply(exprs, 1, function(x) var(x, na.rm = T))
  var_cutoff = quantile(variance, probs = c(quant), names = FALSE)
  rem = (variance >= var_cutoff)
  exprs = exprs[rem,]
  #cluster on selected data
  #number of clusters should equal number of types it is compared to
  num_s = length(unique(subpheno[,contrast]))
  #record clustering results in subpheno
  clust = clustering(exprs, num_s)
  subpheno$clust = paste('c', clust, sep = '')
  
  #creating a contingency table helps determine which variance cutoff gives the best agreement between clustering and tumour type
  cont = as.data.frame.matrix(table(subpheno[,contrast], subpheno$clust))
  #make a square matrix to see how well the clusters agree with the tumour types
  if (num_s < num_c){
    new_cont = matrix(0, ncol = num_c, nrow = num_c)
    rownames(new_cont) = rows
    for (r in rows){
      new_cont[r, ] = c(unlist(cont[r,]), rep(0, num_c - num_s))
    }
    cont = as.data.frame.matrix(new_cont)
    colnames(cont) = paste('c', c(1:num_c), sep = '')
    rownames(cont) = rows
  }
  #aggregate contingency tables
  top = sum(cont)
  keep = cont
  for (i in c(1:dim(options)[1])){
    opt = cont[,c(options[i,])]
    for (j in c(1:num_c)){
      opt[j,j] = 0
    }
    if (sum(opt) < top) {
      keep = cont[,c(options[i,])]
      top = sum(opt)
    }
  }
  #record agreement between clustering
  n = which(hashes == hash)
  if (n == 1){
    #cc is a table recording which clusters the samples were put into during each 'study'
    cc = data.frame(Sample = subpheno$Sample, clust = clust)
    colnames(cc)[2] = paste('clust', n, sep = '.')
  }else{
    add_cc = data.frame(Sample = subpheno$Sample, clust = clust)
    colnames(add_cc)[2] = paste('clust', n, sep = '.')
    cc = merge(cc, add_cc, all = T, by = 'Sample')
  }

  conts[[n]] = keep
  print(n)
}


ave_cont = matrix(0, nrow = num_c, ncol = num_c)
rownames(ave_cont) = rownames(conts[[1]])
colnames(ave_cont) = paste('c', c(1:num_c), sep = '')
for (i in c(1:num_c)){
  for (j in c(1:num_c)){
    for (k in c(1:length(hashes))){
      ave_cont[i,j] = ave_cont[i,j] + conts[k][[1]][i,j]
    }
  }
}

cc = cc[match(pheno$Sample, cc$Sample),]
cc_data = t(as.matrix(cc[,dataRan(cc)]))
colnames(cc_data) = pheno$Sample
rownames(cc_data) = sapply(hashes, function(x) toString(hashBatch(x)))

#create a distance matrix dst using the weighted jaccard distance
dst = matrix(0, nrow = dim(cc_data)[2], ncol = dim(cc_data)[2])
rownames(dst) = pheno$Sample
colnames(dst) = pheno$Sample

weights = sqrt(table(pepData$hash))
for (i in c(1:dim(cc_data)[2])){
  for (j in c(1:i)){
    dst[i,j] = jaccard_weight(cc_data[,i], cc_data[,j], weights)
  }
}


for (i in c(1:dim(cc_data)[2])){
  for (j in c(1:dim(cc_data)[2])){
    if (j>i){
      dst[i,j] = NA
    }
  }
}


batches = sapply(pheno$batch, toString)
annot_col = data.frame(row.names = pheno$Sample, Batch = batches, Type = pheno$type)
batchcols = brewer.pal(length(unique(batches)), 'Dark2')
names(batchcols) = unique(batches)
typecols = brewer.pal(9, 'Set1')[1:length(unique(pheno$type))]
names(typecols) = unique(pheno$type)

mycolors = list(Batch = batchcols, Type = typecols)
#Cluster samples together based on how other the samples clustered together
#To use dst as the distance sample, call fix(pheatmap)
#replace line 92: "tree_col = cluster_mat(t(mat), distance = clustering_distance_cols, 
#   method = clustering_method)" with hclust(as.dist(dst), method = 'average')
pdf(paste('Clusters_weight_stable', num_c,'.pdf', sep = ''), onefile = F)
pheatmap(cc_data,
         col=rev(brewer.pal(num_c,"Set2"))[1:num_c],
         scale = 'none',
         clustering_method = 'complete',
         breaks = seq(1,num_c,length.out = (1+num_c)),
         legend = T,
         legend_breaks = c(1:num_c),
         legend_labels = c(1:num_c),
         annotation_col = annot_col,
         annotation_colors = mycolors,
         cluster_rows = F,
         treeheight_row = 0,
         main = 'Co-clustering Comparison',
         show_rownames = T,
         fontsize_row = 6,
         #cellwidth = 7,
         border_color = NA)

dev.off()

clust = hclust(as.dist(dst), method = 'ward.D2')
data = top_pepData[,dataRan(top_pepData)]


clusters = (rep(list(list()), dim(data)[1]))

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

for (i in 1:(dim(data)[2] - num_c)){
  instruct = clust$merge[i,]
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

clusters = Filter(non_zero, clusters)
singles = Filter(function(x) !(x %in% unlist(clusters)), 1:dim(data)[2])
clusters = append(clusters, singles)

cluster = rep(0, dim(data)[2])
for (i in 1:length(clusters)){
  cluster[unlist(clusters[i])] = i
}

#record final clustering into pheno data 
pheno$fclust = paste('f', cluster, sep= '')
conting_f = as.data.frame.matrix(table(pheno$type, pheno$fclust))





#New Peptide Analysis ================================================
#Run PECA on each 'study', and then use accepted meta-analysis techniques to combine the results

group1<- pheno[pheno$type == 'hgsc',]$Sample %>%
  sapply(., toString)
group2<- pheno[pheno$type == 'ccc',]$Sample %>%
  sapply(., toString)

#Run PECA on each peptide group

made_allPECA = F
for (hash in hashes){
  exprs = hash_out[hash_out$hash == hash,]
  sel = sapply(pheno[pheno$batch %in% hashBatch(hash),]$Sample, toString)
  grp1 = intersect(sel, group1)
  grp2 = intersect(sel, group2)
  if (length(grp1) != 0 & length(grp2) != 0){
    exprs[,dataRan(exprs)][,-c(1)] = exprs[,dataRan(exprs)][-c(1)] - 1
    xnorm = PECA_df(df = exprs,id = 'Accession',grp1,grp2,normalize=FALSE,test='t')	#PECA analysis of peptide-level data
    xnorm$Accession = rownames(xnorm)	# set row names using Gene column
    rownames(xnorm) = NULL	# remove header from row names column
    exprs$pepNum = 1	# add peptide counter
    anno = aggregate(pepNum~Accession,data=exprs,sum,na.action=na.pass,na.rm=TRUE)	# aggregate peptides to proteins and sum peptide counts
    add_PECA = merge(anno,xnorm,by='Accession')	# merge PECA data to protein info 
    if (length(hashBatch(hash)) < 2){
      add_PECA = add_PECA[add_PECA$pepNum > 1, ]
    }
    add_PECA$weight = length(hashBatch(hash))*add_PECA$pepNum
    if (!made_allPECA){
      all_PECA = add_PECA
      made_allPECA = T
    }
    else {
      all_PECA = rbind.fill(all_PECA, add_PECA)
    }
  }
  print(which(hashes == hash))
}
#remove protein results that had NA values
all_PECA = all_PECA[complete.cases(all_PECA$slr),]
accs = unique(all_PECA$Accession)

#create a new q-value (basically p.fdr) based on the p-values from the peca analyses over all peptide groups
all_PECA$q = p.adjust(all_PECA$p, method = 'fdr')


#Run through possible cut-offs for q-values to figure out which cut-off gives the largest number of significant proteins
x = seq(0.05, 1, 0.05)
y = rep(0, length(x))
for (foo in c(1:length(x))){
  madePECA = F
  pb = txtProgressBar(min = 0, max = 1, initial = 0, style = 3)
  for(acc in accs){
    sel = all_PECA[all_PECA$Accession == acc,]
    sel = sel[sel$q < x[foo] |sel$pepNum > 5,]
    if (dim(sel)[1] > 0){
      pepNum = sum(sel$pepNum)
      slr = weighted.mean(sel$slr, sqrt(sel$weight))
      q = combine.test(sel$q, weight = sqrt(sel$weight), method = 'fisher', na.rm = T)
      add_PECA = data.frame(Accession = acc, pepNum = pepNum, slr = slr, q = q)
      if (!madePECA){
        PECA = add_PECA
        madePECA = T
      } else {
        PECA = rbind.fill(PECA, add_PECA)
      }
    }
    setTxtProgressBar(pb, round(which(accs == acc)/length(accs), 2))
  }
  PECA$extent = abs(PECA$slr - mean(PECA$slr))/sd(PECA$slr)
  PECA[PECA$q > 0.05 | PECA$pepNum < 2,]$extent = 0
  cat('\n')
  print(x[foo])
  y[foo] = dim(PECA[PECA$extent>0,])[1]
}



cutoff = x[which(y == max(y))]
if (length(cutoff) > 1){
  cutoff = cutoff[1]
}
madePECA = F
pb = txtProgressBar(min = 0, max = 1, initial = 0, style = 3)
for(acc in accs){
  sel = all_PECA[all_PECA$Accession == acc,]
  sel = sel[sel$p.fdr < cutoff | sel$pepNum > 5,]
  if (dim(sel)[1] > 0){
    pepNum = sum(sel$pepNum)
    slr = weighted.mean(sel$slr, sqrt(sel$weight))
    p = combine.test(sel$p, weight = sqrt(sel$weight), method = 'fisher', na.rm = T)
    p.fdr = combine.test(sel$p.fdr, weight = sqrt(sel$weight), method = 'fisher', na.rm = T)
    q = combine.test(sel$q, weight = sqrt(sel$weight), method = 'fisher', na.rm = T)
    add_PECA = data.frame(Accession = acc, pepNum = pepNum, slr = slr, p = p, p.fdr.comb = p.fdr, q = q)
    if (!madePECA){
      PECA = add_PECA
      madePECA = T
    } else {
      PECA = rbind.fill(PECA, add_PECA)
    }
  }
  setTxtProgressBar(pb, round(which(accs == acc)/length(accs), 2))
}
PECA$extent = abs(PECA$slr - mean(PECA$slr))/sd(PECA$slr)
PECA[PECA$q > 0.05 | PECA$pepNum < 2,]$extent = 0


PECA = merge(MasterData, PECA, by = 'Accession', all.y = T)
plot(x,y, main = 'Significant Proteins vs. cutoff', xlab = 'Cutoff',ylab = '# of Significant Proteins')
dim(PECA[PECA$extent > 0, ])

#Write PECA results as a text file
write_tsv(PECA, 'PECA_HGSC_vs_CCC_full.txt')
### VOLCANO PLOT FROM PECA DATA
peca=PECA
lowest = min(peca[peca$q > 0, ]$q)
#Replace q values of 0 with a very small number to prevent infinity from messing things up
peca[peca$p.fdr == 0, ]$q = lowest/10
peca$adj.p.fdr=-log10(peca$q) # calculate the -log10(p.fdr) values (y-axis)
peca$abs.score = abs(peca$slr) * peca$adj.p.fdr

#only colour proteins that have a high abs.score
cutoff = mean(peca$abs.score) + 1*sd(peca$abs.score)
#peca=as.data.frame(subset(peca,p.fdr<0.05)) # you can remove proteins with p.fdr > 0.05 if you would like
yMax=max(peca$adj.p.fdr) # the largest value in the adjusted p.fdr column

# set the colour parameters for the volcano plot. Red = significantly up-regulated in high risk, 
#Blue = significantly down-regulated in high risk, Grey = not significantly changed
colour_parameters <- ifelse((peca$abs.score < cutoff | peca$q > 0.05), 'gray22', 
                            ifelse(peca$slr > 0,'red4',
                                   ifelse(peca$slr < 0,'blue','gray22')))
colour_opacity <- adjustcolor(colour_parameters, alpha.f = 0.5) # adjust the opacity of the points to 50%

#make curved lines showing where the cutoff is
xlim = c(min(peca$slr), max(peca$slr))
ylim = c(0,yMax*1.1)

y_points = seq(-log10(0.05), max(ylim)*1.1, length.out = 100)
x_points = cutoff/y_points
smoothingLine = smooth.spline(x_points, y_points, spar = 0.1)
smoothingLine2 = smooth.spline(-x_points, y_points, spar = 0.1)

pdf('VolcanoPlot_HGSC_vs_CCC_full.pdf')
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
