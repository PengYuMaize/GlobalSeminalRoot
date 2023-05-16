
# Part I, III, IV were conducted on desktop, and part II was conducted on the cluster


# Part I ------------------------------------------------------------------

library(genio)
library(dplyr)
library(tidyverse)
library(data.table)

#read phenotypic data to get the individual list
SRNMat <- read.csv('output/SRN_Predictions861Accessions_82322_Final.csv', header = TRUE, row.names = 'X', na.strings = '')
colnames(SRNMat)[2] <- 'SampleID'
SRNMat <- SRNMat %>% drop_na(SampleID) # keep SEEDs 1777

SRN<-SRNMat[,'SampleID'] 
names(SRN)<-SRNMat$SampleID

# read genotypic file;genotype for SRN (remove na)
colms<-lapply(1:10,function(chr){
    mat<-data.frame(fread(paste('genotype/Ch',chr,'Merged.hmp.txt',sep="")))
    nms<-sapply(mat[,1],function(x) strsplit(x,'.M')[[1]][1])
    dups<-duplicated(nms)
    xx<-mat[-which(dups==TRUE),]
    rownames(xx)<-sapply(xx[,1],function(x) strsplit(x,'.M')[[1]][1])
    mat<-xx
    mat<-mat[,-1] 
    mat<-mat[names(SRN),]
    mono<-apply(mat,MARGIN=2,function(x) length(table(x)))
    mat<-mat[,-which(mono==1)]
    mat<-t(as.matrix(mat))
    mat<-mat*2
    return(mat)
})

save(colms, file = 'cimmyt_genotype_list_SRN.Rimage')

# Part II -----------------------------------------------------------------

#unlist file
load('cimmyt_genotype_list_SRN.Rimage')

colm<-data.frame(rbind(colms[[1]],colms[[2]],colms[[3]],colms[[4]],colms[[5]],colms[[6]],colms[[7]],colms[[8]],colms[[9]],colms[[10]]))

save(colm, file = 'cimmyt_genotype_file.Rimage')
#library(devtools)

#install_version("vctrs", version = "0.3.8", repos = "http://cran.us.r-project.org", lib = '/storage/home/m/mul826/R/x86_64-redhat-linux-gnu-library/3.6')

#install.packages('genio', repos = "http://cran.us.r-project.org", lib = '/storage/home/m/mul826/R/x86_64-redhat-linux-gnu-library/3.6', INSTALL_opts = '--no-lock')

library(genio)
library(dplyr)
#library(tidyverse)

load('cimmyt_genotype_file.Rimage')
colm[colm ==3] <- NA

colm.m <- as.matrix(colm)

#get id for .fam file
id.SRN <- as.data.frame(colnames(colm)) %>% mutate(fam = 0)
colnames(id.SRN)[1] <- 'id'

#read marker position
marker <- read.csv('genotype/Dan_genotype_site_summary.csv', row.names = 'X')
#rename columns to create .bim file
names(marker) <- c('id','chr','pos')
rownames(marker) <- marker$id

inters<-as.data.frame(intersect(rownames(colm),marker$id))
colnames(inters) <-'inter'
rownames(inters) <- inters$inter

marker <- marker[rownames(inters),] # reorder markers as the same as colm file

# write this data to BED/BIM/FAM files
# output path without extension
file_out <- 'Peng_SRN_Cimmyt'

SRN_bim <- make_bim(marker)

SRN_fam <- make_fam(id.SRN)

write_plink(file_out, colm.m, bim = SRN_bim,fam = SRN_fam)

#-----part III: create a phenotype file -----
# create phenotype file for magma
id.pheno <- SRNMat %>% 
    #filter(SaID %in% cimmyt) %>%
    mutate(fam = 0) %>% 
    dplyr::select(fam, SampleID,SRN_pred)
colnames(id.pheno)[2] <- 'id'
write.table(id.pheno, file = 'cimmyt_pred_SRN_pheno.txt',quote = FALSE, row.names = FALSE )

#------part IV: create a covar file------------------
# create a covar file including population structure
#pop structure:
load('genotype/HighFilteredReduced.Rimage')
inters<-intersect(SRN,rownames(genoMDS)) #n=1614
genoMDS<-genoMDS[inters,]
pops<-cmdscale(dist(genoMDS),k=5)
pops_fam <- as.data.frame(pops) %>% rownames_to_column(var = 'id') %>% mutate(fam=0) %>% select(fam,id,V1:V5)

write.table(pops_fam, file = 'population_structure_n_1614_SRN.txt', quote = FALSE, row.names = FALSE)
