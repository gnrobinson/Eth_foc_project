################################################################################
###### Title: Plotting pangenome curves from presence-absence matrix ###########
###### Author: Guy Robinson                                          ###########
################################################################################

#loading required libraries
library(tidyverse)
library(reshape2)
library(micropan)
library(ggplot2)

#read in data
#### Example input file 
#cluster  FO0090 FO0001 FO0028 FO0046 FO0049 FO0061 FO0066 FO0067 
#<chr>     <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl
#1 Cluster2      1      1      1      1      1      1      1
#2 Cluster3      1      1      1      1      1      1      1  
#3 Cluster4      1      1      1      1      1      1      1  
#4 Cluster5      1      1      1      1      1      1      1  
#5 Cluster6      1      1      1      1      1      1      1  
#6 Cluster7      1      1      1      1      1      1      1  

setwd("/home/guy/Documents/Research/Cook/fusarium_pnas_final/133_isolates/pangenome_plotting/")
tbl <- read.table("binary-data-conversion-selected-genomes-excluded.txt", header = TRUE, stringsAsFactors = FALSE)

finalpan <- as.data.frame(tbl[2:length(tbl)])# transpose data frame to be suitable for rarefaction. Warning: Removes row names too.

####################
cons_gene_count <- list()
acc_gene_count <- list()
for (k in seq(2, 99, 1)) {
  x <- list()
  make_permuts <- function(indf, a) {
    for (i in seq(1, 99, 1)) {
      Ps <- indf[,sample(ncol(indf), a)]
      x[[paste0("perm", i)]] <- Ps
    }
    return(x)
  }
  
  permuts_list <- make_permuts(finalpan, k)
  #make function that adds row sums to each dataframe within a list
  add_rowsums <- function(df) {
    lapply(df, function(y) {y$RowSum <- rowSums(y, na.rm=T); return(y)})
  }
  #make function that makes a list of data frames containing only conserved genes
  conserved_genes <- function(df) { 
    lapply(df, function(x) { x[ x$RowSum == k, ]})
  }
  #make function that makes a list of data frames containing only accessory genes
  accessory_genes <- function(df) {
    lapply(df, function(x) { x[ x$RowSum < k, ]})
  }

  permuts_list_rowsums <- add_rowsums(permuts_list) #Use function add_rowsums()
  
  cons_genes <- conserved_genes(permuts_list_rowsums) #Separate out number of conserved genes in presence-absence matrix
  cons_genes_sum <- sapply(cons_genes, nrow) #create vectors of number of conserved genes for all permutations
  cons_gene_count[[paste0(k)]] <- cons_genes_sum #create list of conserved genes at each genome number
  
  acc_genes <- accessory_genes(permuts_list_rowsums) #Separate out number of conserved genes in presence-absence matrix
  acc_genes_sum <- sapply(acc_genes, nrow) #create vectors of number of conserved genes for all permutations
  acc_gene_count[[paste0(k)]] <- acc_genes_sum #create list of conserved genes at each genome number
}

#calculating the accessory/core genome for genome = 1
sum(permuts_list$perm1$FO0090)
genome_1 <- sapply(permuts_list$perm1, sum)

# adding genomes
conserved_df <- as.data.frame(t(do.call(rbind, cons_gene_count)))
conserved_df2 <- cbind(genome_1, conserved_df)
colnames(conserved_df2)[1] <- '1'
conserved_df2$type <- "core"

accessory_df <- as.data.frame(t(do.call(rbind, acc_gene_count)))
accessory_df2 <- cbind(genome_1, accessory_df)
colnames(accessory_df2)[1] <- '1'
accessory_df2$type <- "accessory"

all_df <- rbind(conserved_df2, accessory_df2)

all_df_long <- melt(all_df, id = "type")
all_df_long$variable <- as.integer(all_df_long$variable)
ggplot(all_df_long, aes(x=variable, y=value)) +
  geom_point(aes(color=type)) +
  stat_smooth(aes(linetype=type), method = "loess", color='black', span = 0.1) +
  scale_y_continuous(limits = c(0, 35000)) +
  theme_classic() +
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        axis.text.x = element_text(size = 21),
        axis.text.y = element_text(size = 21),  
        axis.title.x = element_text(size = 21),
        axis.title.y = element_text(size = 21),) +
  xlab('Number of Genomes') +
  ylab("Number of Genes")
