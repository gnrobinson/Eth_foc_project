#############################################################
# Title: Calculating coverage diversity (H) for all groups
# Author: Guy Robinson
#############################################################

suppressMessages(library(tidyverse))
suppressMessages(library(optparse))

###### Example input files ######
#FO0090 FO0001	FO0028	FO0046	FO0049	FO0061
#1  0 1	0	0	0
#1	1	1	1	1	1
#1	1	1	1	1	1
#1	1	1	1	0	0
#1	0	0	1	1	1
###############

#reading input files with Optparser
option_list = list(
  make_option(c("-a","--populationa"), type="character", default=NULL, 
              help="Pangenome matrix file of population 1 (unzipped) with isolate header", metavar="character"),
  make_option(c("-b","--populationb"), type="character", default=NULL, 
              help="Pangenome matrix file of population 2 (unzipped) with isolate header", metavar="character")
)

opt = parse_args(OptionParser(option_list = option_list))

tbl <- read_delim(file = opt$populationa, delim = "\t", col_names = TRUE, show_col_type = FALSE)
tbl0_1 <- tbl[1:length(tbl)]

tbl <- read_delim(file = opt$populationb, delim = "\t", col_names = TRUE, show_col_type = FALSE)
tbl0_2 <- tbl[1:length(tbl)]

# Calculating Total Haplotype Diversity (Hs) for first population
tbl1 <- as.data.frame(rowSums(tbl0_1 == "0"))
colnames(tbl1)[1] <- "zero"
tbl1$one <- rowSums(tbl0_1 == "1")
tbl1$total <- tbl1$zero + tbl1$one

tbl1$zero_sq <- (tbl1$zero / tbl1$total)^2
tbl1$one_sq <- (tbl1$one / tbl1$total)^2

tbl1$total_sq <- tbl1$zero_sq + tbl1$one_sq

Hs_1 <- 1-mean(tbl1$total_sq, na.rm = TRUE)

# Calculating Total Haplotype Diversity (Hs) for second population
tbl2 <- as.data.frame(rowSums(tbl0_2 == "0"))
colnames(tbl2)[1] <- "zero"
tbl2$one <- rowSums(tbl0_2 == "1")
tbl2$total <- tbl2$zero + tbl2$one

tbl2$zero_sq <- (tbl2$zero / tbl2$total)^2
tbl2$one_sq <- (tbl2$one / tbl2$total)^2

tbl2$total_sq <- tbl2$zero_sq + tbl2$one_sq

Hs_2 <- 1-mean(tbl2$total_sq, na.rm = TRUE)

#Calculating total haplotype diversity (Ht)
tbl_tot <- cbind(tbl0_1, tbl0_2)

tbl_tot2 <- as.data.frame(rowSums(tbl_tot == "0"))
colnames(tbl_tot2)[1] <- "zero"
tbl_tot2$one <- rowSums(tbl_tot == "1")
tbl_tot2$total <- tbl_tot2$zero + tbl_tot2$one

tbl_tot2$zero_sq <- (tbl_tot2$zero / tbl_tot2$total)^2
tbl_tot2$one_sq <- (tbl_tot2$one / tbl_tot2$total)^2

tbl_tot2$total_sq <- tbl_tot2$zero_sq + tbl_tot2$one_sq

Ht <- 1-mean(tbl_tot2$total_sq, na.rm = TRUE)

#calculate:
# Fst with equation (Ht - Hs) / Ht = (Ht - (Hs_1+ Hs_2)/2) / Ht
# Weighted Fst with equation (Ht - (Hs_1*Ns_1 + Hs_2*Ns_2) / (Ns_1 + Ns_2)) / Ht

fst <- (Ht - (Hs_1 + Hs_2) / 2) / Ht
fst_weighted <- (Ht - (Hs_1*ncol(tbl0_1) + Hs_2*ncol(tbl0_2)) / (ncol(tbl0_1)+ncol(tbl0_2))) / Ht

print(paste("Mean Fst for", opt$populationa, "&", opt$populationb, ":", fst))
print(paste("Weighted Fst for", opt$populationa, "&", opt$populationb, ":", fst_weighted))
