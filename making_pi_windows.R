#####################################################################################################
#Title: "Making Nucleotide Diversity (Pi) Windows using VCFtools --window-pi"
#Author: Guy Robinson
#####################################################################################################

#loading libraries
library(tidyverse)
library(ggplot2)
library(rlist)

#loading nucleotide data
df.geno <- read_table2("$FILE",
                       col_names = TRUE)
colnames(df.geno)[c(1,2,3)] = c("contig", "start_bp", "pi")


#making 40Kb windows across all of the supercontigs
mkwindow_avg <- function(infile,midpoint){
  window_size <- subset(infile, start_bp > (midpoint-20000) & start_bp < (midpoint+20000))
  avg_pi <- aggregate(pi ~ contig, window_size, FUN=mean)
  avg_pi_n <- aggregate(pi ~ contig, window_size, FUN=length)
  avg_pi$n <- length(avg_pi_n)
  avg_pi$midpoint <- midpoint
  colnames(avg_pi)[c(1,2,3)] = c("contig", "avg_pi","n")
  return(avg_pi)
}

#Making genomic windows (40Kb in size in 10Kb increments)
```{r}
N <- max(df.geno$start_bp)
x <- list()
for (i in seq(10000, N, 10000)) {
  Ps <- mkwindow_avg(df.geno, i)
  x[[paste0("element", i)]] <- Ps
}
concatenated_df <- list.rbind(x)

#plotting the nucleotide diversity across the genome
prep_for_plotting <- function(infile, supercontig, chrom, start_contig_pos){
  subset_df <- subset(infile, contig == supercontig)
  subset_df$chrom <- chrom
  subset_df$start_contig_pos <- start_contig_pos
  subset_df$midpoint_pos <- subset_df$midpoint + subset_df$start_contig_pos
  return(subset_df)
}

supercontig_2.14_df <- prep_for_plotting(concatenated_df, "Supercontig_2.14", 1, 0)
supercontig_2.1_df <- prep_for_plotting(concatenated_df, "Supercontig_2.1", 1, 1603705)
supercontig_2.6_df <- prep_for_plotting(concatenated_df, "Supercontig_2.6", 2, 5961943)
supercontig_2.10_df <- prep_for_plotting(concatenated_df, "Supercontig_2.10", 2, 8618966)
supercontig_2.31_df <- prep_for_plotting(concatenated_df, "Supercontig_2.31", 2, 10863102)
supercontig_2.8_df <- prep_for_plotting(concatenated_df, "Supercontig_2.8", 4, 11537300)
supercontig_2.4_df <- prep_for_plotting(concatenated_df, "Supercontig_2.4", 4, 13865338)
supercontig_2.26_df <- prep_for_plotting(concatenated_df, "Supercontig_2.26", 5, 16749062)
supercontig_2.2_df <- prep_for_plotting(concatenated_df, "Supercontig_2.2", 5, 17628167)
supercontig_2.5_df <- prep_for_plotting(concatenated_df, "Supercontig_2.5", 7, 21662322)
supercontig_2.13_df <- prep_for_plotting(concatenated_df, "Supercontig_2.13", 7, 24350954)
supercontig_2.3_df <- prep_for_plotting(concatenated_df, "Supercontig_2.3", 8, 26008504)
supercontig_2.29_df <- prep_for_plotting(concatenated_df, "Supercontig_2.29", 8, 29315844)
supercontig_2.11_df <- prep_for_plotting(concatenated_df, "Supercontig_2.11", 9, 29991914)
supercontig_2.17_df <- prep_for_plotting(concatenated_df, "Supercontig_2.17", 9, 31968020)
supercontig_2.20_df <- prep_for_plotting(concatenated_df, "Supercontig_2.20", 10, 33295615)
supercontig_2.15_df <- prep_for_plotting(concatenated_df, "Supercontig_2.15", 10, 34449357)
supercontig_2.45_df <- prep_for_plotting(concatenated_df, "Supercontig_2.45", 10, 35997568)
supercontig_2.35_df <- prep_for_plotting(concatenated_df, "Supercontig_2.35", 11, 36190455)
supercontig_2.12_df <- prep_for_plotting(concatenated_df, "Supercontig_2.12", 11, 36650010)
supercontig_2.19_df <- prep_for_plotting(concatenated_df, "Supercontig_2.19", 12, 38526589)
supercontig_2.23_df <- prep_for_plotting(concatenated_df, "Supercontig_2.23", 12, 39821944)
supercontig_2.16_df <- prep_for_plotting(concatenated_df, "Supercontig_2.16", 13, 40758996)
supercontig_2.39_df <- prep_for_plotting(concatenated_df, "Supercontig_2.39", 13, 42136195)

pi_plot_df <- rbind(supercontig_2.14_df,
                    supercontig_2.1_df, 
                    supercontig_2.6_df, 
                    supercontig_2.10_df,
                    supercontig_2.31_df,
                    supercontig_2.8_df, 
                    supercontig_2.4_df, 
                    supercontig_2.26_df,
                    supercontig_2.2_df, 
                    supercontig_2.5_df, 
                    supercontig_2.13_df,
                    supercontig_2.3_df, 
                    supercontig_2.29_df,
                    supercontig_2.11_df,
                    supercontig_2.17_df,
                    supercontig_2.20_df,
                    supercontig_2.15_df,
                    supercontig_2.45_df,
                    supercontig_2.35_df,
                    supercontig_2.12_df,
                    supercontig_2.19_df,
                    supercontig_2.23_df,
                    supercontig_2.16_df,
                    supercontig_2.39_df)


pi_plot <- ggplot(pi_plot_df, aes(x=midpoint_pos, y=avg_pi)) +
  theme_minimal() +
  geom_rect(data=NULL,aes(xmin=0,xmax=5961943,ymin=-Inf,ymax=Inf),fill="lightblue") +   geom_rect(data=NULL,aes(xmin=5961943,xmax=11537300,ymin=-Inf,ymax=Inf),fill="white") +  geom_rect(data=NULL,aes(xmin=11537300,xmax=16749062,ymin=-Inf,ymax=Inf),fill="lightblue") +  geom_rect(data=NULL,aes(xmin=16749062,xmax=21662322,ymin=-Inf,ymax=Inf),fill="white") +  geom_rect(data=NULL,aes(xmin=21662322,xmax=26008504,ymin=-Inf,ymax=Inf),fill="lightblue") +  geom_rect(data=NULL,aes(xmin=26008504,xmax=29991914,ymin=-Inf,ymax=Inf),fill="white") +  geom_rect(data=NULL,aes(xmin=29991914,xmax=33295615,ymin=-Inf,ymax=Inf),fill="lightblue") +  geom_rect(data=NULL,aes(xmin=33295615,xmax=36190455,ymin=-Inf,ymax=Inf),fill="white") +  geom_rect(data=NULL,aes(xmin=36190455,xmax=38526589,ymin=-Inf,ymax=Inf),fill="lightblue") +  geom_rect(data=NULL,aes(xmin=38526589,xmax=40758996,ymin=-Inf,ymax=Inf),fill="white") +
  geom_rect(data=NULL,aes(xmin=40758996,xmax=Inf,ymin=-Inf,ymax=Inf),fill="lightblue") +
  geom_point(size=0.1, alpha=1) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.major.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.minor.y=element_blank()) +
  ylab(expression(pi))

pi_plot
