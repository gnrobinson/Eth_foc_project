
#loading libraries
library(zoo)
library(rbin)
library(ggplot2)
library(tidyverse)

###### Example input file ######
#  contig           start_bp final_bp       r2  diff
#  <chr>               <dbl>    <dbl>    <dbl> <dbl>
# 1 Supercontig_2.15      216     3219 0.0148    3003
# 2 Supercontig_2.15      216     3228 0.000905  3012
# 3 Supercontig_2.15      216     3949 0.00148   3733
# 4 Supercontig_2.15      216     4220 0.000866  4004
# 5 Supercontig_2.15      216     4307 1         4091
# 6 Supercontig_2.15      216     5080 0.00304   4864
##### This file was originally formatted from the output of VCFtools/PLINK --r2 ####

# ANCESTRALLY-DERIVED SNPS #

# INTRACHROMROSOMAL DATA #

#Loading LD files
ld.geno <- read_table2("/home/guy/Documents/Research/Cook/fusarium_pnas_final/133_isolates/clonality_metrics/r-squared/eth_clonecorrected.vcf_masked.vcf.60000.ld.diff", 
     col_names = FALSE, col_types = cols(X4 = col_number()))
colnames(ld.geno)[c(1,2,3,4,5)] = c("contig", "start_bp", "final_bp", "r2", "diff")
ld.geno$r2 <- as.numeric(ld.geno$r2)
ld.geno_plot <- subset(ld.geno, r2 <= 1)

binning_ld <- function(infile, min_diff, max_diff) {
  out <- subset(infile, diff >= min_diff & diff <= max_diff)
  return(out)
}

#binning for plotting
ld_geno0 <- binning_ld(ld.geno_plot, 0, 500)
ld_geno1 <- binning_ld(ld.geno_plot, 500, 1000)
ld_geno2 <- binning_ld(ld.geno_plot, 1000, 2000)
ld_geno3 <- binning_ld(ld.geno_plot, 2000, 3000)
ld_geno4 <- binning_ld(ld.geno_plot, 3000, 4000)
ld_geno5 <- binning_ld(ld.geno_plot, 4000, 5000)
ld_geno6 <- binning_ld(ld.geno_plot, 5000, 6000)
ld_geno7 <- binning_ld(ld.geno_plot, 6000, 7000)
ld_geno8 <- binning_ld(ld.geno_plot, 7000, 8000)
ld_geno9 <- binning_ld(ld.geno_plot, 8000, 9000)
ld_geno10 <- binning_ld(ld.geno_plot, 9000, 10000)
ld_geno11 <- binning_ld(ld.geno_plot, 10000, 11000)
ld_geno12 <- binning_ld(ld.geno_plot, 11000, 12000)
ld_geno13 <- binning_ld(ld.geno_plot, 12000, 13000)
ld_geno14 <- binning_ld(ld.geno_plot, 13000, 14000)
ld_geno15 <- binning_ld(ld.geno_plot, 14000, 15000)
ld_geno16 <- binning_ld(ld.geno_plot, 15000, 16000)
ld_geno17 <- binning_ld(ld.geno_plot, 16000, 17000)
ld_geno18 <- binning_ld(ld.geno_plot, 17000, 18000)
ld_geno19 <- binning_ld(ld.geno_plot, 18000, 19000)
ld_geno20 <- binning_ld(ld.geno_plot, 19000, 20000)
ld_geno21 <- binning_ld(ld.geno_plot, 20000, 25000)
ld_geno22 <- binning_ld(ld.geno_plot, 25000, 30000)
ld_geno23 <- binning_ld(ld.geno_plot, 30000, 35000)
ld_geno24 <- binning_ld(ld.geno_plot, 35000, 40000)
ld_geno25 <- binning_ld(ld.geno_plot, 40000, 45000)
ld_geno26 <- binning_ld(ld.geno_plot, 45000, 50000)
ld_geno27 <- binning_ld(ld.geno_plot, 50000, 55000)
ld_geno28 <- binning_ld(ld.geno_plot, 55000, 60000)

mean = c(sum(ld_geno0$r2)/length(ld_geno0$r2),
        sum(ld_geno1$r2)/length(ld_geno1$r2),
         sum(ld_geno2$r2)/length(ld_geno2$r2),
         sum(ld_geno3$r2)/length(ld_geno3$r2),
         sum(ld_geno4$r2)/length(ld_geno4$r2),
         sum(ld_geno5$r2)/length(ld_geno5$r2),
         sum(ld_geno6$r2)/length(ld_geno6$r2),
         sum(ld_geno7$r2)/length(ld_geno7$r2),
         sum(ld_geno8$r2)/length(ld_geno8$r2),
         sum(ld_geno9$r2)/length(ld_geno9$r2),
        sum(ld_geno10$r2)/length(ld_geno10$r2),
         sum(ld_geno11$r2)/length(ld_geno11$r2),
         sum(ld_geno12$r2)/length(ld_geno12$r2),
         sum(ld_geno13$r2)/length(ld_geno13$r2),
         sum(ld_geno14$r2)/length(ld_geno14$r2),
         sum(ld_geno15$r2)/length(ld_geno15$r2),
         sum(ld_geno16$r2)/length(ld_geno16$r2),
         sum(ld_geno17$r2)/length(ld_geno17$r2),
         sum(ld_geno18$r2)/length(ld_geno18$r2),
         sum(ld_geno19$r2)/length(ld_geno19$r2),
         sum(ld_geno20$r2)/length(ld_geno20$r2),
         sum(ld_geno21$r2)/length(ld_geno21$r2),
         sum(ld_geno22$r2)/length(ld_geno22$r2),
         sum(ld_geno23$r2)/length(ld_geno23$r2),
         sum(ld_geno24$r2)/length(ld_geno24$r2),
         sum(ld_geno25$r2)/length(ld_geno25$r2),
        sum(ld_geno26$r2)/length(ld_geno26$r2), 
        sum(ld_geno27$r2)/length(ld_geno27$r2),
        sum(ld_geno28$r2)/length(ld_geno28$r2))

sd <- c(sqrt(sum((ld_geno0$r2-mean(ld_geno0$r2))^2/(length(ld_geno0$r2)-1))),
        sqrt(sum((ld_geno1$r2-mean(ld_geno1$r2))^2/(length(ld_geno1$r2)-1))),
        sqrt(sum((ld_geno2$r2-mean(ld_geno2$r2))^2/(length(ld_geno2$r2)-1))),
        sqrt(sum((ld_geno3$r2-mean(ld_geno3$r2))^2/(length(ld_geno3$r2)-1))),
        sqrt(sum((ld_geno4$r2-mean(ld_geno4$r2))^2/(length(ld_geno4$r2)-1))),
        sqrt(sum((ld_geno5$r2-mean(ld_geno5$r2))^2/(length(ld_geno5$r2)-1))),
        sqrt(sum((ld_geno6$r2-mean(ld_geno6$r2))^2/(length(ld_geno6$r2)-1))),
        sqrt(sum((ld_geno7$r2-mean(ld_geno7$r2))^2/(length(ld_geno7$r2)-1))),
        sqrt(sum((ld_geno8$r2-mean(ld_geno8$r2))^2/(length(ld_geno8$r2)-1))),
        sqrt(sum((ld_geno9$r2-mean(ld_geno9$r2))^2/(length(ld_geno9$r2)-1))),
        sqrt(sum((ld_geno10$r2-mean(ld_geno10$r2))^2/(length(ld_geno10$r2)-1))),
        sqrt(sum((ld_geno11$r2-mean(ld_geno11$r2))^2/(length(ld_geno11$r2)-1))),
        sqrt(sum((ld_geno12$r2-mean(ld_geno12$r2))^2/(length(ld_geno12$r2)-1))),
        sqrt(sum((ld_geno13$r2-mean(ld_geno13$r2))^2/(length(ld_geno13$r2)-1))),
        sqrt(sum((ld_geno14$r2-mean(ld_geno14$r2))^2/(length(ld_geno14$r2)-1))),
        sqrt(sum((ld_geno15$r2-mean(ld_geno15$r2))^2/(length(ld_geno15$r2)-1))),
        sqrt(sum((ld_geno16$r2-mean(ld_geno16$r2))^2/(length(ld_geno16$r2)-1))),
        sqrt(sum((ld_geno17$r2-mean(ld_geno17$r2))^2/(length(ld_geno17$r2)-1))),
        sqrt(sum((ld_geno18$r2-mean(ld_geno18$r2))^2/(length(ld_geno18$r2)-1))),
        sqrt(sum((ld_geno19$r2-mean(ld_geno19$r2))^2/(length(ld_geno19$r2)-1))),
        sqrt(sum((ld_geno20$r2-mean(ld_geno20$r2))^2/(length(ld_geno20$r2)-1))),
        sqrt(sum((ld_geno21$r2-mean(ld_geno21$r2))^2/(length(ld_geno21$r2)-1))),
        sqrt(sum((ld_geno22$r2-mean(ld_geno22$r2))^2/(length(ld_geno22$r2)-1))),
        sqrt(sum((ld_geno23$r2-mean(ld_geno23$r2))^2/(length(ld_geno23$r2)-1))),
        sqrt(sum((ld_geno24$r2-mean(ld_geno24$r2))^2/(length(ld_geno24$r2)-1))),
        sqrt(sum((ld_geno25$r2-mean(ld_geno25$r2))^2/(length(ld_geno25$r2)-1))),
        sqrt(sum((ld_geno26$r2-mean(ld_geno26$r2))^2/(length(ld_geno26$r2)-1))),
        sqrt(sum((ld_geno27$r2-mean(ld_geno27$r2))^2/(length(ld_geno27$r2)-1))),
        sqrt(sum((ld_geno28$r2-mean(ld_geno28$r2))^2/(length(ld_geno28$r2)-1))))

ci95 <- c(1.9560*(sqrt(sum((ld_geno0$r2-mean(ld_geno0$r2))^2/(length(ld_geno0$r2)-1))/sqrt(length(ld_geno0$r2)-1))),
        1.9560*((sqrt(sum((ld_geno1$r2-mean(ld_geno1$r2))^2/(length(ld_geno1$r2)-1))/sqrt(length(ld_geno1$r2)-1)))),
        1.9560*((sqrt(sum((ld_geno2$r2-mean(ld_geno2$r2))^2/(length(ld_geno2$r2)-1))/sqrt(length(ld_geno2$r2)-1)))),
        1.9560*((sqrt(sum((ld_geno3$r2-mean(ld_geno3$r2))^2/(length(ld_geno3$r2)-1))/sqrt(length(ld_geno3$r2)-1)))),
        1.9560*((sqrt(sum((ld_geno4$r2-mean(ld_geno4$r2))^2/(length(ld_geno4$r2)-1))/sqrt(length(ld_geno4$r2)-1)))),
        1.9560*((sqrt(sum((ld_geno5$r2-mean(ld_geno5$r2))^2/(length(ld_geno5$r2)-1))/sqrt(length(ld_geno5$r2)-1)))),
        1.9560*((sqrt(sum((ld_geno6$r2-mean(ld_geno6$r2))^2/(length(ld_geno6$r2)-1))/sqrt(length(ld_geno6$r2)-1)))),
        1.9560*((sqrt(sum((ld_geno7$r2-mean(ld_geno7$r2))^2/(length(ld_geno7$r2)-1))/sqrt(length(ld_geno7$r2)-1)))),
        1.9560*((sqrt(sum((ld_geno8$r2-mean(ld_geno8$r2))^2/(length(ld_geno8$r2)-1))/sqrt(length(ld_geno8$r2)-1)))),
        1.9560*((sqrt(sum((ld_geno9$r2-mean(ld_geno9$r2))^2/(length(ld_geno9$r2)-1))/sqrt(length(ld_geno9$r2)-1)))),
        1.9560*((sqrt(sum((ld_geno10$r2-mean(ld_geno10$r2))^2/(length(ld_geno10$r2)-1))/sqrt(length(ld_geno10$r2)-1)))),
        1.9560*((sqrt(sum((ld_geno11$r2-mean(ld_geno11$r2))^2/(length(ld_geno11$r2)-1))/sqrt(length(ld_geno11$r2)-1)))),
        1.9560*((sqrt(sum((ld_geno12$r2-mean(ld_geno12$r2))^2/(length(ld_geno12$r2)-1))/sqrt(length(ld_geno12$r2)-1)))),
        1.9560*((sqrt(sum((ld_geno13$r2-mean(ld_geno13$r2))^2/(length(ld_geno13$r2)-1))/sqrt(length(ld_geno13$r2)-1)))),
        1.9560*((sqrt(sum((ld_geno14$r2-mean(ld_geno14$r2))^2/(length(ld_geno14$r2)-1))/sqrt(length(ld_geno14$r2)-1)))),
        1.9560*((sqrt(sum((ld_geno15$r2-mean(ld_geno15$r2))^2/(length(ld_geno15$r2)-1))/sqrt(length(ld_geno15$r2)-1)))),
        1.9560*((sqrt(sum((ld_geno16$r2-mean(ld_geno16$r2))^2/(length(ld_geno16$r2)-1))/sqrt(length(ld_geno16$r2)-1)))),
        1.9560*((sqrt(sum((ld_geno17$r2-mean(ld_geno17$r2))^2/(length(ld_geno17$r2)-1))/sqrt(length(ld_geno17$r2)-1)))),
        1.9560*((sqrt(sum((ld_geno18$r2-mean(ld_geno18$r2))^2/(length(ld_geno18$r2)-1))/sqrt(length(ld_geno18$r2)-1)))),
        1.9560*((sqrt(sum((ld_geno19$r2-mean(ld_geno19$r2))^2/(length(ld_geno19$r2)-1))/sqrt(length(ld_geno19$r2)-1)))),
        1.9560*((sqrt(sum((ld_geno20$r2-mean(ld_geno20$r2))^2/(length(ld_geno20$r2)-1))/sqrt(length(ld_geno20$r2)-1)))),
        1.9560*((sqrt(sum((ld_geno21$r2-mean(ld_geno21$r2))^2/(length(ld_geno21$r2)-1))/sqrt(length(ld_geno21$r2)-1)))),
        1.9560*((sqrt(sum((ld_geno22$r2-mean(ld_geno22$r2))^2/(length(ld_geno22$r2)-1))/sqrt(length(ld_geno22$r2)-1)))),
        1.9560*((sqrt(sum((ld_geno23$r2-mean(ld_geno23$r2))^2/(length(ld_geno23$r2)-1))/sqrt(length(ld_geno23$r2)-1)))),
        1.9560*((sqrt(sum((ld_geno24$r2-mean(ld_geno24$r2))^2/(length(ld_geno24$r2)-1))/sqrt(length(ld_geno24$r2)-1)))),
        1.9560*((sqrt(sum((ld_geno25$r2-mean(ld_geno25$r2))^2/(length(ld_geno25$r2)-1))/sqrt(length(ld_geno25$r2)-1)))),
        1.9560*((sqrt(sum((ld_geno26$r2-mean(ld_geno26$r2))^2/(length(ld_geno26$r2)-1))/sqrt(length(ld_geno26$r2)-1)))),
        1.9560*((sqrt(sum((ld_geno27$r2-mean(ld_geno27$r2))^2/(length(ld_geno27$r2)-1))/sqrt(length(ld_geno27$r2)-1)))),
        1.9560*((sqrt(sum((ld_geno28$r2-mean(ld_geno28$r2))^2/(length(ld_geno28$r2)-1))/sqrt(length(ld_geno28$r2)-1)))))

bp = c(250, 750, 1500, 2500, 3500, 4500, 5500, 6500, 7500, 8500, 
       9500, 10500, 11500, 12500, 13500, 14500, 15500, 16500, 17500, 18500, 
       19500, 22500, 27500, 32500, 37500, 42500, 47500, 52500, 57500)

average = data.frame(bp, mean, sd, ci95)

# INTERCHROMOSOMAL DATA #

ld.inter <- read_table2("/home/guy/Documents/Research/Cook/fusarium_pnas_final/133_isolates/clonality_metrics/r-squared/eth_clonecorrected.vcf_masked.vcf_all_inter.ld.diff")
colnames(ld.inter)[c(1,2,3,4,5)] = c("contig", "start_bp", "final_bp", "r2", "diff")
ld.inter_plot <- subset(ld.inter, r2 <= 1)

#############################################################################################################

# CLONALLY-DERIVED SNPS #
# INTRACHROMOSOMAL DATA (ONLY) #

# load data
ld.clones <- read_table2("/media/guy/Guy HardDrive/fusarium_ld/final/var_clones/diff_calculated/var_clones.ld", 
     col_types = cols(R2 = col_number()))
colnames(ld.clones)[c(1,2,3,4,5)] = c("contig", "start_bp", "final_bp", "r2", "diff")
ld.clones$r2 <- as.numeric(ld.clones$r2)
ld.clones_plot <- subset(ld.clones, r2 <= 1)

#binning for plotting
binning_ld <- function(infile, min_diff, max_diff) {
  out <- subset(infile, diff >= min_diff & diff <= max_diff)
  return(out)
}

ld_clone0 <- binning_ld(ld.clones_plot, 0, 500)
ld_clone1 <- binning_ld(ld.clones_plot, 500, 1000)
ld_clone2 <- binning_ld(ld.clones_plot, 1000, 2000)
ld_clone3 <- binning_ld(ld.clones_plot, 2000, 3000)
ld_clone4 <- binning_ld(ld.clones_plot, 3000, 4000)
ld_clone5 <- binning_ld(ld.clones_plot, 4000, 5000)
ld_clone6 <- binning_ld(ld.clones_plot, 5000, 6000)
ld_clone7 <- binning_ld(ld.clones_plot, 6000, 7000)
ld_clone8 <- binning_ld(ld.clones_plot, 7000, 8000)
ld_clone9 <- binning_ld(ld.clones_plot, 8000, 9000)
ld_clone10 <- binning_ld(ld.clones_plot, 9000, 10000)
ld_clone11 <- binning_ld(ld.clones_plot, 10000, 11000)
ld_clone12 <- binning_ld(ld.clones_plot, 11000, 12000)
ld_clone13 <- binning_ld(ld.clones_plot, 12000, 13000)
ld_clone14 <- binning_ld(ld.clones_plot, 13000, 14000)
ld_clone15 <- binning_ld(ld.clones_plot, 14000, 15000)
ld_clone16 <- binning_ld(ld.clones_plot, 15000, 16000)
ld_clone17 <- binning_ld(ld.clones_plot, 16000, 17000)
ld_clone18 <- binning_ld(ld.clones_plot, 17000, 18000)
ld_clone19 <- binning_ld(ld.clones_plot, 18000, 19000)
ld_clone20 <- binning_ld(ld.clones_plot, 19000, 20000)
ld_clone21 <- binning_ld(ld.clones_plot, 20000, 25000)
ld_clone22 <- binning_ld(ld.clones_plot, 25000, 30000)
ld_clone23 <- binning_ld(ld.clones_plot, 30000, 35000)
ld_clone24 <- binning_ld(ld.clones_plot, 35000, 40000)
ld_clone25 <- binning_ld(ld.clones_plot, 40000, 45000)
ld_clone26 <- binning_ld(ld.clones_plot, 45000, 50000)
ld_clone27 <- binning_ld(ld.clones_plot, 50000, 55000)
ld_clone28 <- binning_ld(ld.clones_plot, 55000, 60000)

mean_clone = c(sum(ld_clone0$r2)/length(ld_clone0$r2),
        sum(ld_clone1$r2)/length(ld_clone1$r2),
         sum(ld_clone2$r2)/length(ld_clone2$r2),
         sum(ld_clone3$r2)/length(ld_clone3$r2),
         sum(ld_clone4$r2)/length(ld_clone4$r2),
         sum(ld_clone5$r2)/length(ld_clone5$r2),
         sum(ld_clone6$r2)/length(ld_clone6$r2),
         sum(ld_clone7$r2)/length(ld_clone7$r2),
         sum(ld_clone8$r2)/length(ld_clone8$r2),
         sum(ld_clone9$r2)/length(ld_clone9$r2),
        sum(ld_clone10$r2)/length(ld_clone10$r2),
         sum(ld_clone11$r2)/length(ld_clone11$r2),
         sum(ld_clone12$r2)/length(ld_clone12$r2),
         sum(ld_clone13$r2)/length(ld_clone13$r2),
         sum(ld_clone14$r2)/length(ld_clone14$r2),
         sum(ld_clone15$r2)/length(ld_clone15$r2),
         sum(ld_clone16$r2)/length(ld_clone16$r2),
         sum(ld_clone17$r2)/length(ld_clone17$r2),
         sum(ld_clone18$r2)/length(ld_clone18$r2),
         sum(ld_clone19$r2)/length(ld_clone19$r2),
         sum(ld_clone20$r2)/length(ld_clone20$r2),
         sum(ld_clone21$r2)/length(ld_clone21$r2),
         sum(ld_clone22$r2)/length(ld_clone22$r2),
         sum(ld_clone23$r2)/length(ld_clone23$r2),
         sum(ld_clone24$r2)/length(ld_clone24$r2),
         sum(ld_clone25$r2)/length(ld_clone25$r2),
        sum(ld_clone26$r2)/length(ld_clone26$r2), 
        sum(ld_clone27$r2)/length(ld_clone27$r2),
        sum(ld_clone28$r2)/length(ld_clone28$r2))

sd_clone <- c(sqrt(sum((ld_clone0$r2-mean(ld_clone0$r2))^2/(length(ld_clone0$r2)-1))),
        sqrt(sum((ld_clone1$r2-mean(ld_clone1$r2))^2/(length(ld_clone1$r2)-1))),
        sqrt(sum((ld_clone2$r2-mean(ld_clone2$r2))^2/(length(ld_clone2$r2)-1))),
        sqrt(sum((ld_clone3$r2-mean(ld_clone3$r2))^2/(length(ld_clone3$r2)-1))),
        sqrt(sum((ld_clone4$r2-mean(ld_clone4$r2))^2/(length(ld_clone4$r2)-1))),
        sqrt(sum((ld_clone5$r2-mean(ld_clone5$r2))^2/(length(ld_clone5$r2)-1))),
        sqrt(sum((ld_clone6$r2-mean(ld_clone6$r2))^2/(length(ld_clone6$r2)-1))),
        sqrt(sum((ld_clone7$r2-mean(ld_clone7$r2))^2/(length(ld_clone7$r2)-1))),
        sqrt(sum((ld_clone8$r2-mean(ld_clone8$r2))^2/(length(ld_clone8$r2)-1))),
        sqrt(sum((ld_clone9$r2-mean(ld_clone9$r2))^2/(length(ld_clone9$r2)-1))),
        sqrt(sum((ld_clone10$r2-mean(ld_clone10$r2))^2/(length(ld_clone10$r2)-1))),
        sqrt(sum((ld_clone11$r2-mean(ld_clone11$r2))^2/(length(ld_clone11$r2)-1))),
        sqrt(sum((ld_clone12$r2-mean(ld_clone12$r2))^2/(length(ld_clone12$r2)-1))),
        sqrt(sum((ld_clone13$r2-mean(ld_clone13$r2))^2/(length(ld_clone13$r2)-1))),
        sqrt(sum((ld_clone14$r2-mean(ld_clone14$r2))^2/(length(ld_clone14$r2)-1))),
        sqrt(sum((ld_clone15$r2-mean(ld_clone15$r2))^2/(length(ld_clone15$r2)-1))),
        sqrt(sum((ld_clone16$r2-mean(ld_clone16$r2))^2/(length(ld_clone16$r2)-1))),
        sqrt(sum((ld_clone17$r2-mean(ld_clone17$r2))^2/(length(ld_clone17$r2)-1))),
        sqrt(sum((ld_clone18$r2-mean(ld_clone18$r2))^2/(length(ld_clone18$r2)-1))),
        sqrt(sum((ld_clone19$r2-mean(ld_clone19$r2))^2/(length(ld_clone19$r2)-1))),
        sqrt(sum((ld_clone20$r2-mean(ld_clone20$r2))^2/(length(ld_clone20$r2)-1))),
        sqrt(sum((ld_clone21$r2-mean(ld_clone21$r2))^2/(length(ld_clone21$r2)-1))),
        sqrt(sum((ld_clone22$r2-mean(ld_clone22$r2))^2/(length(ld_clone22$r2)-1))),
        sqrt(sum((ld_clone23$r2-mean(ld_clone23$r2))^2/(length(ld_clone23$r2)-1))),
        sqrt(sum((ld_clone24$r2-mean(ld_clone24$r2))^2/(length(ld_clone24$r2)-1))),
        sqrt(sum((ld_clone25$r2-mean(ld_clone25$r2))^2/(length(ld_clone25$r2)-1))),
        sqrt(sum((ld_clone26$r2-mean(ld_clone26$r2))^2/(length(ld_clone26$r2)-1))),
        sqrt(sum((ld_clone27$r2-mean(ld_clone27$r2))^2/(length(ld_clone27$r2)-1))),
        sqrt(sum((ld_clone28$r2-mean(ld_clone28$r2))^2/(length(ld_clone28$r2)-1))))

ci95_clone <- c(1.9560*(sqrt(sum((ld_clone0$r2-mean(ld_clone0$r2))^2/(length(ld_clone0$r2)-1))/sqrt(length(ld_clone0$r2)-1))),
        1.9560*((sqrt(sum((ld_clone1$r2-mean(ld_clone1$r2))^2/(length(ld_clone1$r2)-1))/sqrt(length(ld_clone1$r2)-1)))),
        1.9560*((sqrt(sum((ld_clone2$r2-mean(ld_clone2$r2))^2/(length(ld_clone2$r2)-1))/sqrt(length(ld_clone2$r2)-1)))),
        1.9560*((sqrt(sum((ld_clone3$r2-mean(ld_clone3$r2))^2/(length(ld_clone3$r2)-1))/sqrt(length(ld_clone3$r2)-1)))),
        1.9560*((sqrt(sum((ld_clone4$r2-mean(ld_clone4$r2))^2/(length(ld_clone4$r2)-1))/sqrt(length(ld_clone4$r2)-1)))),
        1.9560*((sqrt(sum((ld_clone5$r2-mean(ld_clone5$r2))^2/(length(ld_clone5$r2)-1))/sqrt(length(ld_clone5$r2)-1)))),
        1.9560*((sqrt(sum((ld_clone6$r2-mean(ld_clone6$r2))^2/(length(ld_clone6$r2)-1))/sqrt(length(ld_clone6$r2)-1)))),
        1.9560*((sqrt(sum((ld_clone7$r2-mean(ld_clone7$r2))^2/(length(ld_clone7$r2)-1))/sqrt(length(ld_clone7$r2)-1)))),
        1.9560*((sqrt(sum((ld_clone8$r2-mean(ld_clone8$r2))^2/(length(ld_clone8$r2)-1))/sqrt(length(ld_clone8$r2)-1)))),
        1.9560*((sqrt(sum((ld_clone9$r2-mean(ld_clone9$r2))^2/(length(ld_clone9$r2)-1))/sqrt(length(ld_clone9$r2)-1)))),
        1.9560*((sqrt(sum((ld_clone10$r2-mean(ld_clone10$r2))^2/(length(ld_clone10$r2)-1))/sqrt(length(ld_clone10$r2)-1)))),
        1.9560*((sqrt(sum((ld_clone11$r2-mean(ld_clone11$r2))^2/(length(ld_clone11$r2)-1))/sqrt(length(ld_clone11$r2)-1)))),
        1.9560*((sqrt(sum((ld_clone12$r2-mean(ld_clone12$r2))^2/(length(ld_clone12$r2)-1))/sqrt(length(ld_clone12$r2)-1)))),
        1.9560*((sqrt(sum((ld_clone13$r2-mean(ld_clone13$r2))^2/(length(ld_clone13$r2)-1))/sqrt(length(ld_clone13$r2)-1)))),
        1.9560*((sqrt(sum((ld_clone14$r2-mean(ld_clone14$r2))^2/(length(ld_clone14$r2)-1))/sqrt(length(ld_clone14$r2)-1)))),
        1.9560*((sqrt(sum((ld_clone15$r2-mean(ld_clone15$r2))^2/(length(ld_clone15$r2)-1))/sqrt(length(ld_clone15$r2)-1)))),
        1.9560*((sqrt(sum((ld_clone16$r2-mean(ld_clone16$r2))^2/(length(ld_clone16$r2)-1))/sqrt(length(ld_clone16$r2)-1)))),
        1.9560*((sqrt(sum((ld_clone17$r2-mean(ld_clone17$r2))^2/(length(ld_clone17$r2)-1))/sqrt(length(ld_clone17$r2)-1)))),
        1.9560*((sqrt(sum((ld_clone18$r2-mean(ld_clone18$r2))^2/(length(ld_clone18$r2)-1))/sqrt(length(ld_clone18$r2)-1)))),
        1.9560*((sqrt(sum((ld_clone19$r2-mean(ld_clone19$r2))^2/(length(ld_clone19$r2)-1))/sqrt(length(ld_clone19$r2)-1)))),
        1.9560*((sqrt(sum((ld_clone20$r2-mean(ld_clone20$r2))^2/(length(ld_clone20$r2)-1))/sqrt(length(ld_clone20$r2)-1)))),
        1.9560*((sqrt(sum((ld_clone21$r2-mean(ld_clone21$r2))^2/(length(ld_clone21$r2)-1))/sqrt(length(ld_clone21$r2)-1)))),
        1.9560*((sqrt(sum((ld_clone22$r2-mean(ld_clone22$r2))^2/(length(ld_clone22$r2)-1))/sqrt(length(ld_clone22$r2)-1)))),
        1.9560*((sqrt(sum((ld_clone23$r2-mean(ld_clone23$r2))^2/(length(ld_clone23$r2)-1))/sqrt(length(ld_clone23$r2)-1)))),
        1.9560*((sqrt(sum((ld_clone24$r2-mean(ld_clone24$r2))^2/(length(ld_clone24$r2)-1))/sqrt(length(ld_clone24$r2)-1)))),
        1.9560*((sqrt(sum((ld_clone25$r2-mean(ld_clone25$r2))^2/(length(ld_clone25$r2)-1))/sqrt(length(ld_clone25$r2)-1)))),
        1.9560*((sqrt(sum((ld_clone26$r2-mean(ld_clone26$r2))^2/(length(ld_clone26$r2)-1))/sqrt(length(ld_clone26$r2)-1)))),
        1.9560*((sqrt(sum((ld_clone27$r2-mean(ld_clone27$r2))^2/(length(ld_clone27$r2)-1))/sqrt(length(ld_clone27$r2)-1)))),
        1.9560*((sqrt(sum((ld_clone28$r2-mean(ld_clone28$r2))^2/(length(ld_clone28$r2)-1))/sqrt(length(ld_clone28$r2)-1)))))

bp = c(250, 750, 1500, 2500, 3500, 4500, 5500, 6500, 7500, 8500, 
       9500, 10500, 11500, 12500, 13500, 14500, 15500, 16500, 17500, 18500, 
       19500, 22500, 27500, 32500, 37500, 42500, 47500, 52500, 57500)

average_clone = data.frame(bp, mean_clone, sd_clone, ci95_clone)

# PLOTTING TOGETHER (FIGURE 4)

#plotting fixed vs variable within lineages
```{r}
average$group = "Total"
average_clone$type = "Clonal"
colnames(average) = c('bp', 'mean', 'sd', 'ci95', 'Group')
colnames(average_clone) = c('bp', 'mean', 'sd', 'ci95', 'Group')
average <- average[,1:5]
average_clone <- average_clone[,1:5]

ld_final_plot <- rbind(average, average_clone)

plot <- ggplot(ld_final_plot, aes(x = bp/1000, y = mean, color=Group)) +
  geom_point(size=1) + 
  geom_line() +
#  stat_smooth(method="loess", span = 0.4) +
  scale_y_continuous(breaks = c(seq(0,1,0.1)), limits = c(0,1)) +
  scale_x_continuous(breaks = c(seq(0,30,5)), limits = c(0,30)) +
  labs(y = bquote('Linkage Disequilibrium'~(r^2)), x = "Physical Distance (Kbp)")+
  theme(axis.text.x = element_text(size = 21),
        axis.text.y = element_text(size = 21),  
        axis.title.x = element_text(size = 21),
        axis.title.y = element_text(size = 21),
        panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.position = "none") +
  geom_errorbar(aes(ymin=mean-ci95, ymax=mean+ci95)) +
  geom_hline(yintercept = 0.083, linetype = 2) #0.083

plot
