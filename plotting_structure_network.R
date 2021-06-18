#####################################################################################################
#Title: "Combining STRUCTURE and Network Analysis to Assess Clonality"
#Author: Guy Robinson
#####################################################################################################
  
#loading required libraries for odata wrangling and calculations
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))
suppressMessages(library(tibble))
suppressMessages(library(reshape2))
suppressMessages(library(optparse))

#reading input files with Optparser
option_list = list(
  make_option(c("-l","--list"), type="character", default=NULL, 
              help="list of input files (with path if not in the current directory)", metavar="character"),
  make_option(c("-b","--label"), type="character", default=NULL, 
              help="list of isolate names (in order!)", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output name", metavar="character"),
  make_option(c("-k", "--title"), type="character", default=NULL, 
              help="title of graph", metavar="character")
)
opt = parse_args(OptionParser(option_list = option_list))
Labels <- read.table(opt$label)
file_list <- scan(opt$list, what=character(), sep="\n", quiet=TRUE)

##################################################################################################
########################### DATA MANIPULATION AND CALCULATIONS ###################################
##################################################################################################
#Make list of data frames from paths contained in file_list (--list flag)
my_data <- list()
for (i in seq_along(file_list)) {
  my_data[[i]] <- read.table(file = file_list[i], header = F)
}

#Function  wrangles data into correct format
data_prep <- function(data_list) {
  INPUT <- cbind(Labels[1], data_list)
  input <- INPUT[c(1,7:length(INPUT))]
  X <- suppressWarnings(dist(input, method = "euclidean"))
  X <- as.matrix(X, labels = T)
  colnames(X) = rownames(X) = input$V1
  K <- melt(X)
  return(K)
}

#apply data_prep function onto list of dataframes
X <- lapply(my_data, FUN = data_prep)

#Merge all dataframes of within list into a single dataframe
merged_data = suppressWarnings(Reduce(function(...) merge(..., all=T, by=c(1,2)), X))
colnames(merged_data) = c("Var1", "Var2", 3:length(merged_data))

#creating new dataframe containing statistics
merged_data_stats <- merged_data[1:2]
#calculate the mean
merged_data_stats$mean = rowSums(merged_data[3:length(merged_data)]) / (length(merged_data)-2)
#calculating the standard deviation
merged_data_stats$SD =  apply(merged_data[3:length(merged_data)], MARGIN = 1, FUN = sd)

#####################################################################################
######################### NETWORK ANALYSIS AND PLOTTING #############################
#####################################################################################

#loading libraries for network analysis
suppressMessages(library(data.tree))
suppressMessages(library(ggraph))
suppressMessages(library(igraph))
suppressMessages(library(RColorBrewer))
suppressMessages(library(graphlayouts))
suppressMessages(library(gridExtra))


#Plotting network analysis
K_grouping <- subset(merged_data_stats, mean - SD <= 0) #grouping based on 1 SD away from 0. Should be only 5% false negative, according to 68:95:99 rule.

nodelist <- Labels
nodelist$Size <- 2 #arbitrarily chose 2, can be changed though

edgelist <- K_grouping[, c(1:3)]
g <- graph_from_data_frame(edgelist, vertices = nodelist, directed = FALSE)


#Export image of network
tiff(file=paste0(opt$out,".tiff"),width = 2000, height = 2000,res=300)
ggraph(g, layout = "fr") +
    geom_edge_link0(width=0.08,colour="black") +
    geom_node_point(col="Black",size=1, position = 'identity') +
    geom_node_text( aes(label=name), size = 3, color='black', repel = TRUE) +
    theme_graph() +
  labs(title = paste(c(opt$title), collapse = " "))
dev.off()


#creating dataframe with all membership groups
K_memb <- graph.data.frame(K_grouping)
K_memb <- as.data.frame(clusters(K_memb)$membership)
colnames(K_memb) <- "Group"

#write out membership tables
write.table(K_memb, file = paste0(opt$out, ".txt"), sep = '\t', row.names = TRUE, quote = FALSE, col.names = NA)

#praise
suppressMessages(library(praise))
praise(template = "${Exclamation}! What a ${adjective} person you are!")
