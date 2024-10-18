####This is the first R script for today##
##Before you begin, you need to sort the bioinfomatics workspace in your local machine##

###Remember the golden rules of R
##Hashtag the R script to prevent running the code or to add comments
##Always check if you have R libraries installed - this doesn't have to be done every time, just once.
##Also remember to load your libraries every time you enter into Rstudio
##Always, always check your file paths and places you save the files to
##Please feel free to play around with the script, and if you have any questions or would like to
##try new things - google R studio resources - there are plenty!

###First - remember to set your working directory. There are two ways to do this
# 1) Click on taskbar - session - set working directory - choose directory 
# 2) Set the path right here in R: e.g.
setwd("/Users/lachlanmiller/code/university/binf7001/Assignment4/Pracs/Prac1_gene_centric_metagenomics")
###or you can set the path yourself.

###Save the script to the "code" folder after you downloaded it from the blackboard.

# Create directories to store data. ##The script below will create directories for you - one for code, input_data, figures and results##
dir.create(file.path(".", "Results"), showWarnings = FALSE)
dir.create(file.path(".", "input_data"), showWarnings = FALSE)
dir.create(file.path(".", "Figures"), showWarnings = FALSE)

# Install and load  packages required for this study

# The BiocManager amd pacman packages are R package management tools. 
# They make it more straightforward to install certain packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("pacman"))
  install.packages("pacman")
if (!require("ComplexHeatmap", quietly = TRUE))
  BiocManager::install("ComplexHeatmap")

# Use the p_load function from the pacman package to check if packages are installed. 
# If not they are not installed, the function will download and and install them for you
pacman::p_load("tidyverse", "vegan", "ggplot2", "reshape2", "openxlsx", "ggnewscale",
               "colorspace", "viridis", "circlize")

###Some libraries need to be loaded separately.
library(scales)

###################################################################################

###2.3)In this script, we will work with read counts. The first code block below will help
#us compare the number of reads for each type - microbial and host. After that we will
#compare how the read counts differ for the superworm groups - Bran, PS and starvation.

###Read the readcount file - this file should be in the "Input_data" folder
read.count <- read.delim("input_data/readcount.txt", header = T)

###Melt the data into a longer format using Sample_ID, microbial and host read counts
reads <- melt(read.count[, c(1, 4, 6)])

##Plot using ggplot2. Please have a look at ggplot2 resources online if you would like to change variables.
read.count.plot <- ggplot(reads, aes(x = Sample_ID, y = value, fill = variable)) +
  geom_col(position = "dodge") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10 ^ x),
              labels = trans_format("log10", math_format(10 ^ .x))) +
  theme_bw() +
  ylab("Log-transformed read count") +
  guides(fill = guide_legend(title = "Read type"))

##Call the plot 
read.count.plot

##Save the plot
ggsave("Figures/Read_count.png", plot = read.count.plot, dpi = 600, units = c("in"), width = 8, height = 5)


##Plot ratio using ggplot2. 
read.ratio.plot <- ggplot(reads, aes(x = Sample_ID, y = value, fill = variable)) +
  geom_col(position = "fill") +
  theme_bw() +
  ylab("read count ratios") +
  guides(fill = guide_legend(title = "Read type")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_discrete(labels = function(x) gsub("_", " ", x))


##Call the plot 
read.ratio.plot
ggsave("Figures/Read_ratio-readtype.png", plot = read.ratio.plot, dpi = 600, units = c("in"), width = 4, height = 3)


###Melt the data into a longer format using sample_type, microbial and host read counts
samples <- melt(read.count[, c(2, 4, 6)])

read.count.box <- ggplot(samples, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(coef = 10) +
#scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
#              labels = trans_format("log10", math_format(10^.x))) +
theme_bw() +
# ylab("Log-transformed read count") +
ylab("read count") +
  guides(fill = guide_legend(title = "Read type"))

read.count.box

ggsave("Figures/Read_count-readtype.png", plot = read.count.box, dpi = 600, units = c("in"), width = 4, height = 5)

######################################################################################################
