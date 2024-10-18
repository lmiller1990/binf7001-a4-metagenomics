####Section 2.4) Calculate alpha diversity in your samples.
###Install and load the libraries

#install.packages("phyloseq")
#install.packages("plyr")
#install.packages("tibble")
library(phyloseq)
library(plyr)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)

###Load the data. This data contains taxa counts for each sample. Explore this file.
#The last row of the file sums up the total count of taxa per sample - just to give you an idea.
data <- read.csv("input_data/species_counts.csv", header = T)
dim(data)
data <- data[-nrow(data),]
tail(data)
data <- as_tibble(data)

##Load the metadata file. Please ensure this file is stored in the Input_data folder.
metadata <- read.delim("input_data/metadata.txt", header = T)
dim(metadata)

#We will use the phyloseq package to create a wrapper before computing diversity indices.
#The way this works is that it stores taxa counts as counts, taxonomy as tax and metadata as sample IDs.
#This package makes it easy to work with and compute alpha diversity indices at once.

otu_mat<- data.frame(data)
colnames(otu_mat)[1] <- "Taxonomy"
tax_mat<- data.frame(otu_mat$Taxonomy)
colnames(tax_mat)[1] <- "Taxonomy"
tax <- separate(
  tax_mat,
  col = 1,
  into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
  sep = ";",
  remove = TRUE,
  convert = FALSE,
  extra = "warn",
  fill = "warn")
tax_mat <- cbind(tax_mat, tax)
samples_df <- as_tibble(metadata)
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("Taxonomy") 
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("Taxonomy")
samples_df <- samples_df %>% 
  tibble::column_to_rownames("Sample_ID") 

####Convert the dataframes to matrices#####
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
dim(otu_mat)
dim(tax_mat)

####Prepare the objects in phyloseq#####
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

####Run phyloseq#####
phy_seq <- phyloseq(OTU, TAX, samples)
phy_seq
sample_names(phy_seq)
rank_names(phy_seq)
sample_variables(phy_seq)
Sample_ID <- row.names(samples_df)
metadata <- cbind(Sample_ID, samples_df)

####Estimate species richness and evenness with phyloseq #######
#####Compute richness#####
species_richness<-estimate_richness(phy_seq, split=TRUE, measures=c("Shannon","Observed", "Simpson", "Chao"))
richness <- cbind(metadata, species_richness)
species_richness.df<-as.data.frame(species_richness)
species_richness.df<-tibble::rownames_to_column(as.data.frame(species_richness), var="Sample_ID")
species_richness_table<-right_join(metadata, species_richness.df, by="Sample_ID")

species_evenness<-as.data.frame((richness$Shannon)/(log(richness$Observed)))
evenness <- cbind(richness, species_evenness)
colnames(evenness)[8] <- "species-evenness"
str(evenness)

#plot Shannon index for all Hosts
pd<-position_dodge(0.3)

ShannonPlot <-ggplot(evenness, aes(x=Type, y=Shannon, color = Type))+
  geom_boxplot() +
  theme_bw()+
  theme(legend.position = "right",
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(),
        axis.title.y = element_text(colour = "black", size = 12),
  ) +
  scale_y_continuous()+
  scale_color_brewer(palette="Dark2") +
  ggtitle("Prokaryotic diversity")

ShannonPlot
ggsave("Figures/Shannonplot.png", plot = ShannonPlot, dpi = 600, units = c("in"), width = 5, height = 6)

SimpsonPlot <-ggplot(evenness, aes(x=Type, y=Simpson, color = Type))+
  geom_boxplot() +
  theme_bw()+
  theme(legend.position = "right",
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(),
        axis.title.y = element_text(colour = "black", size = 12),
  ) +
  scale_y_continuous()+
  scale_color_brewer(palette="Dark2") +
  ggtitle("Prokaryotic diversity")

SimpsonPlot
ggsave("Figures/Simpsonplot.png", plot = SimpsonPlot, dpi = 600, units = c("in"), width = 5, height = 6)

###Feel free to plot other diversity indices!
#########################################################################


# Combine the plots
CombinedPlot <- grid.arrange(ShannonPlot, SimpsonPlot, ncol=2)

# Save the combined plot
ggsave("Figures/CombinedPlot.png", plot = CombinedPlot, dpi = 600, units = c("in"), width = 10, height = 6)