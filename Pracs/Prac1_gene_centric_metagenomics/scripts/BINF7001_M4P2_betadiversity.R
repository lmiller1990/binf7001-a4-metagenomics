###2.5) Determine the most abundant OTUs using DESeq2.

###Install the packages and load the libraries
# install.packages("DESeq2")
# BiocManager::install("DESeq2")
#install.packages("ggforce")
#install.packages("pheatmap")
#install.packages("RColorBrewer")
library(DESeq2)
library(ggplot2)
library(ggforce)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

#Read the counts table - this has all the counts summed up at the genus level.Also, the 
#last row contains grand total of all counts, this needs to be removed.
genus.data <- read.csv("Input_data/genus_counts.csv", header = T, row.names = 1)
dim(genus.data)
genus.data <- genus.data[-nrow(genus.data),]
tail(genus.data)

###Load the metadata file
metadata <- read.delim("Input_data/metadata.txt", header = T, row.names = 1)
dim(metadata)

###Prepare a DESeq object by creating a count matrix input with the count data
#For this, metadata inputs need to be factorised before running the matrix
colData <- metadata

# It is absolutely critical that the columns of the count matrix and the rows 
# of the column data (information about samples) are in the same order. 
# DESeq2 will not make guesses as to which column of the count matrix belongs to 
# which row of the column data, these must be provided to DESeq2 already in consistent order.

all(rownames(colData) == colnames(genus.data))

###Create the count matrix with DESeqDataSetFromMatrix, and run deseq2
count.data.set <- DESeqDataSetFromMatrix(countData= genus.data, 
                                         colData=colData, design= ~Type)
count.data.set
dds <- DESeq(count.data.set)
dim(dds)
dds

####Keep those OTUs with hits >=1 since we have low number of OTUs.
keep <- rowSums(counts(dds)) >=1
dds <- dds[keep,]
dim(dds)

#Normalize with variance stabilizing transformation with blind dispersion
vsd <- varianceStabilizingTransformation(dds)
vsd
head(assay(vsd))
dim(vsd)
vsd_tax <- rownames(vsd)
vsd_norm <- cbind(vsd_tax, assay(vsd))
write.table(vsd_norm, file = "Results/normalized_data.csv", sep = ",", quote = F, row.names = F)

###################################################################
#Explore the variation in the community using Principal component analyses
###################################################################

#How to get PCA scree plot?
pca_res <- prcomp(t(assay(vsd)), scale. = TRUE)
## calculate the variance for each gene
rv <- rowVars(assay(vsd))
## select the ntop genes by variance
select <- order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
## perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(vsd)[select,]))
## the contribution to the total variance for each component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
##plot the "percentVar"
scree_plot=data.frame(percentVar)
scree_plot[,2]<- c(1:9)
colnames(scree_plot)<-c("variance","component_number")
ggplot(scree_plot, mapping=aes(x=component_number, y=variance))+ 
  geom_bar(stat="identity")

ggsave("Figures/scree_plot.png", dpi = 600, units = c("in"), width = 7, height = 5)

###########Plot the PCA plot#############
data.PCA <- plotPCA(vsd, intgroup=c("Type"), returnData=TRUE)
percentVar <- round(100 * attr(data.PCA, "percentVar"))
#data.PCA$Month <- factor(data.PCA$Month, levels = c("February", "May", "July", "August", "October", "December", "Feb-21"))

all_PCA <- ggplot(data.PCA, aes(PC1, PC2, color=Type)) +
  geom_point(size=2)+ 
  ggtitle("Principal Component Analysis - superworm data set") +
  theme(panel.background = element_blank(),
        axis.line = element_line()) +
  scale_color_manual(values = c("#0072B2", "#009E73", "#D55E00")) +
  theme(legend.position = "bottom") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  xlim(-30,30) +
  ylim(-15,15) +
  ggforce::geom_mark_ellipse(aes(color = data.PCA$Type), show.legend = NA)

all_PCA

ggsave("Figures/PCA_plot.png", plot = all_PCA, dpi = 600, units = c("in"), width = 5, height = 5)


##################################################################
######Dissimilarity matrix to get sample to sample distances#####
#################################################################

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$Type
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         filename = "Figures/dissimilar.png",
         col=colors)

#####################################################################
#Look at the most abundant members of the community using a heatmap
####################################################################


#####sort the abundance table in decreasing abundance, and pick top 50 abundant OTUs
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)
DF <- assay(vsd)[select,][1:50,] ##picks the top 50 OTUs
DF <- DF[order(row.names(DF)),] ##arranges them by order

###Plot the pretty heatmap!
pheatmap(DF, cluster_rows=FALSE, show_rownames=TRUE,
         fontsize_col = 12, 
         show_colnames = TRUE, cluster_cols=TRUE,
         annotation_legend = TRUE,main = "Top50 abundant OTUs", 
         legend = TRUE, 
         border_color = "NA" ,annotation_names_col = TRUE,
         gaps_col = c(3,6),
         #gaps_row = c(3,6,13,15,16,18,24,27,30,31,34,40,46,47,48,49),
         filename = "Figures/Top_50_genus-BlRd.png", cellwidth = 15,cellheight = 15,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdBu")))(500))

dev.off()

###Play around with row numbers if you're interested!

#######################################################################################