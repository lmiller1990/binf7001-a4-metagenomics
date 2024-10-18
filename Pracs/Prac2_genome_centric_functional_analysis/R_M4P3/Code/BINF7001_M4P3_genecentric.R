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
#Click on taskbar - session - set working directory - choose directory and copy over your command from the console below, eg
setwd("/Users/lachlanmiller/code/university/binf7001/Assignment4/Pracs/Prac2_genome_centric_functional_analysis/R_M4P3/Code")
###or you can set the path yourself.

###Save the script to the "code" folder after you download it from the blackboard.

# Create directories to store data. ##The script below will create directories for you - one for code, input_data, figures and results##
dir.create(file.path("../", "Results"), showWarnings = FALSE)
#dir.create(file.path(".", "Input_data"), showWarnings = FALSE)
dir.create(file.path("../", "Figures"), showWarnings = FALSE)

# Install and load  packages required for this study

# The BiocManager amd pacman packages are R package management tools. 
# They make it more straightforward to install certain packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(!require("pacman"))
  install.packages("pacman")
if (!require("ComplexHeatmap", quietly = TRUE))
  BiocManager::install("ComplexHeatmap")
if (!require("KEGGREST", quietly = TRUE))
  BiocManager::install("KEGGREST")


# Use the p_load function from the pacman package to check if packages are installed. 
# If not they are not installed, the function will download and and install them for you
pacman::p_load("tidyverse", "vegan", "ggplot2", "reshape2", "openxlsx", "ggnewscale", 
               "colorspace", "viridis", "circlize")

###Some libraries need to be loaded seperately.
library(scales)

###################################################################
#
###################################################################
###3.2) Determine the most abundant KO using DESeq2.

###Install the packages and load the libraries
#install.packages("DESeq2")
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
KO.data <- read.delim("../Input_data/KO_superworms.txt", header = T, row.names = 1)
dim(KO.data)

###Load the metadata file
metadata <- read.delim("../Input_data/metadata.txt", header = T, row.names = 1)
dim(metadata)

###Prepare a DESeq object by creating a count matrix input with the count data
#For this, metadata inputs need to be factorised before running the matrix
colData <- metadata

# It is absolutely critical that the columns of the count matrix and the rows 
# of the column data (information about samples) are in the same order. 
# DESeq2 will not make guesses as to which column of the count matrix belongs to 
# which row of the column data, these must be provided to DESeq2 already in consistent order.

all(rownames(colData) == colnames(KO.data))

###Create the count matrix with DESeqDataSetFromMatrix, and run deseq2
count.data.set <- DESeqDataSetFromMatrix(countData= KO.data, 
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
vsd_kegg_ID <- rownames(vsd)
vsd_norm <- cbind(vsd_kegg_ID, assay(vsd))
write.table(vsd_norm, file = "../Results/normalized_data.csv", sep = ",", quote = F, row.names = F)

############################################################################
#3.3)Look at the variation in the community using Principal component analyses
############################################################################

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
scree_plot[,2]<- c(1:8)
colnames(scree_plot)<-c("variance","component_number")
ggplot(scree_plot, mapping=aes(x=component_number, y=variance))+ 
  geom_bar(stat="identity")

ggsave("../Figures/scree_plot.png", dpi = 600, units = c("in"), width = 7, height = 5)

###########Plot the PCA plot#############
data.PCA <- plotPCA(vsd, intgroup=c("Type"), returnData=TRUE)
percentVar <- round(100 * attr(data.PCA, "percentVar"))
data.PCA$Type<- factor(data.PCA$Type, levels = c("Bran" ,"Polystyrene", "No food" , "Faeces from Bran group", "Faeces from PS group"))

all_PCA <- ggplot(data.PCA, aes(PC1, PC2, color=Type)) +
  geom_point(size=3)+ 
  ggtitle("Principal Component Analysis - superworm data set") +
  theme(panel.background = element_blank(),
        axis.line = element_line()) +
  scale_color_manual(values = c("#0072B2", "#009E73", "#D55E00", "#A52A2A", "#00FF00")) +
  theme(legend.position = "right") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  xlim(-20,20) +
  ylim(-30,10)
  #ggforce::geom_mark_ellipse(aes(color = data.PCA$Type), show.legend = NA)

all_PCA

ggsave("../Figures/PCA_plot.png", plot = all_PCA, dpi = 600, units = c("in"), width = 5, height = 5)


######################################################################
######3.4) Dissimilarity matrix to get sample to sample distances#####
#####################################################################

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$Type
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

##################################################################
###############KEGGREST database in R ############################
##################################################################

###Individual pathway search from kegg databases.
#For getting a list of kegg ids (KO) , description and modules
pacman::p_load("tidyverse", "data.table", "KEGGREST","openxlsx")
library(utils)
library(KEGGREST)
library(raster)
# Function to convert a matrix to a dataframe. Rows become first column. 
# Specify 'column_name' to set the name of the new column (defaults to "Row_variable")
m2df <- function(mymatrix, column_name = "Row_variable"){
  mydf <- as.data.frame(mymatrix)
  cur_names <- names(mydf)
  mydf[, column_name] <- rownames(mydf)
  rownames(mydf) <- NULL
  mydf <- mydf[,c(column_name,cur_names)]
  return(mydf)
}

listDatabases()

# Link modules, EC numbers etc. to KEGG IDs
ko_enzyme.l <-keggLink("ko", "enzyme")
head(ko_enzyme.l)
ko_enzyme.df <- utils::stack(ko_enzyme.l)
ko_enzyme.df$values <- gsub("ko:", "",ko_enzyme.df$values)
ko_enzyme.df$ind <- gsub("ec:", "",ko_enzyme.df$ind)
names(ko_enzyme.df) <- c("KO","EC")
head(ko_enzyme.df)

ko_modules.l <- keggLink("ko", "module")
ko_modules.df <- utils::stack(ko_modules.l)
ko_modules.df$values <- gsub("ko:", "",ko_modules.df$values)
ko_modules.df$ind <- gsub("md:", "",ko_modules.df$ind)
names(ko_modules.df) <- c("KO","Module")
head(ko_modules.df)

ko_pathway.l <- keggLink("ko", "pathway")
ko_pathway.df <- utils::stack(ko_pathway.l)
ko_pathway.df$values <- gsub("ko:", "",ko_pathway.df$values)
ko_pathway.df$ind <- gsub("md:", "",ko_pathway.df$ind)
names(ko_pathway.df) <- c("KO","pathway")
head(ko_pathway.df)


# Get descriptions for modules, kegg IDs, EC numbers etc.
modules_data.df <- keggList("module") %>% utils::stack()
modules_data.df$ind <- gsub("md:", "",modules_data.df$ind)
names(modules_data.df) <- c("Description","Module")
head(modules_data.df)

kegg_data.df <- keggList("ko") %>% utils::stack()
kegg_data.df$ind <- gsub("ko:", "",kegg_data.df$ind)
names(kegg_data.df) <- c("Description","KO")
head(kegg_data.df)

enzyme_data.df <- keggList("enzyme") %>% utils::stack()
enzyme_data.df$ind <- gsub("ec:", "",enzyme_data.df$ind)
names(enzyme_data.df) <- c("Description","EC")
head(enzyme_data.df)

kegg_data.df <- left_join(kegg_data.df, ko_modules.df, by = "KO")
kegg_data.df <- left_join(kegg_data.df, ko_enzyme.df, by = "KO")
kegg_data.df <- left_join(kegg_data.df, enzyme_data.df, by = "EC")
kegg_data.df <- left_join(kegg_data.df, modules_data.df, by = "Module")


#######################################################################
#####################################################################
#Look at the most abundant members of the community using a heatmap
####################################################################

############For now you will learn to group into KOs. 
vsd_df <- m2df(assay(vsd))

anno_DF <-  left_join(vsd_df,
                      kegg_data.df,
                      by = c("Row_variable" = "KO")) ##
anno_DF<- anno_DF %>% distinct(Row_variable, .keep_all= TRUE)
dim(anno_DF)

select <- order(rowMeans((anno_DF[,c(2:9)])),
                decreasing=TRUE)
DF <- anno_DF[select,][1:50,] ##picks the top 50 OTUs
rownames(DF) <- paste(DF$Row_variable, DF$Description.x)

# Add annotation as described above, and change the name of annotation
annotation <- data.frame(Type = factor(1:8, labels = c("Bran", "Bran", "Polystyrene", "Polystyrene", "No food", "No food", "Faeces", "Faeces")))
rownames(annotation) <- colnames(DF[,c(2:9)]) # check out the row names of annotation

###Plot the pretty heatmap!
pheatmap(DF[,c(2:9)], cluster_rows=TRUE, show_rownames=TRUE,
         fontsize_col = 12, annotation = annotation,
         show_colnames = TRUE, cluster_cols=TRUE,
         annotation_legend = TRUE,main = "Top50 abundant KOs", 
         legend = TRUE, 
         border_color = "NA" ,annotation_names_col = TRUE,
         gaps_col = c(2,4,6),
         #gaps_row = c(3,6,13,15,16,18,24,27,30,31,34,40,46,47,48,49),
         filename = "../Figures/Top_50_KO_anno-BlRd.png", cellwidth = 15,cellheight = 15,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdBu")))(500))

dev.off()


############For now you will learn to group into  modules. 

module_anno_df <- anno_DF %>%
                  group_by(anno_DF$Module) %>%
                  dplyr::summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

colnames(module_anno_df)[1] <- "Module"

module_DF <- left_join(module_anno_df,
                      kegg_data.df,
                      by = c("Module" = "Module"))

module_DF<- module_DF %>% distinct(Module, .keep_all= TRUE)
dim(module_DF)

select <- order(rowMeans((module_DF[,c(2:9)])),
                decreasing=TRUE)
DF <- module_DF[select,][1:51,] ##picks the top 50 OTUs
DF <- DF[-1,]

DF <- m2df(DF)
rownames(DF) <- paste(DF$Module, DF$Description)

# Add annotation as described above, and change the name of annotation
annotation <- data.frame(Type = factor(1:8, labels = c("Bran", "Bran", "Polystyrene", "Polystyrene", "No food", "No food", "Faeces", "Faeces")))
rownames(annotation) <- colnames(DF[,c(3:10)]) # check out the row names of annotation



###Plot the pretty heatmap!
pheatmap(DF[,c(3:10)], cluster_rows=TRUE, show_rownames=TRUE,
         fontsize_col = 12, annotation = annotation,
         show_colnames = TRUE, cluster_cols=TRUE,
         annotation_legend = TRUE,main = "Top50 abundant Modules", 
         legend = TRUE, 
         border_color = "NA" ,annotation_names_col = TRUE,
         gaps_col = c(2,4,6),
         #gaps_row = c(3,6,13,15,16,18,24,27,30,31,34,40,46,47,48,49),
         filename = "../Figures/Top_50_module_anno-BlRd.png", cellwidth = 15,cellheight = 15,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdBu")))(500))

dev.off()



#######################################################################################
##############3.5) Determine which KOs are differentially abundant in which sample type
###In this script, we first screen for the differentially abundant KOs in two independent groups at a time.
###Then, we filter the results to get "significant" KOs.###############################
########################################################################################

library(tibble)
###Reorder the deseq2 normalised counts by abundance
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)

###First comparison -  Bran vs Polystyrene
res.BP<- results(dds, c("Type", "Bran", "Polystyrene")) ###extract results from DESeq2 output for the two groups
res.BP <- res.BP[select,][1:25,] ###Look at the top 25 differentially abundant OTUs
head(res.BP)
dim(res.BP)

###Convert the results to a dataframe with tibble
res.BP<- res.BP %>%
  data.frame() %>%
  rownames_to_column()%>% 
  #filter(padj < 1) %>%
  as_tibble()

dim(res.BP)

res.BP <- inner_join(res.BP,
                     kegg_data.df, 
                     by = c("rowname" = "KO"))

res.BP<- res.BP %>% distinct(rowname, .keep_all= TRUE)
res.BP$KO <- factor(paste(res.BP$Description.x, res.BP$rowname))

##Plot the differential/sigificant OTUs for the two groups
BP <- ggplot(data = res.BP,
             aes(x = reorder(KO, log2FoldChange), y = log2FoldChange,
                 fill = log2FoldChange > 0.001))+
  geom_bar(stat = "identity")+
  theme(legend.position = "none",
        axis.text.y = element_text(size = 10))+ 
  geom_hline(yintercept=c(-2,2), linetype="dashed", 
             color = "black", size=1) +
  labs(title = "Bran vs PS", y = "Log2FC", x = "Significant KOs")+
  coord_flip()

###Call the plot function
BP

###Save the plots to file
ggsave("../Figures/BP_Log2FC_KO.png", plot = BP, dpi = 600, units = c("in"), width = 10, height = 5)

###################################################################
##Second comparision: For Bran vs CO
res.BC<- results(dds, c("Type", "Bran", "No food"))
res.BC <- res.BC[select,][1:25,]
head(res.BC)
dim(res.BC)

###Convert the results to a dataframe with tibble
res.BC<- res.BC%>%
  data.frame() %>%
  rownames_to_column()%>% 
  as_tibble() 

dim(res.BC)

res.BC <- inner_join(res.BC,
                     kegg_data.df, 
                     by = c("rowname" = "KO"))

res.BC<- res.BC %>% distinct(rowname, .keep_all= TRUE)
res.BC$KO <- factor(paste(res.BC$Description.x, res.BC$rowname))
##Plot the differential/sigificant OTUs for the two groups
BC <- ggplot(data = res.BC,
             aes(x = reorder(KO, log2FoldChange), y = log2FoldChange,
                 fill = log2FoldChange > 0.001))+
  geom_bar(stat = "identity")+
  theme(legend.position = "none")+ 
  geom_hline(yintercept=c(-2,2), linetype="dashed", 
             color = "black", size=1) +
  labs(title = "Bran vs No food", y = "Log2FC", x = "Significant KOs")+
  coord_flip()

###Call the plot function
BC

###Save the plots to file
ggsave("../Figures/BC_Log2FC.png", plot = BC, dpi = 600, units = c("in"), width = 10, height = 5)

#########################################################################

##Third comparison: For PS vs CO
res.PC<- results(dds, c("Type", "Polystyrene", "No food"))
res.PC <- res.PC[select,][1:25,]
head(res.PC)
dim(res.PC)

###Convert the results to a dataframe with tibble
res.PC<- res.PC%>%
  data.frame() %>%
  rownames_to_column()%>% 
  as_tibble()

dim(res.PC)


res.PC <- inner_join(res.PC,
                     kegg_data.df, 
                     by = c("rowname" = "KO"))

res.PC<- res.PC %>% distinct(rowname, .keep_all= TRUE)
res.PC$KO <- factor(paste(res.PC$Description.x, res.PC$rowname))


##Plot the differential/sigificant OTUs for the two groups
PC <- ggplot(data = res.PC,
             aes(x = reorder(KO,log2FoldChange), y = log2FoldChange,
                 fill = log2FoldChange > 0.001))+
  geom_bar(stat = "identity")+
  theme(legend.position = "none")+ 
  geom_hline(yintercept=c(-2,2), linetype="dashed", 
             color = "black", size=1) +
  labs(title = "Polystyrene vs No food", y = "Log2FC", x = "Significant KOs")+
  coord_flip()

PC

###Save the plots to file
ggsave("../Figures/PC_Log2FC.png", plot = PC, dpi = 600, units = c("in"), width = 10, height = 5)

###############################################################################
##########Investigate individual KOs###########################################
################################################################################

####


ko <- 'K01990'
boxplot(count ~ Type, 
        plotCounts(dds, gene = ko, intgroup = 'Type', 
                  normalized = T, returnData = T), 
                   ylab='Counts', main = ko)

