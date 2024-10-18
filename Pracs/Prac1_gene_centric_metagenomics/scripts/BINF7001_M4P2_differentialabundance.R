##############2.6) Determine which taxa are more abundant in which sample type (Differential Abundance)
###In this script, we first screen for the differentially abundant taxa in two independent groups at a time.
###Then, we filter the results to get "significant" OTUS.

#install.packages("tibble")

###Reorder the deseq2 normalised counts by abundance
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)

###First comparison -  Bran vs Polystyrene
res.BP<- results(dds, c("Type", "Bran", "Polystyrene")) ###extract results from DESeq2 output for the two groups
res.BP <- res.BP[select,][1:25,] ###Look at the top 25 differentially abundant taxa
head(res.BP)
dim(res.BP)

###Convert the results to a dataframe with tibble
res.BP<- res.BP %>%
  data.frame() %>%
  rownames_to_column()%>% 
  as_tibble()

dim(res.BP)

##Plot the differential/sigificant OTUs for the two groups
BP <- ggplot(data = res.BP,
             aes(x = reorder(rowname, log2FoldChange), y = log2FoldChange,
                 fill = log2FoldChange > 1))+
  geom_bar(stat = "identity")+
  theme(legend.position = "none")+ 
  geom_hline(yintercept=c(-2,2), linetype="dashed", 
             color = "black", size=1) +
  labs(title = "Bran vs PS", y = "Log2FC", x = "Significant OTUs")+
  coord_flip()

###Call the plot function
BP

###Save the plots to file
ggsave("Figures/BP_Log2FC.png", plot = BP, dpi = 600, units = c("in"), width = 10, height = 5)

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

##Plot the differential/sigificant OTUs for the two groups
BC <- ggplot(data = res.BC,
             aes(x = reorder(rowname, log2FoldChange), y = log2FoldChange,
                 fill = log2FoldChange > 1))+
  geom_bar(stat = "identity")+
  theme(legend.position = "none")+ 
  geom_hline(yintercept=c(-2,2), linetype="dashed", 
             color = "black", size=1) +
  labs(title = "Bran vs Starvation", y = "Log2FC", x = "Significant OTUs")+
  coord_flip()

###Call the plot function
BC

###Save the plots to file
ggsave("Figures/BC_Log2FC.png", plot = BC, dpi = 600, units = c("in"), width = 10, height = 5)

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

##Plot the differential/sigificant OTUs for the two groups
PC <- ggplot(data = res.PC,
             aes(x = reorder(rowname, log2FoldChange), y = log2FoldChange,
                 fill = log2FoldChange > 1))+
  geom_bar(stat = "identity")+
  theme(legend.position = "none")+ 
  geom_hline(yintercept=c(-2,2), linetype="dashed", 
             color = "black", size=1) +
  labs(title = "Polystyrene vs Starvation", y = "Log2FC", x = "Significant OTUs")+
  coord_flip()

PC

###Save the plots to file
ggsave("Figures/PC_Log2FC.png", plot = PC, dpi = 600, units = c("in"), width = 10, height = 5)

###############################################################################
###Try again to identify significant OTUs! Copy over the entire script below this.
###Hint: add pvalue<1 when converting the results to a dataframe with tibble. This means the results are 99% significant
################################################################################

##############2.6) Determine which taxa are more abundant in which sample type (Differential Abundance)
###In this script, we first screen for the differentially abundant taxa in two independent groups at a time.
###Then, we filter the results to get "significant" OTUS.

#install.packages("tibble")

###Reorder the deseq2 normalised counts by abundance
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)

###First comparison -  Bran vs Polystyrene
res.BP<- results(dds, c("Type", "Bran", "Polystyrene")) ###extract results from DESeq2 output for the two groups
res.BP <- res.BP[select,][1:25,] ###Look at the top 25 differentially abundant taxa
head(res.BP)
dim(res.BP)

###Convert the results to a dataframe with tibble
res.BP<- res.BP %>%
  data.frame() %>%
  rownames_to_column()%>% 
  filter(padj < 1)  %>%
  as_tibble()

dim(res.BP)

##Plot the differential/sigificant OTUs for the two groups
BP <- ggplot(data = res.BP,
             aes(x = reorder(rowname, log2FoldChange), y = log2FoldChange,
                 fill = log2FoldChange > 1))+
  geom_bar(stat = "identity")+
  theme(legend.position = "none")+ 
  geom_hline(yintercept=c(-2,2), linetype="dashed", 
             color = "black", size=1) +
  labs(title = "Bran vs PS", y = "Log2FC", x = "Significant OTUs")+
  coord_flip()

###Call the plot function
BP

###Save the plots to file
ggsave("Figures/BP_Log2FC_sig.png", plot = BP, dpi = 600, units = c("in"), width = 10, height = 5)

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
  filter(padj < 1)  %>%
  as_tibble() 

dim(res.BC)

##Plot the differential/sigificant OTUs for the two groups
BC <- ggplot(data = res.BC,
             aes(x = reorder(rowname, log2FoldChange), y = log2FoldChange,
                 fill = log2FoldChange > 1))+
  geom_bar(stat = "identity")+
  theme(legend.position = "none")+ 
  geom_hline(yintercept=c(-2,2), linetype="dashed", 
             color = "black", size=1) +
  labs(title = "Bran vs Starvation", y = "Log2FC", x = "Significant OTUs")+
  coord_flip()

###Call the plot function
BC

###Save the plots to file
ggsave("Figures/BC_Log2FC_sig.png", plot = BC, dpi = 600, units = c("in"), width = 10, height = 5)

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
  filter(padj < 1)  %>%
  as_tibble()

dim(res.PC)

##Plot the differential/sigificant OTUs for the two groups
PC <- ggplot(data = res.PC,
             aes(x = reorder(rowname, log2FoldChange), y = log2FoldChange,
                 fill = log2FoldChange > 1))+
  geom_bar(stat = "identity")+
  theme(legend.position = "none")+ 
  geom_hline(yintercept=c(-2,2), linetype="dashed", 
             color = "black", size=1) +
  labs(title = "Polystyrene vs Starvation", y = "Log2FC", x = "Significant OTUs")+
  coord_flip()

PC

###Save the plots to file
ggsave("Figures/PC_Log2FC_sig.png", plot = PC, dpi = 600, units = c("in"), width = 10, height = 5)