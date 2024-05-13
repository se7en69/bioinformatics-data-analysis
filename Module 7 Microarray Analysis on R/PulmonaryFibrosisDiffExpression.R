packages <- c("GEOquery", "dplyr", "limma",'ggplot2','pheatmap','ggrepel','data.table')

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      BiocManager::install(x)
      library(x, character.only = TRUE)
    }
  }
)

#data from this study - https://pubmed.ncbi.nlm.nih.gov/23783374/
#reference - https://sbc.shef.ac.uk/geo_tutorial/tutorial.nb.html

#retrieve the microarray data and store the first dataset
my_id <- "GSE32537"
gse <- getGEO(my_id)
gse <- gse[[1]]


sampleInfo <- pData(gse) ## get sample information
anno <- fData(gse) ## get gene annotation data

#check if the data is normalized, distributions should be similar
boxplot(exprs(gse),outline=FALSE)

#There looks to be some clustering of samples with control diagnosis
heatMapSampleInfo <-  dplyr::select(sampleInfo, sample=geo_accession, diagnosis=`final diagnosis:ch1`)
corMatrix <- cor(exprs(gse),use="c")
rownames(heatMapSampleInfo) <- colnames(corMatrix)
pheatmap(corMatrix, annotation_col=heatMapSampleInfo)  

#I hypothesize age and pack years of smoking have the biggest impact on diagnosis, so I'd like to include them in PCA
sampleInfo <- dplyr::select(sampleInfo, sample=geo_accession, age=`age:ch1`, pack_years=`pack years:ch1`, diagnosis=`final diagnosis:ch1`)

#create principle component analysis to see how much each principle component explains diagnosis
#control diagnoses cluster to the right and IPF/UIP cluster to the left on PC1
pca <- prcomp(t(exprs(gse)))
cbind(sampleInfo, pca$x) %>%
  ggplot(aes(x = PC1, y= PC2, col=diagnosis, label = paste("Sample", sample))) + geom_point() + geom_text_repel(size = 2) 



#create model and rename columns for diagnoses
design <- model.matrix(~0+sampleInfo$diagnosis)
colnames(design) <- c("control", "COP", "DIP", "IPF.UIP","NSIP", "RB.ILD", "UF")

#calculate median expression level
cutoff <- median(exprs(gse))
is_expressed <- exprs(gse) > cutoff

#identify genes expressed in 2 or more samples
keep <- rowSums(is_expressed) > 2

#subset those genes
gse <- gse[keep,]

geneAssignment <- data.table(fData(gse)$gene_assignment)
gse@featureData@data$Symbol<-geneAssignment[,tstrsplit(V1, " // ")[2]]
anno <- fData(gse) #update anno with gene symbols

#create coefficient matrix for contrasts
x <- c("control-COP", "COP-DIP", "DIP-IPF.UIP", "IPF.UIP-NSIP",
       "NSIP-RB.ILD", "RB.ILD-UF", "UF-control")
contrasts <- makeContrasts(contrasts = x, levels=design)

## calculate relative array weights, fit linear model, and rank genes in order of evidence using emperical Bayes
aw <- arrayWeights(exprs(gse),design)
fit <- lmFit(exprs(gse), design,
             weights = aw)
contrasts <- makeContrasts(x, levels=design)
weightedfit2 <- contrasts.fit(fit, contrasts)
weightedfit2 <- eBayes(weightedfit2)

#add gene symbols for labeling volcano plot
weightedfit2$genes <- anno$Symbol

results <- topTable(weightedfit2, number=Inf)
results <- tibble::rownames_to_column(results,"ID")

#create volcano plot, labeling top 20 genes, turning significant (based on p-value and fold change) results blue
p_cutoff <- 0.005
fc_cutoff <- 1
topN <- 21

volcanoPlot <- results %>% 
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  mutate(Rank = 1:n(), Label = ifelse(Rank < topN, V1,"")) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant,label=Label)) + geom_point() + geom_text_repel(col="black") + ylab('log-odds (B)')

print(volcanoPlot)

# Print top 20 most significant differentially expressed genes
head(volcanoPlot$data, 20)

