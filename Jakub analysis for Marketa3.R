#These are the libraries usued to import the idat files and process them

library(ggplot2)
library(ggrepel)
library(ChAMP)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(ggfortify)
library(minfi)
library(PCAtools)
library(limma)

#I have imported all samples that were run alongside the ovarian cell lines
#this is so that I can effectively determine the effect on the data on the samples being run on different microarray chips
#once I have corrected for chip-to-chip variation, I remove all the samples I am not interested in
setwd("D:/01 Datasets/Moran cell line 450K methylation/champ/idats") #the location of the idat files
myLoad <- champ.load(directory=getwd(),arraytype="450K", filterXY=FALSE) #importing the idat files
champ.QC() #visualising some quality scores
myNorm <- champ.norm(plotBMIQ=TRUE,arraytype="450K") #normalising data using BMIQ method which outperforms other methods
champ.SVD() #SVD is singular value decomposition to see if any of the PCs are related to the microarray chip the sample was run on i.e. a batch effect
myCombat<- champ.runCombat(beta=myNorm,pd=myLoad$pd,batchname=c("Slide")) #as SVD detected PCs associated with chip to chip variation, I correct this using combat -- it adjusts the values to remove variation while preserving biological diffeences
champ.SVD(beta=myCombat) #repeating SVD on the batch corrected data. There is now much less PCs associated with chip-to-chip variation
targets <- myLoad[["pd"]] #extracting targets/annotation file from initial import
colnames(myCombat) <- targets$GSM_number #Adding sample names to the m value matrix
targets <- targets[targets$Sample_Group == "ovary", ] #reducing the annoation file down from all samples to just ovarian ones
sample.names <- targets$GSM_number #making a list of ovarian cell lines
beta_vals <- myCombat[,sample.names] #selecting only ovarian cell lines from m value matrix
m_vals <-logit2(beta_vals)
rm(Anno,hm450.manifest.hg19,multi.hit,myCombat,myLoad,myNorm,probeInfoALL.lv) #remove things I no longer need from the work environment
#save.image(file="imported_idats_champ.RData")

pca_prcomp <- prcomp(t(m_vals), center = TRUE, scale = FALSE)
autoplot(pca_prcomp, label = TRUE, label.size = 3)

#all ovarian cell lines ~45 cell lines
write.csv(m_vals, file="ovarian_mvals.csv")
write.csv(beta_vals,file="ovarian_betavals.csv")
write.csv(targets, file="ovarian_targets.csv")

#select lines I experimented on 
targets.decitabine <-read.csv(file='ovarian_targets.csv')
targets.decitabine <-targets.decitabine[targets.decitabine$Combination!= "Unknown", ]
sample.names <- targets.decitabine$GSM_number

m_vals.decitabine <- read.csv(file='ovarian_mvals.csv')
m_vals.decitabine <- m_vals.decitabine[,sample.names]

beta_vals.decitabine <- read.csv(file='ovarian_betavals.csv')
beta_vals.decitabine <- beta_vals.decitabine[,sample.names]

write.csv(m_vals.decitabine, file="pargi_mvals.csv")
write.csv(beta_vals.decitabine,file="pargi_betavals.csv")
write.csv(targets.decitabine, file="pargi_targets.csv")

rownames(targets.decitabine) <- targets.decitabine$GSM_number #make colnames and rownames between mvals and targets the same

#PCA and look at pairs plot for decitabine sensitivity
p <- pca(m_vals.decitabine, metadata = targets.decitabine, removeVar = 0.1) #pca using the pcatools package
pairsplot(p,  colby = 'Combination')
#look at pc1 versus 2 coloured by decitabine sensitivity
biplot(p, x = 'PC1', y = 'PC2', lab = NULL, colby = 'Combination', pointSize = 5, legendPosition = 'bottom')

#limma code to determine which CpG sites are differentially methylated
Combination <- factor(targets.decitabine$Combination) #need to create a factor stating which cell lines are decitabine sensitive/resistant
design <- model.matrix(~0+Combination) #create a design matrix, which is like a binary matrix with 1 and o
colnames(design) <- c("Resistant", "Sensitive") #tidy up design matrix, 
fit <- lmFit(m_vals.decitabine, design) #fit the model
cont.matrix <- makeContrasts(Combination=Sensitive-Resistant, levels=design) #define what contrast I want to make
fit2 <- contrasts.fit(fit, cont.matrix) #make contrast I'm interested in within model
fit2 <- eBayes(fit2) #use eBayes() to moderate the estimated error variances
topTable(fit2, adjust="BH") #look at top hits, gives  an idea of p values but not very informative as all probes just have cg codes
Combination_sens <- topTable(fit2, coef="Combination", number=Inf, p.value=1, sort.by="p") #this extracts the stats for all probes even the non-significantly differentially methylated ones

anno <- read.csv(file="HumanMethylation450_15017482_v1-2.csv", row.names=1) #read in an annotation file to convert the cpg ids to gene names
anno_hits <- merge(Combination_sens, anno, by=0) #then I merge the statistics I just did with the annotation
#merging did not work as there are some differences in naming structure so I adjust the names

rownames(Combination_sens) <- paste0("cg", rownames(Combination_sens))
rownames(Combination_sens) <- substr(rownames(Combination_sens), 1, 8)

unique_ids <- make.unique(substr(rownames(Combination_sens), 1, 6))
rownames(Combination_sens) <- unique_ids

unique_ids <- make.unique(substr(rownames(anno), 1, 6))
rownames(anno) <- unique_ids

Combination_sens$id <- substr(rownames(Combination_sens), 1, 6)
anno$id <- substr(rownames(anno), 1, 6)

anno_hits <- merge(Combination_sens, anno, by="id")


write.csv(anno_hits, file="Combination_senstive_dmps.csv")

#This way I got a file with genes which had altered methylation status in cell lines that DNMTi treatment sensitised to the standard treatment, versus cell lines that DNMTi did not change standrad treatment response. 
#Now I proceed to some visualisations of result.

#Volcano plot
significant_points <- anno_hits[anno_hits$P.Value < 0.05 & abs(anno_hits$logFC) > 5, ]
top_points <- head(anno_hits[order(anno_hits$P.Value), ], 10)
anno_hits$negLogPval <- -log10(anno_hits$P.Value)

volcano_plot <- ggplot(anno_hits, aes(x = logFC, y = -log10(P.Value))) +
  geom_point() +
  geom_point(data = top_points, aes(x = logFC, y = -log10(P.Value)), color = "red") +
  geom_text_repel(data = top_points,
                  aes(x = logFC, y = -log10(P.Value), label = UCSC_RefGene_Name),
                  size = 3,
                  max.overlaps = 10) 

print(volcano_plot)

#heatmap

data$logFC <- as.numeric(data$logFC)  
top40 <- data[order(data$adj.P.Val, na.last = TRUE), ][1:40, ] 

heatmap_data <- matrix(top40$logFC, nrow = 40, ncol = 1) 
rownames(heatmap_data) <- top40$Row.names

gene_annotations <- data.frame(Gene = top40$UCSC_RefGene_Name)
rownames(gene_annotations) <- top40$Row.names

pheatmap(heatmap_data,
         annotation_row = gene_annotations,
         show_rownames = TRUE,
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         cluster_cols = FALSE,  
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Top 40 Significant Methylation Changes")

#Gene function enrichment analysis

library(clusterProfiler)
library(org.Hs.eg.db)

genes <- unique(data$UCSC_RefGene_Name[!is.na(data$UCSC_RefGene_Name)])

ego <- enrichGO(gene         = genes,
                OrgDb        = org.Hs.eg.db,
                keyType      = 'SYMBOL',
                ont          = "MF",  #"MF" for Molecular Functions, "CC" for Cellular Components
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable     = TRUE)

dotplot(ego)

