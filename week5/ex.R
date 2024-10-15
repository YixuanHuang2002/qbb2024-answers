# Load libraries we'll need
library(DESeq2)
library(vsn)
library(matrixStats)
library(readr)
library(dplyr)
library(tibble)
library(ggfortify)

# Step 3.1: Loading and filtering the data

# Load tab-separated data file
data = readr::read_tsv("~/qbb2024-answers/week5/salmon.merged.gene_counts.tsv")
# Change gene names into row names
data = column_to_rownames(data, var="gene_name")
# Get rid of gene id column
data = data %>% select(-gene_id)
# Change data to integers
data = data %>% mutate_if(is.numeric, as.integer)
# Keep rows with at least 100 reads
data = data[rowSums(data) > 100,]
# select only the “narrow region” samples
narrow_region= dplyr::select(data,"A1_Rep1":"P2-4_Rep3")

# Step 3.2: Creating DESeq2 model and batch-correction
metadata = tibble(tissue=as.factor(c("A1", "A1", "A1",
                                     "A2-3","A2-3","A2-3",
                                     "Cu","Cu","Cu",
                                     "LFC-Fe","LFC-Fe","LFC-Fe",
                                     "Fe","Fe","Fe",
                                     "P1","P1","P1",
                                     "P2-4","P2-4","P2-4")),
                  rep=as.factor(c("Rep1","Rep2","Rep3",
                                  "Rep1","Rep2","Rep3",
                                  "Rep1","Rep2","Rep3",
                                  "Rep1","Rep2","Rep3",
                                  "Rep1","Rep2","Rep3",
                                  "Rep1","Rep2","Rep3",
                                  "Rep1","Rep2","Rep3")))
# Create a DESeq data
DESqdata = DESeqDataSetFromMatrix(countData=as.matrix(narrow_region),
                                    colData=metadata,
                                    design=~tissue)
# correct the data for batch-effects
Vstdata = vst(DESqdata)
meanSdPlot(assay(Vstdata))

# Step 3.3: PCA analysis

# Perform PCA on the corrected data and plot it
pca_data = plotPCA(Vstdata,
                   intgroup=c("rep","tissue"),
                   returnData=TRUE)
ggplot(pca_data, aes(PC1, PC2, color=tissue, shape=rep)) +
  geom_point(size=5)+
  labs(title = "PCA Plot (wrong)")
ggsave("~/qbb2024-answers/week5/PCA_Plot_Wrong.png")


## corrected the data by make changes in metadata

# Step 3.2: Creating DESeq2 model and batch-correction
metadata = tibble(tissue=as.factor(c("A1", "A1", "A1",
                                     "A2-3","A2-3","A2-3",
                                     "Cu","Cu","Cu",
                                     "LFC-Fe","LFC-Fe","Fe", ## change1
                                     "LFC-Fe","Fe","Fe", ## change2
                                     "P1","P1","P1",
                                     "P2-4","P2-4","P2-4")),
                  rep=as.factor(c("Rep1","Rep2","Rep3",
                                  "Rep1","Rep2","Rep3",
                                  "Rep1","Rep2","Rep3", 
                                  "Rep1","Rep2","Rep3",
                                  "Rep1","Rep2","Rep3",
                                  "Rep1","Rep2","Rep3",
                                  "Rep1","Rep2","Rep3")))
# Create a DESeq data
DESqdata = DESeqDataSetFromMatrix(countData=as.matrix(narrow_region),
                                  colData=metadata,
                                  design=~tissue)
# correct the data for batch-effects
Vstdata = vst(DESqdata)
meanSdPlot(assay(Vstdata))

# Step 3.3: PCA analysis

# Perform PCA on the corrected data and plot it
pca_data = plotPCA(Vstdata,
                   intgroup=c("rep","tissue"),
                   returnData=TRUE)
ggplot(pca_data, aes(PC1, PC2, color=tissue, shape=rep)) +
  geom_point(size=5)+
  labs(title = "PCA Plot (wrong)")
ggsave("~/qbb2024-answers/week5/PCA_Plot_Correct.png")

# Step 3.4: Filtering genes by variance

mat = as.matrix(assay(Vstdata))
# standard deviations
combined = mat[,seq(1, 21, 3)]
combined = combined + mat[,seq(2, 21, 3)]
combined = combined + mat[,seq(3, 21, 3)]
combined = combined / 3
sds = rowSds(combined)
# keep only genes with a standard deviation greater than 1
mat = mat[rowSds(combined) > 1,]

# Step 3.5: K-means clustering genes

# first set the seed to 42
set.seed(42)
# then perform the clustering
k=kmeans(mat, centers=12)$cluster
# Find ordering of samples to put them into their clusters
ordering = order(k)
# Reorder genes
k = k[ordering]
# Plot heatmap of expressions and clusters
png("~/qbb2024-answers/week5/3.5 heatmap.png", width = 800, height = 800)
heatmap(mat[ordering,],
        Rowv=NA,
        Colv=NA,
        RowSideColors = RColorBrewer::brewer.pal(12,"Paired")[k])
# Close the PNG device to save the file
dev.off()

# Step 3.6: Gene ontology enrichment analysis

# Get the gene names using rownames for all genes in cluster 1.
genes = rownames(mat[k == 1,])
# save the txt file
write.table(genes, file = "~/qbb2024-answers/week5/cluster1_genes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)







