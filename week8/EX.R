### Load packages
BiocManager::install("zellkonverter")
library("zellkonverter")
library( "scater" )        
library( "scran" ) 
library( "scuttle" )
library("ggplot2")

# Load Gut data from flycellatlas.org
gut <- readH5AD("~/qbb2024-answers/week8/v2_fca_biohub_gut_10x_raw.h5ad")
# Change the assay name from X to counts
assayNames(gut) <- "counts"
# Normalize counts
gut <- logNormCounts(gut)
gut
# Question 1
# How many genes are quantitated? 13407
# How many cells are in the dataset? 11788
# What dimension reduction datasets are present? X_pca X_tsne X_umap

# Question 2
num_columns <- ncol(colData(gut))
# How many columns are there in colData(gut)? 39

colnames(colData(gut))
# Which three column names reported by colnames() seem most interesting? Briefly explain why.
# 1.	annotation: This column likely contains detailed annotations for each cell, such as specific cell types or states.
# 2.	percent_mito: This column represents the percentage of mitochondrial reads for each cell. Mitochondrial content is often used as a quality control metric.
# 3.	n_genes: This column shows the number of unique genes detected in each cell, which helps assess cell quality. 

# Plot cells according to X_umap using plotReducedDim() and colouring by broad_annotation
set.seed(1234)
umap_plot <- plotReducedDim(gut,"X_umap",color="broad_annotation")
# Save the plot
ggsave("umap_plot.png", plot = umap_plot, width = 8, height = 6, dpi = 300)

### 2. Explore data

# Explore gene-level statistics (1 pt)
genecounts = rowSums(assay(gut))

# Question 3: Explore the genecounts distribution (1 pt)

# What is the mean and median genecount according to summary()? What might you conclude from these numbers?
summary(genecounts)
# mean: 3185, median: 254
# Most of the genes have low counts, but a few genes have extremely high counts.

# What are the three genes with the highest expression after using sort()? What do they share in common?
View(sort(genecounts,TRUE))
# IncRNA:Hsromega, pre-rRNA:CR45845, IncRNA:roX1
# They are all non-coding RNA.

# Question 4a: Explore the total expression in each cell across all genes (0.5 pt)

# Create a vector named cellcounts using colSums()
cellcounts <- colSums(assay(gut))
# Create a histogram of cellcounts using hist()
hist(cellcounts)
# What is the mean number of counts per cell?
summary(cellcounts)
# the mean number is 3622
# How would you interpret the cells with much higher total counts (>10,000)?
# may represent cells with a higher transcriptional activity or possibly doublets, where two cells were inadvertently captured together.

# Question 4b: Explore the number of genes detected in each cell (0.5 pt)

# Create a vector named celldetected using colSums() but this time on assay(gut)>0
celldetected <- colSums(assay(gut)>0)
# Create a histogram of celldetected using hist()
hist(celldetected)
# What is the mean number of genes detected per cell?
summary(celldetected)
# mean number: 1059
# What fraction of the total number of genes does this represent?
# 8%

# Explore mitochondrial reads

mito <- grep("^mt:",rownames(gut), value=TRUE)
df = perCellQCMetrics(gut,subsets=list(Mito=mito))
df = as.data.frame(df)
summary(df)
colData(gut) <- cbind( colData(gut), df )

# Question 5: Visualize percent of reads from mitochondria
# Question 5: Visualize percent of reads from mitochondria
# Generate the plot and assign it to a variable
plot <- plotColData(gut, 
                    y = "subsets_Mito_percent", 
                    x = "broad_annotation") +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(x = "Cell Type",
       y = "Percent of Reads from Mitochondria",
       title = "Mitochondrial reads by cell type")
# Save the plot
ggsave(filename = "~/qbb2024-answers/week8/mitochondrial_reads_by_cell_type.png", 
       plot = plot, 
       width = 8, 
       height = 6, 
       dpi = 300)
# Which cell types may have a higher percentage of mitochondrial reads? Why might this be the case?
# Glands, epithelial cells, gut cells, and muscle system have higher mitochondria reads.
# This is because those cells types are metabolically active and have higher energy needs.

### 3. Identify marker genes

# Question 6a: Subset cells annotated as “epithelial cell”
coi <- colData(gut)$broad_annotation == "epithelial cell"
epi <- gut[,coi]
plotUMAP <- plotReducedDim(epi, dimred = "X_umap", colour_by = "annotation") +
  labs(title = "Epithelial cells UMAP",
       x = "UMAP 1",
       y = "UMAP 2")
# Save the UMAP plot
ggsave(filename = "~/qbb2024-answers/week8/umap_epithelial_cells.png", 
       plot = plotUMAP, 
       width = 8, 
       height = 6, 
       dpi = 300)

# find pairwise comparisons between all annotation categories
marker.info <- scoreMarkers(epi, group = colData(epi)$annotation)
# Identify the top marker genes in the anterior midgut 
chosen <- marker.info[["enterocyte of anterior adult midgut epithelium"]]
ordered <- chosen[order(chosen$mean.AUC, decreasing = TRUE),]
head(ordered[, 1:4])

# Question 6b: Evaluate top marker genes

# What are the six top marker genes in the anterior midgut? 
# Mal-A6, Men-b, vnd, betaTry, Mal-A1, Nhe2
# Based on their functions at flybase.org, what macromolecule does this region of the gut appear to specialize in metabolizing?
# saccharide digestion and homostasis. 
# Plot the expression of the top marker gene across cell types
plot_expr <- plotExpression(epi, 
                            features = "Mal-A6", 
                            x = "annotation") +
  labs(title = "Expression of Mal-A6 Across Cell Types",
       x = "Cell Type",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 90))
# Save the plot to the specified directory
ggsave(filename = "~/qbb2024-answers/week8/expression_Mal-A6_across_cell_types.png", 
       plot = plot_expr, 
       width = 8, 
       height = 6, 
       dpi = 300)


# Analyze somatic precursor cells 

# Create a vector indicating cells annotated as "somatic precursor cell"
somatic_coi <- colData(gut)$broad_annotation == "somatic precursor cell"
# subsetting 'gut' with somatic_coi
somatic <- gut[, somatic_coi]
somatic_marker.info <- scoreMarkers(somatic, group = colData(somatic)$annotation)
# Extract and order the top markers for "intestinal stem cell" 
intestinal_stem_chosen <- somatic_marker.info[["intestinal stem cell"]]
ordered_stem <- intestinal_stem_chosen[order(intestinal_stem_chosen$mean.AUC, decreasing = TRUE),]


# Create a vector with the top six marker genes for intestinal stem cells
goi <- rownames(ordered_stem)[1:6]

# Plot the expression of the top six marker genes across cell types
plot_stem_expr <- plotExpression(somatic, 
                                 features = goi, 
                                 x = "annotation") +
  labs(title = "Expression of Top Six Marker Genes for Intestinal Stem Cells",
       x = "Cell Type",
       y = "Expression Level") +
  theme(axis.text.x = element_text(angle = 90))

# Display the plot
print(plot_stem_expr)

# Save the plot to the specified directory
ggsave(filename = "~/qbb2024-answers/week8/expression_top6_intestinal_stem_markers.png", 
       plot = plot_stem_expr, 
       width = 10, 
       height = 8, 
       dpi = 300)
# Which two cell types have more similar expression based on these markers? enterblast and intestinal stem cell
# Which marker looks most specific for intestinal stem cells? DI
