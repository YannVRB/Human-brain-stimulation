#Load librairies
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVAR)
library(motifmatchr)
library(ggplot2)

set.seed(1234)

#Load the RNA and ATAC data
counts4A <- Read10X_h5("sample4A_filtered_feature_bc_matrix.h5")
counts4B <- Read10X_h5("sample4B_filtered_feature_bc_matrix.h5")
counts5A <- Read10X_h5("sample5A_filtered_feature_bc_matrix.h5")
counts5B <- Read10X_h5("sample5B_filtered_feature_bc_matrix.h5")
fragpath4A <- "/Shared/NEURO/AbelLab/Yann/Multiomic_humanBrainStimulation/sample4A/outs/atac_fragments.tsv.gz"
fragpath4B <- "/Shared/NEURO/AbelLab/Yann/Multiomic_humanBrainStimulation/sample4B/outs/atac_fragments.tsv.gz"
fragpath5A <- "/Shared/NEURO/AbelLab/Yann/Multiomic_humanBrainStimulation/sample5A/outs/atac_fragments.tsv.gz"
fragpath5B <- "/Shared/NEURO/AbelLab/Yann/Multiomic_humanBrainStimulation/sample5B/outs/atac_fragments.tsv.gz"

#Get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

#Create a Seurat object containing the RNA adata
sample4A <- CreateSeuratObject(
  counts = counts4A$`Gene Expression`,
  assay = "RNA"
)
sample4B <- CreateSeuratObject(
  counts = counts4B$`Gene Expression`,
  assay = "RNA"
)
sample5A <- CreateSeuratObject(
  counts = counts5A$`Gene Expression`,
  assay = "RNA"
)
sample5B <- CreateSeuratObject(
  counts = counts5B$`Gene Expression`,
  assay = "RNA"
)

#create ATAC assay and add it to the object
sample4A[["ATAC"]] <- CreateChromatinAssay(
  counts = counts4A$Peaks,
  sep = c(":", "-"),
  fragments = fragpath4A,
  annotation = annotation
)
sample4B[["ATAC"]] <- CreateChromatinAssay(
  counts = counts4B$Peaks,
  sep = c(":", "-"),
  fragments = fragpath4B,
  annotation = annotation
)
sample5A[["ATAC"]] <- CreateChromatinAssay(
  counts = counts5A$Peaks,
  sep = c(":", "-"),
  fragments = fragpath5A,
  annotation = annotation
)
sample5B[["ATAC"]] <- CreateChromatinAssay(
  counts = counts5B$Peaks,
  sep = c(":", "-"),
  fragments = fragpath5B,
  annotation = annotation
)

DefaultAssay(sample4A) <- "ATAC"
DefaultAssay(sample4B) <- "ATAC"
DefaultAssay(sample5A) <- "ATAC"
DefaultAssay(sample5B) <- "ATAC"

#Quality control
sample4A <- NucleosomeSignal(sample4A)
sample4B <- NucleosomeSignal(sample4B)
sample5A <- NucleosomeSignal(sample5A)
sample5B <- NucleosomeSignal(sample5B)
sample4A <- TSSEnrichment(sample4A)
sample4B <- TSSEnrichment(sample4B)
sample5A <- TSSEnrichment(sample5A)
sample5B <- TSSEnrichment(sample5B)

DefaultAssay(sample4A) <- "RNA"
DefaultAssay(sample4B) <- "RNA"
DefaultAssay(sample5A) <- "RNA"
DefaultAssay(sample5B) <- "RNA"

sample4A[["percent.mt"]] <- PercentageFeatureSet(sample4A, pattern = "^MT-")
sample4B[["percent.mt"]] <- PercentageFeatureSet(sample4B, pattern = "^MT-")
sample5A[["percent.mt"]] <- PercentageFeatureSet(sample5A, pattern = "^MT-")
sample5B[["percent.mt"]] <- PercentageFeatureSet(sample5B, pattern = "^MT-")

#Filter out low quality cells
sample4A <- subset(
  x = sample4A,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1 &
    percent.mt < 20
)
sample4B <- subset(
  x = sample4B,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1 &
    percent.mt < 20
)
sample5A <- subset(
  x = sample5A,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1 &
    percent.mt < 20
)
sample5B <- subset(
  x = sample5B,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1 &
    percent.mt < 20
)

DefaultAssay(sample4A) <- "ATAC"
DefaultAssay(sample4B) <- "ATAC"
DefaultAssay(sample5A) <- "ATAC"
DefaultAssay(sample5B) <- "ATAC"

#Call peaks using MACS2
peaks4A <- CallPeaks(sample4A, macs2.path = "/Shared/NEURO/AbelLab/Yann/bcbio/anaconda/bin/macs2")
peaks4B <- CallPeaks(sample4B, macs2.path = "/Shared/NEURO/AbelLab/Yann/bcbio/anaconda/bin/macs2")
peaks5A <- CallPeaks(sample5A, macs2.path = "/Shared/NEURO/AbelLab/Yann/bcbio/anaconda/bin/macs2")
peaks5B <- CallPeaks(sample5B, macs2.path = "/Shared/NEURO/AbelLab/Yann/bcbio/anaconda/bin/macs2")

#Remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks4A <- keepStandardChromosomes(peaks4A, pruning.mode = "coarse")
peaks4B <- keepStandardChromosomes(peaks4B, pruning.mode = "coarse")
peaks5A <- keepStandardChromosomes(peaks5A, pruning.mode = "coarse")
peaks5B <- keepStandardChromosomes(peaks5B, pruning.mode = "coarse")
peaks4A <- subsetByOverlaps(x = peaks4A, ranges = blacklist_hg38_unified, invert = TRUE)
peaks4B <- subsetByOverlaps(x = peaks4B, ranges = blacklist_hg38_unified, invert = TRUE)
peaks5A <- subsetByOverlaps(x = peaks5A, ranges = blacklist_hg38_unified, invert = TRUE)
peaks5B <- subsetByOverlaps(x = peaks5B, ranges = blacklist_hg38_unified, invert = TRUE)

#Quantify counts in each peak
macs2_counts4A <- FeatureMatrix(
  fragments = Fragments(sample4A),
  features = peaks4A,
  cells = colnames(sample4A)
)
macs2_counts4B <- FeatureMatrix(
  fragments = Fragments(sample4B),
  features = peaks4B,
  cells = colnames(sample4B)
)
macs2_counts5A <- FeatureMatrix(
  fragments = Fragments(sample5A),
  features = peaks5A,
  cells = colnames(sample5A)
)
macs2_counts5B <- FeatureMatrix(
  fragments = Fragments(sample5B),
  features = peaks5B,
  cells = colnames(sample5B)
)

#Create a new assay using the MACS2 peak set and add it to the Seurat object
sample4A[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts4A,
  fragments = fragpath4A,
  annotation = annotation
)
sample4B[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts4B,
  fragments = fragpath4B,
  annotation = annotation
)
sample5A[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts5A,
  fragments = fragpath5A,
  annotation = annotation
)
sample5B[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts5B,
  fragments = fragpath5B,
  annotation = annotation
)

#Switching to RNA assay
DefaultAssay(sample4A) <- "RNA"
DefaultAssay(sample4B) <- "RNA"
DefaultAssay(sample5A) <- "RNA"
DefaultAssay(sample5B) <- "RNA"

#Normalize gene expression
sample4A <- NormalizeData(sample4A)
sample4B <- NormalizeData(sample4B)
sample5A <- NormalizeData(sample5A)
sample5B <- NormalizeData(sample5B)

#Find most variable features in each cell
sample4A <- FindVariableFeatures(sample4A, selection.method = "vst", nfeatures = 2000)
sample4B <- FindVariableFeatures(sample4B, selection.method = "vst", nfeatures = 2000)
sample5A <- FindVariableFeatures(sample5A, selection.method = "vst", nfeatures = 2000)
sample5B <- FindVariableFeatures(sample5B, selection.method = "vst", nfeatures = 2000)

#Assign condition before integration
sample4A$dataset <- "Baseline"
sample4B$dataset <- "Stimulated"
sample5A$dataset <- "Baseline"
sample5B$dataset <- "Stimulated"

#Integration
#Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list(sample4A,sample4B,sample5A,sample5B))

#Perform integration
anchors <- FindIntegrationAnchors(object.list = list(sample4A,sample4B,sample5A,sample5B), anchor.features = features)

#This command creates an 'integrated' data assay
combined <- IntegrateData(anchorset = anchors)

#Specify that we will perform downstream analysis on the corrected data note that the
#Original unmodified data still resides in the 'RNA' assay
DefaultAssay(combined) <- "integrated"

#Run the standard workflow for visualization and clustering
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)

#Annotating cell types using a reference scRNA-seq from Allen
#For each matrix: read in, make Seurat and remove from memory
sc.list = lapply(subclass, function(i) {
  #Inner join to ensure all cells have metadata
  dt = merge(counts, meta[subclass_label == i], by="sample_name")
  #Transpose counts (g x sn), get rownames and meta.data from dt3
  x = CreateSeuratObject(counts=transpose(dt[, 1:gn], make.names="sample_name"), 
                         row.names=names(dt)[2:gn], 
                         project="AIBS_mouse_cort-hipp_10x",
                         meta.data=data.frame(dt[, names(meta), with=F], row.names=1))
  #Do you wish to randomly downsample? (https://github.com/satijalab/seurat/issues/3108)
  x = x[, sample(colnames(x), size=min(20000, nrow(dt)), replace=F)]
  remove(dt)
  #Update and exit
  print(paste0("I finished ", i))
  return(x)
})

#Merge and save
print(paste("sc.list size:", format(object.size(sc.list), units="auto")))
#Quick&dirty: takes huge amounts of memory (something about merge)
allen = Reduce(merge, sc.list) 
print(paste("allen size:", format(object.size(allen), units="auto")))
saveRDS(allen, file="allen_cortex.rds")
remove(sc.list)

#Normalizaing and scaling the scRNA-seq reference
allen <- NormalizeData(allen)
allen <- FindVariableFeatures(allen)
allen <- RunPCA(allen)
allen <- RunUMAP(allen, dims = 1:30)

#Transfer labels with reference scRNA-seq data
transfer.anchors <- FindTransferAnchors(
  reference = allen,
  query = combined,
  reference.reduction = 'pca',
  dims = 1:30
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = allen$subclass_label,
  dims = 1:30
)

combined <- AddMetaData(combined, metadata = predicted.labels)

#Visualization
DimPlot(combined, group.by = 'predicted.id', pt.size = 0.1, label = TRUE) + NoLegend()

#Set the cell identities to the cell type predictions
for(i in levels(combined)) {
  cells_to_reid <- WhichCells(combined, idents = i)
  newid <- names(sort(table(combined$predicted.id[cells_to_reid]),decreasing=TRUE))[1]
  Idents(combined, cells = cells_to_reid) <- newid
}

#Rename cluster identification to group major celltypes together
combined_maincelltypeslabeled <- RenameIdents(combined, `L2/3 IT` = "excitatory neurons", `L5 IT` = "excitatory neurons", `L6 IT Car3` = "excitatory neurons", `Pvalb` = "Pvalb-Sst inhibitory neurons", `L5/6 NP` = "excitatory neurons", `Lamp5` = "Vip-Sncg-Lamp5 inhibitory neurons", `Sncg` = "Vip-Sncg-Lamp5 inhibitory neurons", `Sst` = "Pvalb-Sst inhibitory neurons", `Vip` = "Vip-Sncg-Lamp5 inhibitory neurons")

#Switching to RNA assay for differential gene expression analysis
DefaultAssay(combined_maincelltypeslabeled) <- 'RNA'
combined_maincelltypeslabeled$celltype.dataset <- paste(Idents(combined_maincelltypeslabeled), combined_maincelltypeslabeled$dataset, sep = "_")
combined_maincelltypeslabeled$celltype <- Idents(combined_maincelltypeslabeled)
Idents(combined_maincelltypeslabeled) <- "celltype.dataset"

#Example of differential gene expression analysis in excitory neurons in stimulated vs baseline conditions
exc.neur.DEG <- FindMarkers(combined_maincelltypeslabeled, ident.1 = "excitatory neurons_Stimulated", ident.2 = "excitatory neurons_Baseline", min.pct=0, logfc.threshold=0)
write.csv(exc.neur.DEG,'exc.neur.DEG.csv')

#Example of violin plot for cytokine activity genes
VlnPlot(combined_maincelltypeslabeled, features = c("IL1B", "CCL2", "IL1A", "OSM", "RGS1", "CCL3", "CCL4"), pt.size = 1, y.max=5, ncol = 4, group.by = "celltype", split.by="dataset")

#Example of violin plot for DNA-binding transcription activation activity
VlnPlot(combined_maincelltypeslabeled, features = c("BTG2", "FOSB", "FOS", "NR4A1", "EGR1", "ATF3", "JUNB", "EGR2"), pt.size = 1, y.max=4, ncol = 4, group.by = "celltype", split.by="dataset")

#DNA accessibility data processing
DefaultAssay(combined_maincelltypeslabeled) <- "peaks"
combined_maincelltypeslabeled <- FindTopFeatures(combined_maincelltypeslabeled, min.cutoff = 5)
combined_maincelltypeslabeled <- RunTFIDF(combined_maincelltypeslabeled)
combined_maincelltypeslabeled <- RunSVD(combined_maincelltypeslabeled)

#Linking peaks to genes
#First compute the GC content for each peak
combined_maincelltypeslabeled_splitdatasetforcoverageplot <- RegionStats(combined_maincelltypeslabeled_splitdatasetforcoverageplot, genome = BSgenome.Hsapiens.UCSC.hg38)
#Link peaks to genes
combined_maincelltypeslabeled <- LinkPeaks(
  object = combined_maincelltypeslabeled,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = c("CCL4", "NPAS4")
)

da_peaks_excNeur <- FindMarkers(
  object = combined_maincelltypeslabeled,
  ident.1 = "excitatory neurons_Stimulated",
  ident.2 = "excitatory neurons_Baseline",
  min.pct = 0.05,
  logfc.threshold=0
)

#Setup windows of 1000bp before and after the TSS of each gene
window1.1000TSS1000 <- GeneActivity(combined_maincelltypeslabeled, extend.upstream = 1000, extend.downstream = 1000, biotypes = "protein_coding")
#Attach to the combined seurat object
combined_maincelltypeslabeled[['X1000TSS1000']] <- CreateAssayObject(counts = window1.1000TSS1000)
#Nornalize the new assay
combined_maincelltypeslabeled <- NormalizeData(
  object = combined_maincelltypeslabeled,
  assay = 'X1000TSS1000',
  normalization.method = 'LogNormalize',
  scale.factor = median(combined_maincelltypeslabeled$nCount_peaks)
)

#Switch to the new assay
DefaultAssay(combined_maincelltypeslabeled) <- "X1000TSS1000"

#Differentially accessible analysis around the TSS
da_peaks_TSSmicrogliaLRtest <- FindMarkers(
  object = combined_maincelltypeslabeled,
  ident.1 = "Micro-PVM_Stimulated",
  ident.2 = "Micro-PVM_Baseline",
  min.pct = 0.05,
  logfc.threshold=0, latent.vars = 'nCount_peaks', test.use = 'LR'
)

#Extract position frequency matrices for the motifs
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

#Add motif information
combined_maincelltypeslabeled <- AddMotifs(
  object = combined_maincelltypeslabeled,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

#Get top differentially accessible peaks in microglia cluster
top.da.peak.microglia <- rownames(da_peaks_TSSmicrogliaLRtest[da_peaks_TSSmicrogliaLRtest$p_val < 0.005, ])

#Test enrichment
enriched.motifs <- FindMotifs(
  object = mouse_brain,
  features = top.da.peak
)
