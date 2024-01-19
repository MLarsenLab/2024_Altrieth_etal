
if (!require("pacman")) install.packages("pacman")
pacman::p_load(SingleCellExperiment, celda, Seurat, ggplot2, cowplot, scater, patchwork, 
               dplyr, sctransform, glmGamPoi, stringr, scCustomize, scDblFinder, singleCellTK, biomaRt,
               clustree, hdf5r, rhdf5, sctransform, remotes)

devtools::install_github("sqjin/CellChat")
library("CellChat")
install_github("chris-mcginnis-ucsf/DoubletFinder")
library(DoubletFinder)

dataPath <- "c:/cellbenderOutput"
projectPath <- "c:/JDRFigures"

setwd(projectPath)
set.seed(12354)
options(future.globals.maxSize = 10 * 1024^3)
indexList <- list(mock = 1, ligated = 2, deligated = 3, immune = 1, endothelial = 2, fibroblast = 3, 
                  epithelial = 4, pericyte = 5, glia = 6)

Read_CellBender_h5_Mat <- function(file_name, use.names = TRUE, unique.features = TRUE) {
  #Loads a cellbender output file and returns a matrix that 10x can use to create a seurat object
  #https://rdrr.io/cran/scCustomize/man/Read_CellBender_h5_Mat.html
  # Check hdf5r installed
  if (!requireNamespace('hdf5r', quietly = TRUE)) {
    cli_abort(message = c("Please install hdf5r to read HDF5 files",
                          "i" = "`install.packages('hdf5r')`")
    )
  }
  # Check file
  if (!file.exists(file_name)) {
    stop("File not found")
  }
  
  if (use.names) {
    feature_slot <- 'features/name'
  } else {
    feature_slot <- 'features/id'
  }
  
  # Read file
  infile <- hdf5r::H5File$new(filename = file_name, mode = "r")
  
  counts <- infile[["matrix/data"]]
  indices <- infile[["matrix/indices"]]
  indptr <- infile[["matrix/indptr"]]
  shp <- infile[["matrix/shape"]]
  features <- infile[[paste0("matrix/", feature_slot)]][]
  barcodes <- infile[["matrix/barcodes"]]
  
  
  sparse.mat <- sparseMatrix(
    i = indices[] + 1,
    p = indptr[],
    x = as.numeric(x = counts[]),
    dims = shp[],
    repr = "T"
  )
  
  if (unique.features) {
    features <- make.unique(names = features)
  }
  
  rownames(x = sparse.mat) <- features
  colnames(x = sparse.mat) <- barcodes[]
  sparse.mat <- as(object = sparse.mat, Class = "dgCMatrix")
  
  infile$close_all()
  
  return(sparse.mat)
}
removeDoublets <- function(seuratObject){
  #uses scDblFinder and DoubletFinder to find droplets that may contain mutiple cells
  #Returns a seruat object where cells that both scDblFinder and DoubletFinder agree on as doublets have been removed
  #debug assignments
  if(0){
    seuratObject <- aSeurat
  }
  
  #scDblFinder needs a single cell experement as input and outputs a score and class (doublet or singlet)
  #Assume 20% doublets, this is an over estimation of number of doublets as we are looking for agreement 
  aSCE <- as.SingleCellExperiment(seuratObject)
  aSCE <- scDblFinder(aSCE, dbr = 0.2)
  
  #copy the cell asignemtns to new metadata in the seurat object
  seuratObject$scDblFinder.score <- aSCE$scDblFinder.score
  seuratObject$scDblFinder.class <- aSCE$scDblFinder.class
  
  #DoubletFinder requires data to be normalized before it can identify doublets
  #The object passed in may not be ready to be normalized or may need to be normalized in a differt way so
  #a copy is made and the copy is normalized for DoubletFinder and the cell assignments copied to the original
  dblSeurat <- seuratObject
  dblSeurat <- NormalizeData(dblSeurat)
  dblSeurat <- FindVariableFeatures(dblSeurat, selection.method = "vst", nfeatures = 2000)
  geneList <- rownames(dblSeurat)
  dblSeurat <- ScaleData(dblSeurat, feature = geneList)
  dblSeurat <- RunPCA(dblSeurat, features = VariableFeatures(object = dblSeurat))
  dblSeurat <- RunUMAP(dblSeurat, dims = 1:40)
  
  #DoubletFinder needs to know how many doublets to look for, using 20% of total number of cells
  dblCount <- as.integer(length(rownames(dblSeurat@meta.data))*.2)
  dblSeurat <- doubletFinder_v3(dblSeurat, PCs = 1:40, pN = 0.25, pK = 0.1, nExp = dblCount, reuse.pANN = FALSE, sct = FALSE)
  #DoubletFinder returns the classification in a new metadata with the name based on the parameters used
  DF <- sprintf("DF.classifications_0.25_0.1_%i", dblCount)
  DFIndex <- which(colnames(dblSeurat@meta.data) == DF)
  
  #Copy cell assigment to the seurat object metadata
  seuratObject$DF.classifications <- dblSeurat@meta.data$DF
  
  #Find cells that are clasified as a doublet by both scDblFinder and doubletFinder
  doublets1 <- subset(seuratObject, scDblFinder.class == "doublet")
  trueDoublets <- subset(doublets1, DF.classifications == "Doublet")
  seuratObject$True_Doublet <- ifelse(colnames(seuratObject)%in% colnames(trueDoublets), "Doublet", "Singlet")
  
  #Print the percentage of doublets found
  print(sprintf("%i doublets found in %i total cells (%.2f%%)", ncol(trueDoublets), ncol(seuratObject), 
                (ncol(trueDoublets)/ncol(seuratObject)*100)))
  
  #Make a subset of only singlets and return the cleaned object
  seuratObject <- subset(seuratObject, True_Doublet == "Singlet")
  return(seuratObject)
}

loadForIntegrate <- function(fileNamePath, cellBenderP = TRUE, projectName = "Project"){
  #Create a seurat object from a data file and do basic processing of the seurat object before retruning the object
  #Accepts either cellranger output file or cellbender output file, defaults to cellbender output
  #debug assignments
  if(0)
  {
    fileNamePath <- "C:/outs/filtered_feature_bc_matrix"
    cellBenderP <- FALSE
    projectName <- "Project"
  }
  
  if(cellBenderP){
    aMatrix <- Read_CellBender_h5_Mat(fileNamePath)
    aName <- str_sub(basename(fileNamePath), 1, (str_length(basename(fileNamePath))-12))
    aSeurat <- CreateSeuratObject(aMatrix, project = aName, min.cells = 3, min.genes = 200)
  }else{
    a10X <- Read10X(data.dir = fileNamePath)
    aName <- projectName
    aSeurat <- CreateSeuratObject(a10X, project = aName, min.cells = 3, min.genes = 200)
  }
  #find mitochondria gene percentage and filter on to few (emtpy) to many (doubletts) features and over 25% mitochondria genes
  aSeurat[["percent.mt"]] <- PercentageFeatureSet(aSeurat, pattern = "^mt-")
  aSeurat <- subset(aSeurat, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
  
  #remove doublets and return processed seurat object
  aSeurat <- removeDoublets(aSeurat)
  return(aSeurat)
}
cellChatFromSeurat <- function(aSeurat, protocol){
  #Create a cellchat object from a seurat object
  #debug assigmnet
  if(0)
  {
    aSeurat = aSeuratIntegrated
    protocol <- "All"
  }
  print(aSeurat)
  print(sprintf("Createing cellChat for %s", protocol))
  data.input <- GetAssayData(aSeurat, assay = "SCT", slot = "data")
  labels <- Idents(aSeurat)
  meta <- data.frame(group = labels, row.names = names(labels))
  
  acellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
  
  acellchat <- addMeta(acellchat, meta = meta) 
  acellchat <- setIdent(acellchat, ident.use = "group") # set "group" as default cell identity
  levels(acellchat@idents) # show factor levels of the cell labels
  groupSize <- as.numeric(table(acellchat@idents)) # number of cells in each cell group
  
  CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data

  # use all CellChatDB for cell-cell communication analysis
  CellChatDB.use <- CellChatDB # simply use the default CellChatDB
  
  #Fixes random error that cellchat has with its own DB
  #https://github.com/sqjin/CellChat/issues/45
  removeID <- which(CellChatDB.use[["interaction"]]$ligand == "H2-BI") # 1887
  CellChatDB.use[["interaction"]] <- CellChatDB.use[["interaction"]][-1887,]
  removeID <- which(CellChatDB.use[["interaction"]]$ligand == "H2-Ea-ps") #1900
  CellChatDB.use[["interaction"]] <- CellChatDB.use[["interaction"]][-removeID,]
  
  acellchat@DB <- CellChatDB.use
  
  # subset the expression data of signaling genes for saving computation cost
  acellchat <- subsetData(acellchat) # This step is necessary even if using the whole database
  future::plan("multisession", workers = 8) # do parallel
  
  acellchat <- identifyOverExpressedGenes(acellchat)
  acellchat <- identifyOverExpressedInteractions(acellchat)
  #acellchat <- projectData(acellchat, PPI.mouse)   #possable cleanup step, optional, need to use raw.use = FALSE in computeComunProb()
  
  acellchat <- computeCommunProb(acellchat)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  acellchat <- filterCommunication(acellchat, min.cells = 10)
  
  acellchat <- computeCommunProbPathway(acellchat)
  
  acellchat <- aggregateNet(acellchat)
  
  # Compute the network centrality scores
  # the slot 'netP' means the inferred intercellular communication network of signaling pathways
  acellchat <- netAnalysis_computeCentrality(acellchat, slot.name = "netP") 
  #Sets 'communication patterns'
  acellchat <- identifyCommunicationPatterns(acellchat, pattern = "outgoing", k = 3)
  acellchat <- identifyCommunicationPatterns(acellchat, pattern = "incoming", k = 3)
  
  #add metadata for the cellchat object and return the cellchat object
  acellchat@meta$protocol <- protocol
  
  return(acellchat)
}
saveClusterMarkers <- function(seuratObject, projName = "project"){
  #Save all markers, positive gene averages and cluster cell count for a seuratObject
  #Debug assignment
  if(0){
    seuratObject = aSeurat
    projName = "cellbenderSeurat"
  }
  
  dir.create(projName)
  
  apath <- sprintf("%s/geneAverages_%s.csv", projName, projName)
  write.csv(FindAllMarkers(seuratObject, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25), file = apath)
  apath <- sprintf("%s/positveGeneAverages_%s.csv", projName, projName)           
  write.csv(FindAllMarkers(seuratObject, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25), file = apath)
  apath <- sprintf("%s/cellsPerCluster%s.csv", projName, projName)
  write.csv(table(Idents(seuratObject)), file = apath)
}
netVisual_bubble_MLD <- function(aCellChat, source = c(1), target = c(1), apath = "netVisualBubble") {
  #Outup netVisual bubble for Mock v Ligated, Mock v Deligated and Deligated v Ligated samples
  
  netVisual_bubble_invert(aCellChat, source, target, comparison = c(indexList$mock, indexList$ligated),
                          label = c("Mock", "Ligated"), apath = apath)
  netVisual_bubble_invert(aCellChat, source, target, comparison = c(indexList$mock, indexList$deligated),
                          label = c("Mock", "Deligated"), apath = apath)
  netVisual_bubble_invert(aCellChat, source, target, comparison = c(indexList$deligated, indexList$ligated),
                          label = c("Deligated", "Ligated"), apath = apath)
}
netVisual_bubble_invert <- function(aCellChat, source = c(1), target = c(1), comparison, label, apath = "netVisualBubble") {
  #Generate the netVisual bubble for a v b and b v a.
  
  dir.create(apath)
  dir.create(sprintf("%s/To_%s", apath, levels(aCellChat@meta$group)[source[1]]))
  dir.create(sprintf("%s/To_%s", apath, levels(aCellChat@meta$group)[target[1]]))
  
  fileNamePath <- sprintf("%s/To_%s/netVisual_bubble_%sv%s_From_%s_to_%s.pdf", apath, 
                          levels(aCellChat@meta$group)[target[1]], label[1], label[2], 
                          levels(aCellChat@meta$group)[source[1]], levels(aCellChat@meta$group)[target[1]])
  pdf(fileNamePath)
  #If there are no pathways in the comparison, then netVisual_bubble halts, this will prevent the function(s) from halting prematurely
  tryCatch({
    print(netVisual_bubble(aCellChat, sources.use = source, targets.use = target,  
                           comparison = comparison, 
                           max.dataset = comparison[1], 
                           title.name = sprintf("Increased signaling in %s vs %s", label[1], label[2]), 
                           angle.x = 45, 
                           remove.isolate = T))
  }, error=function(e){cat(fileNamePath, conditionMessage(e), "\n")})
  dev.off()
  
  fileNamePath <- sprintf("%s/To_%s/netVisual_bubble_%sv%s_From_%s_to_%s.pdf", apath, 
                          levels(aCellChat@meta$group)[target[1]], label[2], label[1], 
                          levels(aCellChat@meta$group)[source[1]], levels(aCellChat@meta$group)[target[1]])
  pdf(fileNamePath)
  tryCatch({
    print(netVisual_bubble(aCellChat, sources.use = source, targets.use = target,  
                           comparison = comparison, 
                           max.dataset = comparison[2], 
                           title.name = sprintf("Increased signaling in %s vs %s", label[2], label[1]), 
                           angle.x = 45, 
                           remove.isolate = T))
  }, error=function(e){cat(fileNamePath ,conditionMessage(e), "\n")})
  dev.off()
}
convert_human_to_mouse <- function(gene_list) {
  #Adapted from https://www.biostars.org/p/9567892/
  #Converts a list of human gene names to a list of mouse gene equivalents
  output = c()
  mouse_human_genes = read.csv("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
  
  for(gene in gene_list) {
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name == "human"))[['DB.Class.Key']]
    if( !identical(class_key, integer(0)) ) {
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="mouse, laboratory"))[,"Symbol"]
      for(human_gene in human_genes) {
        output = rbind(c(gene, human_gene), output)
      }
    }
  }
  return (output[,2])
}

#Cellbender cleanup, all used datasets SCTransform integrated----------------------------
seuratList <- c()
aName <- "ligated_14_repeat"
afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
aSeurat <- loadForIntegrate(afile)
aSeurat@meta.data[, "protocol"] <- "Ligated"
seuratList <- c(seuratList, aSeurat)

aName <- "mock_14_repeat"
afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
aSeurat <- loadForIntegrate(afile)
aSeurat@meta.data[, "protocol"] <- "Mock"
seuratList <- c(seuratList, aSeurat)

aName <- "mock_14_new1"
afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
aSeurat <- loadForIntegrate(afile)
aSeurat@meta.data[, "protocol"] <- "Mock"
seuratList <- c(seuratList, aSeurat)

aName <- "mock_14_new2"
afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
aSeurat <- loadForIntegrate(afile)
aSeurat@meta.data[, "protocol"] <- "Mock"
seuratList <- c(seuratList, aSeurat)

aName <- "deligated"
afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
aSeurat <- loadForIntegrate(afile)
aSeurat@meta.data[, "protocol"] <- "Deligated"
seuratList <- c(seuratList, aSeurat)

aName <- "deligated_repeat"
afile <- sprintf("%s/%s/%s_CellBender_filtered.h5", dataPath, aName, aName)
aSeurat <- loadForIntegrate(afile)
aSeurat@meta.data[, "protocol"] <- "Deligated"
seuratList <- c(seuratList, aSeurat)

seuratList <- lapply(X = seuratList, FUN = function(x) {
  x <- SCTransform(x, vst.flavor = "v2", verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = seuratList, nfeatures = 3000)
seuratList <- PrepSCTIntegration(object.list = seuratList, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = seuratList, anchor.features = features, normalization.method = "SCT")
aSeuratIntegrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
DefaultAssay(aSeuratIntegrated) <- "integrated"

aSeuratIntegrated <- ScaleData(aSeuratIntegrated, verbose = FALSE)
aSeuratIntegrated <- RunPCA(aSeuratIntegrated, npcs = 40, verbose = FALSE)
aSeuratIntegrated <- RunUMAP(aSeuratIntegrated, reduction = "pca", dims = 1:40)
aSeuratIntegrated <- FindNeighbors(aSeuratIntegrated, reduction = "pca", dims = 1:40)

clusterResolution <- 2
seuratTree <- FindClusters(aSeuratIntegrated, resolution = (seq(0, (clusterResolution*2), by = (clusterResolution/5))))
pdf("All_integrated_clustree.pdf", width = 40, height = 30)
print(clustree(seuratTree, prefix = "integrated_snn_res.", node_colour = "sc3_stability"))
dev.off()

aSeuratIntegrated <- FindClusters(aSeuratIntegrated, resolution = clusterResolution)

Idents(aSeuratIntegrated) <- "seurat_clusters"
levels(aSeuratIntegrated)
pdf("UMAP_integrate_All.pdf", width = 40, height = 30)
print(DimPlot(aSeuratIntegrated, reduction = "umap", label = TRUE, label.size = 16, pt.size = 0.5, repel = TRUE) +
        #NoLegend() +
        theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.text = element_text(size = 40)))
dev.off()

Idents(aSeuratIntegrated) <- "protocol"
levels(aSeuratIntegrated)
pdf("UMAP_integrate_All_protocol.pdf", width = 40, height = 30)
print(DimPlot(aSeuratIntegrated, reduction = "umap", label = TRUE, label.size = 16, pt.size = 0.5, repel = TRUE) +
        #NoLegend() +
        theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.text = element_text(size = 40)))
dev.off()

DefaultAssay(aSeuratIntegrated) <- "SCT"
Idents(aSeuratIntegrated) <- "seurat_clusters"
aSeuratIntegrated <- PrepSCTFindMarkers(aSeuratIntegrated, assay = "SCT", verbose = TRUE)

saveRDS(aSeuratIntegrated, file = "JDR2023Seurat-2.rds")
aSeuratIntegrated <- readRDS("JDR2023Seurat-2.rds")

saveClusterMarkers(aSeuratIntegrated, projName = "ClusterMarkers")

#Metadata setupfor aSeuratIntegrated------------------------------
#Load seurat object of Integrated datasets created in E034 and clutered with a resolution of 2
aSeuratIntegrated <- readRDS(sprintf("%s/JDR2023Seurat-2.rds", dataPath))

#Print UMAP of all the integrated data
Idents(aSeuratIntegrated) <- "seurat_clusters"
levels(aSeuratIntegrated)
pdf("UMAP_integrate_All.pdf", width = 40, height = 30)
print(DimPlot(aSeuratIntegrated, reduction = "umap", label = TRUE, label.size = 16, pt.size = 0.5, repel = TRUE) +
        #NoLegend() +
        theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.text = element_text(size = 40)))
dev.off()

#Reorder the protocols for beter visulization
Idents(aSeuratIntegrated) <- "protocol"
levels(aSeuratIntegrated)
aSeuratIntegrated$protocol <- factor(aSeuratIntegrated$protocol, levels = c("Mock", "Ligated", "Deligated"))
Idents(aSeuratIntegrated) <- "protocol"
levels(aSeuratIntegrated)
#Print UMAP broken down by protocol
pdf("UMAP_integrate_All_protocol.pdf", width = 40, height = 30)
print(DimPlot(aSeuratIntegrated, reduction = "umap", label = FALSE, label.size = 16, pt.size = 0.5, repel = TRUE) +
        #NoLegend() +
        theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.text = element_text(size = 120)))
dev.off()

#Rename cluster IDs to cell type
Idents(aSeuratIntegrated) <- "seurat_clusters"
levels(aSeuratIntegrated)
newClusterIDsLong <- c("0_Immune", "1_Endothelial", "2_Endothelial", "3_Endothelial", "4_Endothelial", "5_Immune", 
                       "6_Endothelial", "7_Immune", "8_Immune", "9_Fibroblast", "10_Endothelial", 
                       "11_Endothelial", "12_Immune", "13_Endothelial", "14_Immune", "15_Endothelial", 
                       "16_Endothelial", "17_Endothelial", "18_Epithelial", "19_Fibroblast", "20_Pericyte", 
                       "21_Endothelial", "22_Immune", "23_Endothelial", "24_Pericyte", "25_Pericyte",
                       "26_Epithelial", "27_Immune", "28_Endothelial", "29_Glia", "30_Immune",
                       "31_Epithelial", "32_Immune", "33_Immune", "34_Endothelial", "35_Endothelial",
                       "36_Endothelial", "37_Immune", "38_Endothelial", "39_Epithelial", "40_Endothelial",
                       "41_Immune", "42_Immune", "43_Immune", "44_Fibroblast", "45_Immune",
                       "46_Fibroblast", "47_Endothelial", "48_Immune")
newClusterIDsShort <- c("Immune", "Endothelial", "Endothelial", "Endothelial", "Endothelial", "Immune", 
                        "Endothelial", "Immune", "Immune", "Fibroblast", "Endothelial", 
                        "Endothelial", "Immune", "Endothelial", "Immune", "Endothelial", 
                        "Endothelial", "Endothelial", "Epithelial", "Fibroblast", "Pericyte", 
                        "Endothelial", "Immune", "Endothelial", "Pericyte", "Pericyte",
                        "Epithelial", "Immune", "Endothelial", "Glia", "Immune",
                        "Epithelial", "Immune", "Immune", "Endothelial", "Endothelial",
                        "Endothelial", "Immune", "Endothelial", "Epithelial", "Endothelial",
                        "Immune", "Immune", "Immune", "Fibroblast", "Immune",
                        "Fibroblast", "Endothelial", "Immune")
Idents(aSeuratIntegrated) <- "seurat_clusters"
newClusterIDs <- newClusterIDsShort
names(newClusterIDs) <- levels(aSeuratIntegrated)
aSeuratIntegrated <- RenameIdents(aSeuratIntegrated, newClusterIDs)
levels(aSeuratIntegrated)
aSeuratIntegrated$ID <- aSeuratIntegrated@active.ident

#Print UMAP using the new cluster IDs
pdf("UMAP_integrate_All_NewID.pdf", width = 40, height = 30)
print(DimPlot(aSeuratIntegrated, reduction = "umap", label = TRUE, label.size = 40, pt.size = 0.5, repel = TRUE) +
        NoLegend() +
        theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.text = element_text(size = 40)))
dev.off()

#Save protocl + cluster as a new metadata field
aSeuratIntegrated$protocol_cluster <- paste0(aSeuratIntegrated$protocol, "_", aSeuratIntegrated$ID)
#Reorder the lables for better figures
Idents(aSeuratIntegrated) <- "protocol_cluster"
levels(aSeuratIntegrated) <- c("Deligated_Endothelial", "Ligated_Endothelial", "Mock_Endothelial", 
                               "Deligated_Epithelial", "Ligated_Epithelial", "Mock_Epithelial", 
                               "Deligated_Fibroblast", "Ligated_Fibroblast", "Mock_Fibroblast", 
                               "Deligated_Glia", "Ligated_Glia", "Mock_Glia", 
                               "Deligated_Immune", "Ligated_Immune", "Mock_Immune",
                               "Deligated_Pericyte", "Ligated_Pericyte", "Mock_Pericyte")
levels(aSeuratIntegrated)

#Save the seurat object with the metadata fields added
saveRDS(aSeuratIntegrated, file = sprintf("%s/JDR2023Seurat-Meta.rds", dataPath))
aSeuratIntegrated <- readRDS(sprintf("%s/JDR2023Seurat-Meta.rds", dataPath))


#Generate the files for NATMI, NATMI needs to be run on the commandline-----------------------
#NATMI needs a data file of counts and a metadata file with cluster ID per cell
#Create the files for each of the 3 conditions

aSeuratIntegrated <- readRDS(sprintf("%s/JDR2023Seurat-Meta.rds", dataPath))
seuratList <- c()
seuratList <- c(seuratList, subset(aSeuratIntegrated, protocol == "Mock"))
seuratList <- c(seuratList, subset(aSeuratIntegrated, protocol == "Ligated"))
seuratList <- c(seuratList, subset(aSeuratIntegrated, protocol == "Deligated"))

for(i in seuratList){
  #From NATMI
  write.csv(100 * (exp(as.matrix(GetAssayData(object = i, assay = "SCT", slot = "data"))) - 1), 
            sprintf("data-%s.csv", i@meta.data$protocol[1]), row.names = T)
  write.csv(Idents(object = i), sprintf("metadata-%s.csv", i@meta.data$protocol[1]), row.names = T)
  
}

#Generate cellchat objects---------------------------------------
#Generates cellchat object for each of the 3 conditions and a merged cell chat object and does the pre-processing for the
#merged object

aSeuratIntegrated <- readRDS(sprintf("%s/JDR2023Seurat-Meta.rds", dataPath))
seuratList <- c()
seuratList <- c(seuratList, subset(aSeuratIntegrated, protocol == "Mock"))
seuratList <- c(seuratList, subset(aSeuratIntegrated, protocol == "Ligated"))
seuratList <- c(seuratList, subset(aSeuratIntegrated, protocol == "Deligated"))

#Create the cellchat object, save it as a .rds file with the name Cellchat<protocl> and save it in a list
cellchatList <- c()
for(i in 1:3){
  aCellChat <- cellChatFromSeurat(seuratList[[i]], seuratList[[i]]$protocol[1])
  saveRDS(aCellChat, file = sprintf("CellChat%s.rds", seuratList[[i]]$protocol[1]))
  cellchatList <- c(cellchatList, aCellChat)
}

aCellchatMerge <- mergeCellChat(cellchatList, cell.prefix = TRUE, add.names = c(cellchatList[[1]]@meta$protocol[1], 
                                                            cellchatList[[2]]@meta$protocol[1], 
                                                            cellchatList[[3]]@meta$protocol[1]))

aCellchatMerge <- computeNetSimilarityPairwise(aCellchatMerge, type = "structural")
aCellchatMerge <- netEmbedding(aCellchatMerge, type = "structural", umap.method = "uwot")
aCellchatMerge <- netClustering(aCellchatMerge, type = "structural", do.parallel = FALSE)

aCellchatMerge <- computeNetSimilarityPairwise(aCellchatMerge, type = "functional")
aCellchatMerge <- netEmbedding(aCellchatMerge, type = "functional", umap.method = "uwot")
aCellchatMerge <- netClustering(aCellchatMerge, type = "functional", do.parallel = FALSE)

#Save the merged cellchat object
saveRDS(aCellchatMerge, file = "CellChatMerge.rds")

#Differanl gene expression---------------------

Idents(aSeuratIntegrated) <- "protocol_cluster"
levels(aSeuratIntegrated)

apath <- "Pos_Markers_Mock_Endo_Ligated_Endo.csv"
write.csv(FindMarkers(aSeuratIntegrated, ident.1 = "Mock_Endothelial", ident.2 = "Ligated_Endothelial", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25), file = apath)
apath <- "Pos_Markers_Ligated_Endo_Mock_Endo.csv"
write.csv(FindMarkers(aSeuratIntegrated, ident.1 = "Ligated_Endothelial", ident.2 = "Mock_Endothelial", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25), file = apath)
apath <- "Pos_Markers_Mock_Endo_Deligated_Endo.csv"
write.csv(FindMarkers(aSeuratIntegrated, ident.1 = "Mock_Endothelial", ident.2 = "Deligated_Endothelial", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25), file = apath)
apath <- "Pos_Markers_Deligated_Endo_Mock_Endo.csv"
write.csv(FindMarkers(aSeuratIntegrated, ident.1 = "Deligated_Endothelial", ident.2 = "Mock_Endothelial", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25), file = apath)
apath <- "Pos_Markers_Ligated_Endo_Deligated_Endo.csv"
write.csv(FindMarkers(aSeuratIntegrated, ident.1 = "Ligated_Endothelial", ident.2 = "Deligated_Endothelial", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25), file = apath)
apath <- "Pos_Markers_Deligated_Endo_Ligated_Endo.csv"
write.csv(FindMarkers(aSeuratIntegrated, ident.1 = "Deligated_Endothelial", ident.2 = "Ligated_Endothelial", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25), file = apath)


#Netvisual bubble show signaling increase from cluster to cluster between conditions (Mock, Ligated, Deligated)---------------------------------
#From Endotheal cells (cluster 2)

aCellchatMerge <- readRDS("CellChatMerge.rds")
levels(aCellchatMerge@meta$group)

#Create netVisual bubbles from endothelial cells to all clusters in mock v ligated, ligated v mock, mock v deligated, deligated v mock, ligated v deligated, deligated v ligated
for(i in 1:nlevels(aCellchatMerge@meta$group)){
  netVisual_bubble_MLD(aCellchatMerge, source = 2, target = i, apath = "netVisualBubble")
}

#General cellchat figures---------------------------

cellchatList <- c()
aCellChat <- readRDS("CellChatMock.rds")
cellchatList <- c(cellchatList, aCellChat)
aCellChat <- readRDS("CellChatLigated.rds")
cellchatList <- c(cellchatList, aCellChat)
aCellChat <- readRDS("CellChatDeligated.rds")
cellchatList <- c(cellchatList, aCellChat)
aCellchatMerge <- readRDS("CellChatMerge.rds")

groupSize <- as.numeric(table(aCellChat@idents))

pdf("netVisual_circle_End_weight_mock.pdf")
netVisual_circle(cellchatList[[indexList$mock]]@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions", sources.use = "Endothelial")
dev.off()

pdf("netVisual_circle_End_weight_ligated.pdf")
netVisual_circle(cellchatList[[indexList$ligated]]@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions", sources.use = "Endothelial")
dev.off()

pdf("netVisual_circle_End_weight_deligated.pdf")
netVisual_circle(cellchatList[[indexList$deligated]]@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions", sources.use = "Endothelial")
dev.off()

pdf("netVisual_circle_End_count_mock.pdf")
netVisual_circle(cellchatList[[indexList$mock]]@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions", sources.use = "Endothelial")
dev.off()

pdf("netVisual_circle_End_count_ligated.pdf")
netVisual_circle(cellchatList[[indexList$ligated]]@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions", sources.use = "Endothelial")
dev.off()

pdf("netVisual_circle_End_count_deligated.pdf")
netVisual_circle(cellchatList[[indexList$deligated]]@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions", sources.use = "Endothelial")
dev.off()

#List of all pathways in all conditions
pathway.union <- union(cellchatList[[indexList$mock]]@netP$pathways, 
                       cellchatList[[indexList$ligated]]@netP$pathways)
pathway.union <- union(pathway.union, 
                       cellchatList[[indexList$deligated]]@netP$pathways)

pdf("netAnalysis_signalingRole_heatmap_mock.pdf", height = 9)
  print(netAnalysis_signalingRole_heatmap(cellchatList[[indexList$mock]], pattern = "all", 
                                        signaling = pathway.union, title = "Mock", 
                                        width = 10, height = 18, font.size = 8, color.heatmap = "OrRd"))
dev.off()
pdf("netAnalysis_signalingRole_heatmap_ligated.pdf", height = 9)
  print(netAnalysis_signalingRole_heatmap(cellchatList[[indexList$ligated]], pattern = "all", 
                                        signaling = pathway.union, title = "Ligated", 
                                        width = 10, height = 18, font.size = 8, color.heatmap = "OrRd"))
dev.off()
pdf("netAnalysis_signalingRole_heatmap_deligated.pdf", height = 9)
  print(netAnalysis_signalingRole_heatmap(cellchatList[[indexList$deligated]], pattern = "all", 
                                        signaling = pathway.union, title = "Deligated", 
                                        width = 10, height = 18, font.size = 8, color.heatmap = "OrRd"))
dev.off()

pdf("compareInteractions_count.pdf")
compareInteractions(aCellchatMerge, show.legend = F, group = c(1,2,3), width = .2)
dev.off()

pdf("compareInteractions_strength.pdf")
compareInteractions(aCellchatMerge, show.legend = F, group = c(1,2,3), measure = "weight", width = .2, digits=2)
dev.off()

#reordering the clusters so that netVisual_diffInteraction works. netVisual_diffInteraction is bugged and will only work for clusters 1 and <max>
tmp <- c()
aCellChat <- updateClusterLabels(cellchatList[[1]], new.order = c("Endothelial", "Immune", "Fibroblast", "Epithelial", "Pericyte", "Glia"))
tmp <- c(tmp, aCellChat)
aCellChat <- updateClusterLabels(cellchatList[[2]], new.order = c("Endothelial", "Immune", "Fibroblast", "Epithelial", "Pericyte", "Glia"))
tmp <- c(tmp, aCellChat)
aCellChat <- updateClusterLabels(cellchatList[[3]], new.order = c("Endothelial", "Immune", "Fibroblast", "Epithelial", "Pericyte", "Glia"))
tmp <- c(tmp, aCellChat)

Mergetmp <- mergeCellChat(tmp, cell.prefix = TRUE, add.names = c(tmp[[1]]@meta$protocol[1], 
                                                                                tmp[[2]]@meta$protocol[1], 
                                                                                tmp[[3]]@meta$protocol[1]))

Mergetmp <- computeNetSimilarityPairwise(Mergetmp, type = "structural")
Mergetmp <- netEmbedding(Mergetmp, type = "structural", umap.method = "uwot")
Mergetmp <- netClustering(Mergetmp, type = "structural", do.parallel = FALSE)

Mergetmp <- computeNetSimilarityPairwise(Mergetmp, type = "functional")
Mergetmp <- netEmbedding(Mergetmp, type = "functional", umap.method = "uwot")
Mergetmp <- netClustering(Mergetmp, type = "functional", do.parallel = FALSE)

pdf("cellchat/netVisual_diffInteraction_mockvligated.pdf")
netVisual_diffInteraction(Mergetmp, comparison = c(indexList$mock, indexList$ligated), measure = "weight",
                          sources.use = 1)
dev.off()
pdf("cellchat/netVisual_diffInteraction_ligatedvmock.pdf")
netVisual_diffInteraction(Mergetmp, comparison = c(indexList$ligated, indexList$mock), measure = "weight",
                          sources.use = 1)
dev.off()
pdf("cellchat/netVisual_diffInteraction_mockvdeligated.pdf")
netVisual_diffInteraction(Mergetmp, comparison = c(indexList$mock, indexList$deligated), measure = "weight",
                          sources.use = 1)
dev.off()
pdf("cellchat/netVisual_diffInteraction_deligatedvmock.pdf")
netVisual_diffInteraction(Mergetmp, comparison = c(indexList$deligated, indexList$mock), measure = "weight",
                          sources.use = 1)
dev.off()
pdf("cellchat/netVisual_diffInteraction_ligatedvdeligated.pdf")
netVisual_diffInteraction(Mergetmp, comparison = c(indexList$ligated, indexList$deligated), measure = "weight",
                          sources.use = 1)
dev.off()
pdf("cellchat/netVisual_diffInteraction_deligatedvligated.pdf")
netVisual_diffInteraction(Mergetmp, comparison = c(indexList$deligated, indexList$ligated), measure = "weight",
                          sources.use = 1)
dev.off()

#EndMT and Nolan factors dotplots------------------------------------
#Nolan et al. 2013 angiocrine factors (PMID: 23871589)

aSeuratIntegrated <- readRDS(sprintf("%s/JDR2023Seurat-Meta.rds", dataPath))
Idents(aSeuratIntegrated) <- aSeuratIntegrated$ID
levels(aSeuratIntegrated)

NolanList <- c("Bmp2", "Ccl2", "Ccl3", "Ccl6", "Ccl9", "Col14a1", "Col4a1", "Cxcl12", "Cxcl13", "Cxcl2", "Dkk2", 
               "Edn1", "Egfl7", "Eln", "Esm1", "Fgf7", "Fn1", "Gja1", "Gja4", "Igf1", "Il1a", "Il33", "Jag1", 
               "Kitl", "Lama4", "Lamb1", "Lamc1", "Mmp13", "Mmp9", "Pdgfa", "Tgfb2", "Tgfb3", "Tnf", "Vegfc")

pdf("Nolan-Endo-DotPlot.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratIntegrated, features = NolanList, cols = c("lightblue", "red3"), 
              col.min = 0, dot.scale = 20, idents = c("Endothelial"), group.by = "protocol") +
        theme(axis.title = element_text(size = 70), axis.text = element_text(size = 70), 
              legend.text = element_text(size = 70), legend.key.size = unit(2, "cm"), 
              legend.title  = element_text(size = 70), legend.spacing = unit(2, "cm"),
              axis.text.x = element_text(angle=90, vjust = 0.45)))
dev.off()

#ENdoMT list(s)
#Emerging roles of inflammation-mediated endothelial–mesenchymal transition in health and disease
#Yasuhiro Yoshimatsu & Tetsuro Watabe (2022)
EndoMT_Inducer <- c("Bmp9", "Mmp14", "Tgfb1", "Tgfb2")
EndoMT_Inflammation <- c("C1qa", "C1qb", "C5ar1", "Ccl2", "Ccl3", "Cxcl2", "Icam1", "Ifng", "Il1b", "Il4", "Il6", 
                         "Jak2", "Rela", "Stat3", "Tnf", "Vcam1")
EndoMT_EC_Gene <- c("Apln", "Cdh5", "Erg", "Fgfr1", "Flt4", "Fstl3", "Kdr", "H19", "Lyve1", 
                    "Nos3", "Pdpn", "Pecam1", "Prox1", "Tal1", "Tie1", 
                    "Tek", "Vwf", "Wnt5a")

EndoMT_Mesen_Gene <- c("Acta2", "Cdh2", "Cnn1", "Col12a1", "Col1a1", "Col1a2", "Col3a1", 
                       "Col6a1", "Ctgf", "Fap", "Fn1", "Lamc2", "Ly6a", "Atxn1", "Mmp14", "Mmp2", "Mmp9", "Notch3", 
                       "Pcolce", "Pdgfrb", "S100a4", "Aifm2", "Serpine1", "Tagln", "Tgfbi", "Tnc", "Vim")

EndoMT_Mesen_Gene_short <- c("Acta2", "Col12a1", "Col1a1", "Col1a2", "Col3a1", 
                       "Col6a1", "Fn1", "Lamc2", "Ly6a", "Atxn1", "Notch3", 
                       "Pdgfrb", "S100a4", "Vim")

EndoMT_TF <- c("Hes1", "Hey1", "Hey2", "Mkl1", "Mrtfa", "Snai1", "Snai2", "Twist1", "Twist2", "Zeb1", "Zeb2")

EndoMT_subset <- c("Tgfb1", "Tgfb2",  #Inducer
                   "Icam1", "Jak2", "Rela", "Stat3", #Inflammation
                   "Cdh5", "Erg", "Flt4", "Kdr", "Nos3", "Pecam1", "Tal1", "Tie1", "Tek", #EC Genes
                   "Ly6a", "Atxn1", "Vim", #Mes Genes
                   "Hes1", "Hey1", "Snai1", "Snai2", "Zeb1") #Transcription factors

pdf("EndoMT-Endo-Inducer-DotPlot.pdf", width = 40, height = 30)
print(DotPlot(aSeuratIntegrated, features = EndoMT_Inducer, cols = c("lightblue", "red3"), 
              col.min = 0, dot.scale = 60, idents = c("Endothelial"), group.by = "protocol") +
        theme(axis.title = element_text(size = 60), axis.text = element_text(size = 50), 
              legend.text = element_text(size = 40), legend.key.size = unit(3, 'cm'), 
              legend.title  = element_text(size = 40), legend.spacing = unit(3, 'cm')) +
        RotatedAxis())
dev.off()

pdf("EndoMT-Endo-Inflammation-DotPlot.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratIntegrated, features = EndoMT_Inflammation, cols = c("lightblue", "red3"), 
              col.min = 0, dot.scale = 60, idents = c("Endothelial"), group.by = "protocol") +
        theme(axis.title = element_text(size = 100), axis.text = element_text(size = 100), 
              legend.text = element_text(size = 100), legend.key.size = unit(3, "cm"), 
              legend.title  = element_text(size = 100), legend.spacing = unit(3, "cm"), 
              axis.text.x = element_text(angle=90, vjust = 0.45)))
dev.off()

pdf("EndoMT-Endo-EC_Gene-DotPlot.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratIntegrated, features = EndoMT_EC_Gene, cols = c("lightblue", "red3"), 
              col.min = 0, dot.scale = 60, idents = c("Endothelial"), group.by = "protocol") +
        theme(axis.title = element_text(size = 80), axis.text = element_text(size = 100), 
              legend.text = element_text(size = 60), legend.key.size = unit(2, "cm"), 
              legend.title  = element_text(size = 60), legend.spacing = unit(2, "cm"), 
              axis.text.x = element_text(angle=90)))
dev.off()

pdf("EndoMT-Endo-Mesen_Gene_short-DotPlot.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratIntegrated, features = EndoMT_Mesen_Gene_short, cols = c("lightblue", "red3"), 
              col.min = 0, dot.scale = 30, idents = c("Endothelial"), group.by = "protocol") +
        theme(axis.title = element_text(size = 100), axis.text = element_text(size = 100), 
              legend.text = element_text(size = 100), legend.key.size = unit(3, "cm"), 
              legend.title  = element_text(size = 100), legend.spacing = unit(3, "cm"), 
              axis.text.x = element_text(angle=90)))
dev.off()

pdf("EndoMT-Endo-Mesen_Gene-DotPlot.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratIntegrated, features = EndoMT_Mesen_Gene, cols = c("lightblue", "red3"), 
              col.min = 0, dot.scale = 30, idents = c("Endothelial"), group.by = "protocol") +
        theme(axis.title = element_text(size = 50), axis.text = element_text(size = 90), 
              legend.text = element_text(size = 40), legend.key.size = unit(1, "cm"), 
              legend.title  = element_text(size = 50), legend.spacing = unit(1, "cm"), 
              axis.text.x = element_text(angle=90, vjust = 0.45)))
dev.off()

pdf("EndoMT-Endo-TF_Gene-DotPlot.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratIntegrated, features = EndoMT_TF, cols = c("lightblue", "red3"), 
              col.min = 0, dot.scale = 60, idents = c("Endothelial"), group.by = "protocol") +
        theme(axis.title = element_text(size = 100), axis.text = element_text(size = 100), 
              legend.text = element_text(size = 100), legend.key.size = unit(3, "cm"), 
              legend.title  = element_text(size = 100), legend.spacing = unit(3, "cm"),
              axis.text.x = element_text(angle=90, vjust = 0.45)))
dev.off()

pdf("EndoMT-Endo-Subset-DotPlot.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratIntegrated, features = EndoMT_EC_Gene, cols = c("lightblue", "red3"), 
              col.min = 0, dot.scale = 40, idents = c("Endothelial"), group.by = "protocol") +
        theme(axis.title = element_text(size = 100), axis.text = element_text(size = 100), 
              legend.text = element_text(size = 100), legend.key.size = unit(3, "cm"), 
              legend.title  = element_text(size = 100), legend.spacing = unit(3, "cm"),
              axis.text.x = element_text(angle=90, vjust = 0.45)))
dev.off()

#Endothelial subsetting---------------------------

aSeuratEndo <- subset(aSeuratIntegrated, idents = "Endothelial")
DefaultAssay(aSeuratEndo) <- "integrated"

#Re normalize the subsetted data
aSeuratEndo <- SCTransform(aSeuratEndo, vst.flavor = "v2", verbose = FALSE)
aSeuratEndo <- RunPCA(aSeuratEndo, npcs = 40, verbose = FALSE)
aSeuratEndo <- RunUMAP(aSeuratEndo, reduction = "pca", dims = 1:40)
aSeuratEndo <- FindNeighbors(aSeuratEndo, reduction = "pca", dims = 1:40)

#Cluster the subsetted data for endo subset ID
clusterResolution <- 0.5
seuratTree <- FindClusters(aSeuratEndo, resolution = (seq(0, (clusterResolution*2), by = (clusterResolution/5))))

pdf("Endo-clusttree.pdf", width = 40, height = 30)
print(clustree(seuratTree, prefix = "SCT_snn_res."))
dev.off()

aSeuratEndo <- FindClusters(aSeuratEndo, resolution = clusterResolution)
aSeuratEndo <- RunTSNE(aSeuratEndo, dims = 1:40, check_duplicates=FALSE)
aSeuratEndo <- RunUMAP(aSeuratEndo, dims = 1:40)

Idents(aSeuratEndo) <- "seurat_clusters"
levels(aSeuratEndo)
#Print UMAP of Endothelial subcluster clusters
pdf("Endo-umap.pdf", width = 40, height = 30)
print(DimPlot(aSeuratEndo, reduction = "umap", label = TRUE, label.size = 40, pt.size = 0.5, repel = TRUE) +
        NoLegend() +
        theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40)))
dev.off()

DefaultAssay(aSeuratEndo) <- "SCT"
saveClusterMarkers(aSeuratEndo, projName = "Endo-markers")
#Try to ID the Endothelial subclusters based on Arter, Vein, Capillary and Lymph markers  
Artery_Markers_Kaluka <- c("8430408G22Rik", "Clu", "Crip1", "Fbln2", "Gja4", "Hey1", "Mecom", 
                           "Sat1", "Sema3g", "Sox17", "Tm4sf1", "Tsc22d1")

Vein_Markers_Kaluka <- c("Apoe", "Bgn", "Ctla2a", "Icam1", "Il6st", "Ptgs1", "Tmsb10", "Vcam1", "Vwf")

Capillary_Markers_Kaluka <- c("AW112010", "BC028528", "Car4", "Cd200", "Cd300lg", "Gpihbp1", "Kdr", "Rgcc", 
                              "Sgk1", "Sparc")

Lymphatic_Markers_Kaluka <- c("Ccl21a", "Cd63", "Cp", "Fgl2", "Flt4", "Fth1", "Fxyd6", "Maf", "Marcks", "Mmrn1", 
                              "Pard6g", "Pdpn", "Prelp", "Reln", "Rnase4", "Scn1b", "Stab1", "Thy1", "Timp2", "Timp3")         

pdf("Markers_Endo_Artery.pdf", width = 40, height = 30)
print(DotPlot(aSeuratEndo, features = Artery_Markers_Kaluka, cols = c("lightblue", "red3"), 
              col.min = 0, dot.scale = 30) +
          theme(axis.title = element_text(size = 60), axis.text = element_text(size = 50), 
          legend.text = element_text(size = 40), legend.key.size = unit(3, 'cm'), 
          legend.title  = element_text(size = 40), legend.spacing = unit(3, 'cm')) +
          RotatedAxis())
dev.off()

pdf("Markers_Endo_Vein.pdf", width = 40, height = 30)
print(DotPlot(aSeuratEndo, features = Vein_Markers_Kaluka, cols = c("lightblue", "red3"), 
              col.min = 0, dot.scale = 30) +
        theme(axis.title = element_text(size = 60), axis.text = element_text(size = 50), 
              legend.text = element_text(size = 40), legend.key.size = unit(3, 'cm'), 
              legend.title  = element_text(size = 40), legend.spacing = unit(3, 'cm')) +
        RotatedAxis())
dev.off()

pdf("Markers_Endo_Capillary.pdf", width = 40, height = 30)
print(DotPlot(aSeuratEndo, features = Capillary_Markers_Kaluka, cols = c("lightblue", "red3"), 
              col.min = 0, dot.scale = 30) +
        theme(axis.title = element_text(size = 60), axis.text = element_text(size = 50), 
              legend.text = element_text(size = 40), legend.key.size = unit(3, 'cm'), 
              legend.title  = element_text(size = 40), legend.spacing = unit(3, 'cm')) +
        RotatedAxis())
dev.off()

pdf("Markers_Endo_Lymphatic.pdf", width = 40, height = 30)
print(DotPlot(aSeuratEndo, features = Lymphatic_Markers_Kaluka, cols = c("lightblue", "red3"), 
              col.min = 0, dot.scale = 30) +
        theme(axis.title = element_text(size = 60), axis.text = element_text(size = 50), 
              legend.text = element_text(size = 40), legend.key.size = unit(3, 'cm'), 
              legend.title  = element_text(size = 40), legend.spacing = unit(3, 'cm')) +
        RotatedAxis())
dev.off()

endo_lables_long_dot <- c("0_Capillary", "1_Capillary", "2_Capillary", "3_Capillary", "4_Artery", "5_Capillary", 
                      "6_Vein", "7_Artery", "8_Capillary", "9_Capillary", "10_Capillary", "11_Capillary", "12_RBC", 
                      "13_Capillary", "14_Capillary", "15_Lymphatic")

endo_lables_short_dot <- c("Capillary", "Capillary", "Capillary", "Capillary", "Artery", "Capillary", 
                          "Vein", "Artery", "Capillary", "Capillary", "Capillary", "Capillary", "RBC", 
                          "Capillary", "Capillary", "Lymphatic")

#Subset of markers that can show all the subclusters have been IDed
subsetIDMarkers<- c("Apoe", "Il6st", "Sema3g", "Clu", "Ccl21a", "Mmrn1", "Gpihbp1", "Rgcc", "Kdr", "Hbb-bs", "Hba-a1")
pdf("ClusterMarkers_subset_clusterID.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratEndo, features = subsetIDMarkers, cols = c("lightblue", "red3"), col.min = 0, dot.scale = 30) +
        theme(axis.title = element_text(size = 100), axis.text = element_text(size = 100), 
              legend.text = element_text(size = 100), legend.key.size = unit(3, "cm"), 
              legend.title  = element_text(size = 100), legend.spacing = unit(3, "cm"), 
              axis.text.x = element_text(angle= 90, vjust = 0.45, hjust = 0.2)) +
        labs(x = "Vein Artery Lymph Capillary RBC"))
dev.off()

#Rename the clusters based on the ID name
Idents(aSeuratEndo) <- "seurat_clusters"
newClusterIDs <- endo_lables_short_dot
names(newClusterIDs) <- levels(aSeuratEndo)
aSeuratEndo <- RenameIdents(aSeuratEndo, newClusterIDs)
levels(aSeuratEndo)
pdf("ClusterMarkers_subsetID.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratEndo, features = subsetIDMarkers, cols = c("lightblue", "red3"), col.min = 0, dot.scale = 30) +
        theme(axis.title = element_text(size = 100), axis.text = element_text(size = 100), 
              legend.text = element_text(size = 100), legend.key.size = unit(3, "cm"), 
              legend.title  = element_text(size = 100), legend.spacing = unit(3, "cm"), 
              axis.text.x = element_text(angle= 90, vjust = 0.45)))
dev.off()

#Print UMAP using the cluster IDs
pdf("UMAP_Endo_NewID_dot.pdf", width = 40, height = 30)
print(DimPlot(aSeuratEndo, reduction = "umap", label = TRUE, label.size = 16, pt.size = 0.5, repel = TRUE) +
        #NoLegend() +
        theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.text = element_text(size = 40)))
dev.off()

#dotplot of markers after cluster lableing
aSeuratEndo$ID <- aSeuratEndo@active.ident
aSeuratEndo$protocol_subcluster <- paste0(aSeuratEndo$protocol, '_', aSeuratEndo$ID)
Idents(aSeuratEndo) <- "protocol_subcluster"
levels(aSeuratEndo)
table(Idents(aSeuratEndo))
levels(aSeuratEndo) <- c("Deligated_Vein", "Ligated_Vein", "Mock_Vein", 
                        "Deligated_Lymphatic", "Ligated_Lymphatic", "Mock_Lymphatic", 
                        "Deligated_Capillary", "Ligated_Capillary", "Mock_Capillary", 
                        "Deligated_Artery", "Ligated_Artery", "Mock_Artery", 
                        "Mock_RBC", "Ligated_RBC", "Deligated_RBC")

levels(aSeuratEndo)
#Reorder the clusters to be figure friendly
endoSubsetIdents <- c("Deligated_Vein", "Ligated_Vein", "Mock_Vein", 
                      "Deligated_Lymphatic", "Ligated_Lymphatic", "Mock_Lymphatic", 
                      "Deligated_Capillary", "Ligated_Capillary", "Mock_Capillary", 
                      "Deligated_Artery", "Ligated_Artery", "Mock_Artery")

NolanListSubset <- c("Kitl", "Col14a1", "Egfl7", "Gja1", "Gja4", "Lama4", "Cxcl12")

pdf("./Endo-markers/NolanMarkers_subset_short.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratEndo, features = NolanListSubset, cols = c("lightblue", "red3"), col.min = 0, dot.scale = 30,
              idents = endoSubsetIdents) +
        theme(axis.title = element_text(size = 100), axis.text = element_text(size = 100), 
              legend.text = element_text(size = 100), legend.key.size = unit(3, "cm"), 
              legend.title  = element_text(size = 100), legend.spacing = unit(3, "cm"),
              axis.text.x = element_text(angle=90, vjust = 0.45)))
dev.off()

TFlist_short <- c("Hes1", "Snai1", "Snai2", "Zeb1")
pdf("./Endo-markers/TFMarkers_subset_short.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratEndo, features = TFlist_short, cols = c("lightblue", "red3"), col.min = 0, dot.scale = 30,
              idents = endoSubsetIdents) +
        theme(axis.title = element_text(size = 100), axis.text = element_text(size = 100), 
              legend.text = element_text(size = 100), legend.key.size = unit(3, "cm"), 
              legend.title  = element_text(size = 100), legend.spacing = unit(3, "cm"),
              axis.text.x = element_text(angle=90, vjust = 0.45)))
dev.off()

Meselist_short <- c("Ly6a", "Vim")
pdf("./Endo-markers/MeseMarkers_subset_short.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratEndo, features = Meselist_short, cols = c("lightblue", "red3"), col.min = 0, dot.scale = 30,
              idents = endoSubsetIdents) +
        theme(axis.title = element_text(size = 100), axis.text = element_text(size = 100), 
              legend.text = element_text(size = 100), legend.key.size = unit(3, "cm"), 
              legend.title  = element_text(size = 100), legend.spacing = unit(3, "cm"),
              axis.text.x = element_text(angle=90, vjust = 0.45)))
dev.off()

Endolist_short <- c("Cdh5", "Fgfr1", "Flt4", "Kdr", "Nos3", "Pecam1", "Tal1", "Tie1", "Tek", "Vwf")
pdf("./Endo-markers/EndoMarkers_subset_short.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratEndo, features = Endolist_short, cols = c("lightblue", "red3"), col.min = 0, dot.scale = 30,
              idents = endoSubsetIdents) +
        theme(axis.title = element_text(size = 100), axis.text = element_text(size = 100), 
              legend.text = element_text(size = 80), legend.key.size = unit(2, "cm"), 
              legend.title  = element_text(size = 80), legend.spacing = unit(2, "cm"),
              axis.text.x = element_text(angle=90, vjust = 0.45)))
dev.off()

Infllist_short <- c("Jak2", "Rela", "Stat3")
pdf("./Endo-markers/InflMarkers_subset_short.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratEndo, features = Infllist_short, cols = c("lightblue", "red3"), col.min = 0, dot.scale = 30,
              idents = endoSubsetIdents) +
        theme(axis.title = element_text(size = 100), axis.text = element_text(size = 100), 
              legend.text = element_text(size = 100), legend.key.size = unit(3, "cm"), 
              legend.title  = element_text(size = 100), legend.spacing = unit(3, "cm"),
              axis.text.x = element_text(angle=90, vjust = 0.45)))
dev.off()

#Save the Endothelial subset seurat object
saveRDS(aSeuratEndo, file = "aSeuratEndo.rds")
aSeuratEndo <- readRDS("aSeuratEndo.rds")

#Create a seurat object that containes the original clusters without Endothelial cluster and the Endothelial with subclusters

aSeuratIntegratedNoEndo <- subset(aSeuratIntegrated, 
                                  ID == "Immune" | 
                                  ID == "Fibroblast" | 
                                  ID == "Epithelial" | 
                                  ID == "Pericyte" | 
                                  ID == "Glia")

aSeurat_subset <- merge(aSeuratIntegratedNoEndo, y = c(aSeuratEndo))

DefaultAssay(aSeurat_subset) <- "RNA"
aSeurat_subset <- SCTransform(aSeurat_subset, vst.flavor = "v2", verbose = FALSE)
DefaultAssay(aSeurat_subset) <- "SCT"

saveRDS(aSeurat_subset, file = "aSeuratSubset.rds")
aSeurat_subset <- readRDS("aSeuratSubset.rds")
levels(aSeurat_subset)

Cell_Cycle_Genes_of_Interest <- c("Cenpa", "Cenpe", "Cenpf", "Mki67",
                                  "Hist1h1a", "Hist1h1b", "Top2a", "Hist1h2ae", "Pcna", "Mcm2")
Idents(aSeuratEndo) <- "protocol_subcluster"
levels(aSeuratEndo)
pdf("./Endo-markers/Cell_Cycle_Endo_Subset.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratEndo, features = Cell_Cycle_Genes_of_Interest, cols = c("lightblue", "red3"), col.min = 0, dot.scale = 30,
              idents = endoSubsetIdents) +
        theme(axis.title = element_text(size = 100), axis.text = element_text(size = 100), 
              legend.text = element_text(size = 100), legend.key.size = unit(3, "cm"), 
              legend.title  = element_text(size = 100), legend.spacing = unit(3, "cm"), 
              axis.text.x = element_text(angle=90, vjust = 0.45)))
dev.off()

#Find what stage of cell cycle the cells are in using seurat cellcyclescoring
#Genes from seurat are human gene names, convert to mouse gene names
s.genes <- convert_human_to_mouse(cc.genes$g2m.genes)
g2m.genes <- convert_human_to_mouse(cc.genes$s.genes)

aSeuratEndo <- CellCycleScoring(aSeuratEndo, s.features = s.genes, g2m.features = g2m.genes)
levels(aSeuratEndo)
#Print the number of cells in each cluster by cell cycle phase
for(i in levels(aSeuratEndo)){
  tmp <- subset(aSeuratEndo, idents = i)
  Idents(tmp) <- "Phase"
  print(i)
  print(table(Idents(tmp)))
}

#Cellchat and NATMI of sub-clustering--------------

aSeurat_subset <- readRDS("aSeuratSubset.rds")

seuratList <- c()
seuratList <- c(seuratList, subset(aSeurat_subset, protocol == "Mock"))
seuratList <- c(seuratList, subset(aSeurat_subset, protocol == "Ligated"))
seuratList <- c(seuratList, subset(aSeurat_subset, protocol == "Deligated"))

cellchatList <- c()
for(i in 1:3){
  aCellChat <- cellChatFromSeurat(seuratList[[i]], seuratList[[i]]$protocol[1])
  saveRDS(aCellChat, file = sprintf("CellChat_subset%s.rds", seuratList[[i]]$protocol[1]))
  cellchatList <- c(cellchatList, aCellChat)
}

aCellchatMerge <- mergeCellChat(cellchatList, cell.prefix = TRUE, add.names = c(cellchatList[[1]]@meta$protocol[1], 
                                                                                cellchatList[[2]]@meta$protocol[1], 
                                                                                cellchatList[[3]]@meta$protocol[1]))

aCellchatMerge <- computeNetSimilarityPairwise(aCellchatMerge, type = "structural")
aCellchatMerge <- netEmbedding(aCellchatMerge, type = "structural", umap.method = "uwot")
aCellchatMerge <- netClustering(aCellchatMerge, type = "structural", do.parallel = FALSE)

aCellchatMerge <- computeNetSimilarityPairwise(aCellchatMerge, type = "functional")
aCellchatMerge <- netEmbedding(aCellchatMerge, type = "functional", umap.method = "uwot")
aCellchatMerge <- netClustering(aCellchatMerge, type = "functional", do.parallel = FALSE)

saveRDS(aCellchatMerge, file = "CellChat_subsetMerge.rds")

cellchatSubList <- c()
aCellChat <- readRDS("CellChat_subsetMock.rds")
cellchatSubList <- c(cellchatSubList, aCellChat)
aCellChat <- readRDS("CellChat_subsetLigated.rds")
cellchatSubList <- c(cellchatSubList, aCellChat)
aCellChat <- readRDS("CellChat_subsetDeligated.rds")
cellchatSubList <- c(cellchatSubList, aCellChat)
aCellchatSubMerge <- readRDS("CellChat_subsetMerge.rds")

levels(cellchatSubList[[1]]@idents)
nlevels(cellchatSubList[[1]]@idents)

#Source clusters of Artery, Capillary, Lymph or Vein to all clusters
sourceCluster <- c(1, 2, 7, 10)
for(i in sourceCluster){
  for(j in 1:nlevels(cellchatSubList[[1]]@idents)){
    netVisual_bubble_MLD(aCellchatSubMerge, source = c(i), target = c(j), apath = "netVisualBubble2_subset")
  }
}

#Generate the files for NATMI on seurat with sub-endotheal cluster IDs, NATMI needs to be run on the commandline
seuratList <- c()
seuratList <- c(seuratList, subset(aSeurat_subset, protocol == "Mock"))
seuratList <- c(seuratList, subset(aSeurat_subset, protocol == "Ligated"))
seuratList <- c(seuratList, subset(aSeurat_subset, protocol == "Deligated"))

for(i in seuratList){
  #From NATMI
  write.csv(100 * (exp(as.matrix(GetAssayData(object = i, assay = "SCT", slot = "data"))) - 1), 
            sprintf("data-%s.csv", i@meta.data$protocol[1]), row.names = T)
  write.csv(Idents(object = i), sprintf("metadata-%s.csv", i@meta.data$protocol[1]), row.names = T)
}

  seuratList <- c()
  seuratList <- c(seuratList, subset(aSeurat_subsetTMP, protocol == "Mock"))
  seuratList <- c(seuratList, subset(aSeurat_subsetTMP, protocol == "Ligated"))
  seuratList <- c(seuratList, subset(aSeurat_subsetTMP, protocol == "Deligated"))
  
  for(i in seuratList){
    #From NATMI
    write.csv(100 * (exp(as.matrix(GetAssayData(object = i, assay = "SCT", slot = "data"))) - 1), 
              sprintf("data2-%s.csv", i@meta.data$protocol[1]), row.names = T)
    write.csv(Idents(object = i), sprintf("meta2data-%s.csv", i@meta.data$protocol[1]), row.names = T)
}

#Create a list of all pathways for all protocols
pathway.union <- union(cellchatSubList[[indexList$mock]]@netP$pathways, 
                       cellchatSubList[[indexList$ligated]]@netP$pathways)
pathway.union <- union(pathway.union, 
                       cellchatSubList[[indexList$deligated]]@netP$pathways)

#Create a list of only the top 10 pathways for all protocols
pathway.union.short <- union(head(cellchatSubList[[indexList$mock]]@netP$pathways, 10), 
                       head(cellchatSubList[[indexList$ligated]]@netP$pathways, 10))
pathway.union.short <- union(pathway.union.short, 
                       head(cellchatSubList[[indexList$deligated]]@netP$pathways, 10))
  
pdf("netAnalysis_signalingRole_heatmap_sub_mock.pdf", height = 9)
  print(netAnalysis_signalingRole_heatmap(cellchatSubList[[indexList$mock]], pattern = "all", 
                                        signaling = pathway.union, title = "Mock", 
                                        width = 10, height = 18, font.size = 8, color.heatmap = "OrRd"))
dev.off()
pdf("netAnalysis_signalingRole_heatmap_sub_ligated.pdf", height = 9)
  print(netAnalysis_signalingRole_heatmap(cellchatSubList[[indexList$ligated]], pattern = "all", 
                                        signaling = pathway.union, title = "Ligated", 
                                        width = 10, height = 18, font.size = 8, color.heatmap = "OrRd"))
dev.off()
pdf("netAnalysis_signalingRole_heatmap_sub_deligated.pdf", height = 9)
  print(netAnalysis_signalingRole_heatmap(cellchatSubList[[indexList$deligated]], pattern = "all", 
                                        signaling = pathway.union, title = "Deligated", 
                                        width = 10, height = 18, font.size = 8, color.heatmap = "OrRd"))
dev.off()

pdf("netAnalysis_signalingRole_heatmap_sub_short_mock.pdf", height = 9)
print(netAnalysis_signalingRole_heatmap(cellchatSubList[[indexList$mock]], pattern = "all", 
                                        signaling = pathway.union.short, title = "Mock", 
                                        width = 10, height = 18, font.size = 8, color.heatmap = "OrRd"))
dev.off()
pdf("netAnalysis_signalingRole_heatmap_sub_short_ligated.pdf", height = 9)
print(netAnalysis_signalingRole_heatmap(cellchatSubList[[indexList$ligated]], pattern = "all", 
                                        signaling = pathway.union.short, title = "Ligated", 
                                        width = 10, height = 18, font.size = 8, color.heatmap = "OrRd"))
dev.off()
pdf("netAnalysis_signalingRole_heatmap_sub_short_deligated.pdf", height = 9)
print(netAnalysis_signalingRole_heatmap(cellchatSubList[[indexList$deligated]], pattern = "all", 
                                        signaling = pathway.union.short, title = "Deligated", 
                                        width = 10, height = 18, font.size = 8, color.heatmap = "OrRd"))
dev.off()
  


#Reading and subsetting NATMI data---------------------
#NATMI files are rather large, helper functions to search for individual ligand/target pairs from Endothelial subsets

NATMIFolder <- "JDRSubsetMockvLigated"
fileName <- sprintf("%s/%s/Delta_edges_lrc2p/DOWN-regulated_mean.csv", projectPath, NATMIFolder)
NATMImockvligated <- read.csv(fileName)

fileName <- sprintf("%s/%s/Delta_edges_lrc2p/UP-regulated_mean.csv", projectPath, NATMIFolder)
NATMIligatedvmock <- read.csv(fileName)

NATMIFolder <- "JDRSubsetMockvDeligated"
fileName <- sprintf("%s/%s/Delta_edges_lrc2p/DOWN-regulated_mean.csv", projectPath, NATMIFolder)
NATMImockvdeligated <- read.csv(fileName)

fileName <- sprintf("%s/%s/Delta_edges_lrc2p/UP-regulated_mean.csv", projectPath, NATMIFolder)
NATMIdeligatedvmock <- read.csv(fileName)

NATMIFolder <- "JDRSubsetLigatedvDeligated"
fileName <- sprintf("%s/%s/Delta_edges_lrc2p/DOWN-regulated_mean.csv", projectPath, NATMIFolder)
NATMIligatedvdeligated <- read.csv(fileName)

fileName <- sprintf("%s/%s/Delta_edges_lrc2p/UP-regulated_mean.csv", projectPath, NATMIFolder)
NATMIligatedvdeligated <- read.csv(fileName)

#Create several tables of Ligand-Receptor sets from Endothelial subsets
ligand <- "Sema7a"
target <- "Fibroblast"
NATMITable <- NATMIligatedvmock
aNATMISubset <- subset(NATMITable, Sending.cluster == "Artery" & Ligand.symbol == ligand & Target.cluster == target)[,c(1, 2, 3, 4, 39)]
aNATMIArtery <- aNATMISubset[order(aNATMISubset$Target.cluster, aNATMISubset$Ligand.symbol),]
aNATMISubset <- subset(NATMITable, Sending.cluster == "Capillary" & Ligand.symbol == ligand & Target.cluster == target)[,c(1, 2, 3, 4, 39)]
aNATMICapillary <- aNATMISubset[order(aNATMISubset$Target.cluster, aNATMISubset$Ligand.symbol),]
aNATMISubset <- subset(NATMITable, Sending.cluster == "Lymphatic" & Ligand.symbol == ligand & Target.cluster == target)[,c(1, 2, 3, 4, 39)]
aNATMILymph <- aNATMISubset[order(aNATMISubset$Target.cluster, aNATMISubset$Ligand.symbol),]
aNATMISubset <- subset(NATMITable, Sending.cluster == "Vein" & Ligand.symbol == ligand & Target.cluster == target)[,c(1, 2, 3, 4, 39)]
aNATMIVein <- aNATMISubset[order(aNATMISubset$Target.cluster, aNATMISubset$Ligand.symbol),]



#Homeostatic-------------------------
#Data has already been processed, loading, cluster IDs and figure generation
#Multi organ data from Kaluca et al. 2020 (PMID: 32059779)

homeostaticFile <- "c:/homeostatic_Merged.rds"
aSeuratHomeostatic <- readRDS(homeostaticFile)

#Print UMAP of homeostatic clusters
Idents(aSeuratHomeostatic) <- "seurat_clusters"
pdf("UMAP_homeostatic.pdf", width = 40, height = 30)
print(DimPlot(aSeuratHomeostatic, reduction = "umap", label = TRUE, label.size = 40, pt.size = 0.5, repel = TRUE) +
        NoLegend() +
        theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.text = element_text(size = 40)))
dev.off()

#set ID for clusters
Idents(aSeuratHomeostatic) <- "seurat_clusters"
levels(aSeuratHomeostatic)
newClusterIDsLong <- c("0_Epithelial", "1_Endothelial", "2_Endothelial", "3_Epithelial", "4_Immune", "5_Immune", 
                       "6_Immune", "7_Endothelial", "8_Epithelial", "9_Endothelial", "10_Fibroblast", 
                       "11_Immune", "12_Immune", "13_Immune", "14_Endothelial", "15_RBC", 
                       "16_Immune", "17_Pericyte", "18_Immune", "19_Immune", "20_Endothelial", 
                       "21_Epithelial", "22_Stroma", "23_Immune", "24_Immune", "25_Endothelial",
                       "26_Fibroblast")
newClusterIDsShort <- c("Epithelial", "Endothelial", "Endothelial", "Epithelial", "Immune", "Immune", 
                        "Immune", "Endothelial", "Epithelial", "Endothelial", "Fibroblast", 
                        "Immune", "Immune", "Immune", "Endothelial", "RBC", 
                        "Immune", "Pericyte", "Immune", "Immune", "Endothelial", 
                        "Epithelial", "Stroma", "Immune", "Immune", "Endothelial",
                        "Fibroblast")

Idents(aSeuratHomeostatic) <- "seurat_clusters"
newClusterIDs <- newClusterIDsShort
names(newClusterIDs) <- levels(aSeuratHomeostatic)
aSeuratHomeostatic <- RenameIdents(aSeuratHomeostatic, newClusterIDs)
levels(aSeuratHomeostatic)

#Print UMAP using cluster names
pdf("UMAP_homeostatic_NewID.pdf", width = 40, height = 30)
print(DimPlot(aSeuratHomeostatic, reduction = "umap", label = TRUE, label.size = 40, pt.size = 0.5, repel = TRUE) +
        NoLegend() +
        theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.text = element_text(size = 40)))
dev.off()

#Load seurat object of Kaluka organ data
kaluckaFile <- "c:/Kalucka_SMG_Integrated.rds"
aSeuratkalucka <- readRDS(kaluckaFile)
levels(aSeuratkalucka)

#Rename cluster IDs to be more human friendly
newClusterIDsShort <- c("Salivary Gland", "Brain", "Heart", "Colon", "Small Intestine", "Kidney", "Liver",
                        "Lung", "EDL Muscle", "Soleus Muscle", "Spleen", "Testis")

newClusterIDs <- newClusterIDsShort
names(newClusterIDs) <- levels(aSeuratkalucka)
aSeuratkalucka <- RenameIdents(aSeuratkalucka, newClusterIDs)
levels(aSeuratkalucka)

#Print UMAP of Kalucka mutli organ with new organ names
pdf("UMAP_kalucka.pdf", width = 40, height = 30)
print(DimPlot(aSeuratkalucka, reduction = "umap", label = FALSE, label.size = 40, pt.size = 0.5, 
              repel = TRUE, order = "Salivary Gland") +
        #NoLegend() +
        theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40), legend.text = element_text(size = 120)))
dev.off()

#CLDN5 MCP1----------------------------
#Look for tight junction genes in all data and Endothelial subsets

aSeuratIntegrated <- readRDS(sprintf("%s/JDR2023Seurat-Meta.rds", dataPath))
aSeuratEndo <- readRDS("aSeuratEndo.rds")

Idents(aSeuratIntegrated) <- aSeuratIntegrated$ID
Idents(aSeuratEndo) <- aSeuratEndo$ID

levels(aSeuratIntegrated)
levels(aSeuratEndo)

#Dotplot of Cldn5 in Endothelial subset
pdf("DotPlot_ENdosubset_CLDN5.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratEndo, features = "Cldn5", cols = c("lightblue", "red3"), 
              idents = endoSubsetIdents, col.min = 0, dot.scale = 60) +
        theme(axis.title = element_text(size = 100), axis.text = element_text(size = 100), 
              legend.text = element_text(size = 100), legend.key.size = unit(3, "cm"), 
              legend.title  = element_text(size = 100), legend.spacing = unit(3, "cm"), 
              axis.text.x = element_text(angle=90, vjust = 0.45)))
dev.off()

#Dotplot of Mcp1 (Ccl2) in all clusters
pdf("DotPlot_ENdosubset_MCP1.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratIntegrated, features = "Ccl2", cols = c("lightblue", "red3"), 
              col.min = 0, dot.scale = 50) +
        theme(axis.title = element_text(size = 100), axis.text = element_text(size = 80), 
              legend.text = element_text(size = 100), legend.key.size = unit(3, "cm"), 
              legend.title  = element_text(size = 100), legend.spacing = unit(3, "cm"), 
              axis.text.x = element_text(angle=90, vjust = 0.45)))
dev.off()

tight_JunctionList <- c("Cldn2", "Cldn3", "Cldn4", "Cldn5", "Cldn6", "Cldn7", "Cldn8", "Cldn9" ,"Cldn10",
                        "Cldn11", "Cldn12", "Cldn13", "Cldn14", "Cldn15", "Cldn16", "Cldn17", "Cldn18",
                        "Ocln", "Marveld3", "F11r", "Jam2", "Jam3")

#Dotplot of tight juncion genes in all clusters (F11r is Jam1)
pdf("DotPlot_ENdosubset_tightJunctions.pdf", width = 40, height = 30, bg = "white")
print(DotPlot(aSeuratEndo, features = tight_JunctionList, cols = c("lightblue", "red3"), 
              idents = endoSubsetIdents, col.min = 0, dot.scale = 40) +
        theme(axis.title = element_text(size = 90), axis.text = element_text(size = 90), 
              legend.text = element_text(size = 80), legend.key.size = unit(2, "cm"), 
              legend.title  = element_text(size = 80), legend.spacing = unit(2, "cm"), 
              axis.text.x = element_text(angle=90, vjust = 0.45)))
dev.off()



#SessionInfo---------------------
# R version 4.3.1 (2023-06-16 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19044)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    
# 
# time zone: America/New_York
# tzcode source: internal
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] DoubletFinder_2.0.4         CellChat_1.6.1              igraph_1.5.0                remotes_2.4.2.1            
# [5] rhdf5_2.44.0                hdf5r_1.3.8                 clustree_0.5.0              ggraph_2.1.0               
# [9] biomaRt_2.56.1              singleCellTK_2.10.0         DelayedArray_0.26.6         S4Arrays_1.0.4             
# [13] scDblFinder_1.14.0          scCustomize_1.1.3           stringr_1.5.0               glmGamPoi_1.11.7           
# [17] sctransform_0.3.5           dplyr_1.1.2                 patchwork_1.1.3             scater_1.28.0              
# [21] scuttle_1.10.1              cowplot_1.1.1               ggplot2_3.4.3               SeuratObject_4.1.3         
# [25] Seurat_4.3.0.1              celda_1.16.1                Matrix_1.6-0                SingleCellExperiment_1.22.0
# [29] SummarizedExperiment_1.30.2 Biobase_2.60.0              GenomicRanges_1.52.0        GenomeInfoDb_1.36.1        
# [33] IRanges_2.34.1              S4Vectors_0.38.1            BiocGenerics_0.46.0         MatrixGenerics_1.12.2      
# [37] matrixStats_1.0.0          
# 
# loaded via a namespace (and not attached):
#   [1] R.methodsS3_1.8.2          progress_1.2.2             urlchecker_1.0.1           goftest_1.2-3             
# [5] Biostrings_2.68.1          HDF5Array_1.28.1           vctrs_0.6.3                spatstat.random_3.1-5     
# [9] digest_0.6.33              png_0.1-8                  shape_1.4.6                registry_0.5-1            
# [13] ggrepel_0.9.3              deldir_1.0-9               parallelly_1.36.0          combinat_0.0-8            
# [17] magick_2.7.4               MASS_7.3-60                reshape2_1.4.4             httpuv_1.6.11             
# [21] foreach_1.5.2              withr_2.5.0                ggrastr_1.0.2              ggpubr_0.6.0              
# [25] ellipsis_0.3.2             survival_3.5-5             memoise_2.0.1              ggbeeswarm_0.7.2          
# [29] janitor_2.2.0              profvis_0.3.8              systemfonts_1.0.4          zoo_1.8-12                
# [33] GlobalOptions_0.1.2        pbapply_1.7-2              R.oo_1.25.0                prettyunits_1.1.1         
# [37] rematch2_2.1.2             KEGGREST_1.40.0            promises_1.2.0.1           httr_1.4.6                
# [41] rstatix_0.7.2              restfulr_0.0.15            globals_0.16.2             fitdistrplus_1.1-11       
# [45] rhdf5filters_1.12.1        ps_1.7.5                   rstudioapi_0.15.0          miniUI_0.1.1.1            
# [49] generics_0.1.3             ggalluvial_0.12.5          processx_3.8.2             curl_5.0.1                
# [53] zlibbioc_1.46.0            ScaledMatrix_1.8.1         polyclip_1.10-4            GenomeInfoDbData_1.2.10   
# [57] RcppEigen_0.3.3.9.3        xtable_1.8-4               doParallel_1.0.17          BiocFileCache_2.8.0       
# [61] hms_1.1.3                  irlba_2.3.5.1              colorspace_2.1-0           filelock_1.0.2            
# [65] ggnetwork_0.5.12           ROCR_1.0-11                reticulate_1.30            spatstat.data_3.0-1       
# [69] magrittr_2.0.3             lmtest_0.9-40              snakecase_0.11.0           later_1.3.1               
# [73] viridis_0.6.4              lattice_0.21-8             spatstat.geom_3.2-2        NMF_0.26                  
# [77] future.apply_1.11.0        scattermore_1.2            XML_3.99-0.14              GSVAdata_1.36.0           
# [81] assertive.numbers_0.0-2    RcppAnnoy_0.0.21           pillar_1.9.0               nlme_3.1-162              
# [85] sna_2.7-1                  iterators_1.0.14           gridBase_0.4-7             compiler_4.3.1            
# [89] beachmat_2.16.0            RSpectra_0.16-1            stringi_1.7.12             eds_1.2.0                 
# [93] assertive.properties_0.0-5 tensor_1.5                 lubridate_1.9.2            devtools_2.4.5            
# [97] GenomicAlignments_1.36.0   MCMCprecision_0.4.0        plyr_1.8.8                 crayon_1.5.2              
# [101] abind_1.4-7                BiocIO_1.10.0              gridGraphics_0.5-1         locfit_1.5-9.8            
# [105] sp_2.0-0                   graphlayouts_1.0.0         bit_4.0.5                  codetools_0.2-19          
# [109] BiocSingular_1.16.0        paletteer_1.5.0            GetoptLong_1.0.5           plotly_4.10.2             
# [113] mime_0.12                  splines_4.3.1              circlize_0.4.15            Rcpp_1.0.11               
# [117] dbplyr_2.3.3               sparseMatrixStats_1.12.2   blob_1.2.4                 utf8_1.2.3                
# [121] clue_0.3-64                WriteXLS_6.4.0             fs_1.6.2                   listenv_0.9.0             
# [125] DelayedMatrixStats_1.22.1  pkgbuild_1.4.2             ggsignif_0.6.4             tibble_3.2.1              
# [129] assertive.base_0.0-9       callr_3.7.3                statmod_1.5.0              svglite_2.1.1             
# [133] network_1.18.1             tweenr_2.0.2               pkgconfig_2.0.3            tools_4.3.1               
# [137] cachem_1.0.8               RSQLite_2.3.1              viridisLite_0.4.2          DBI_1.1.3                 
# [141] fastmap_1.1.1              scales_1.2.1               grid_4.3.1                 usethis_2.2.2             
# [145] ica_1.0-3                  Rsamtools_2.16.0           broom_1.0.5                coda_0.19-4               
# [149] FNN_1.1.3.2                BiocManager_1.30.22        ggprism_1.0.4              carData_3.0-5             
# [153] RANN_2.6.1                 farver_2.1.1               tidygraph_1.2.3            yaml_2.3.7                
# [157] rtracklayer_1.60.0         cli_3.6.1                  assertive.types_0.0-3      purrr_1.0.1               
# [161] leiden_0.4.3               lifecycle_1.0.3            uwot_0.1.16                backports_1.4.1           
# [165] bluster_1.10.0             sessioninfo_1.2.2          assertive.files_0.0-2      DropletUtils_1.20.0       
# [169] BiocParallel_1.34.2        timechange_0.2.0           gtable_0.3.4               rjson_0.2.21              
# [173] ggridges_0.5.4             progressr_0.13.0           parallel_4.3.1             limma_3.56.2              
# [177] jsonlite_1.8.7             edgeR_3.42.4               bitops_1.0-7               bit64_4.0.5               
# [181] xgboost_1.7.5.1            Rtsne_0.16                 spatstat.utils_3.0-3       BiocNeighbors_1.18.0      
# [185] metapod_1.8.0              dqrng_0.3.0                enrichR_3.2                R.utils_2.12.2            
# [189] lazyeval_0.2.2             shiny_1.7.4.1              htmltools_0.5.5            rappdirs_0.3.3            
# [193] glue_1.6.2                 XVector_0.40.0             RCurl_1.98-1.12            scran_1.28.1              
# [197] gridExtra_2.3              R6_2.5.1                   tidyr_1.3.0                forcats_1.0.0             
# [201] rngtools_1.5.2             cluster_2.1.4              pkgload_1.3.2.1            Rhdf5lib_1.22.0           
# [205] statnet.common_4.9.0       tidyselect_1.2.0           vipor_0.4.5                ggforce_0.4.1             
# [209] xml2_1.3.5                 car_3.1-3                  AnnotationDbi_1.62.2       future_1.33.0             
# [213] rsvd_1.0.5                 munsell_0.5.0              KernSmooth_2.23-21         multipanelfigure_2.1.2    
# [217] data.table_1.14.8          ComplexHeatmap_2.16.0      htmlwidgets_1.6.2          RColorBrewer_1.1-3        
# [221] rlang_1.1.1                spatstat.sparse_3.0-2      spatstat.explore_3.2-1     fansi_1.0.4               
# [225] beeswarm_0.4.0  
