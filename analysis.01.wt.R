###################################################################################################

                                 ##### Processing of WT Data #####

###################################################################################################

{
  library(Seurat)
  library(Signac)
  library(GenomeInfoDb)
  library(AnnotationHub)
  library(dplyr)
  library(ggplot2)
}

{
  INPUT.FILE1 = '~/10x/data/wt/filtered_feature_bc_matrix.h5'
  INPUT.FILE2 = '~/10x/data/wt/atac_fragments.tsv.gz'
  
  OUTPUT.FILE1 = 'ANL001a.wt.data.processed.rds'
}

data.input = Read10X_h5(INPUT.FILE1)

rna.counts = data.input$`Gene Expression`
atac.counts = data.input$Peaks

sobj = CreateSeuratObject(counts = rna.counts)
sobj[['percent.mt']] = PercentageFeatureSet(sobj, pattern = '^mt-')

grange.counts = StringToGRanges(rownames(atac.counts), sep = c(':', '-'))
if ( !all(standardChromosomes(grange.counts) %in% as.character(1:25)) ) {
  grange.use = seqnames(grange.counts) %in% as.character(1:25)
  atac.counts = atac.counts[as.vector(grange.use), ]
}

EDB = query(AnnotationHub(), pattern = c('Danio rerio', 'EnsDb', 103))[[1]]
annotations = GetGRangesFromEnsDb(ensdb = EDB)
genome(annotations) = 'danRer11'

sobj[['ATAC']] = CreateChromatinAssay(
  counts = atac.counts,
  sep = c(':', '-'),
  genome = 'danRer11',
  fragments = INPUT.FILE2,
  min.cells = 10,
  annotation = annotations
)

{ ### QC and filtering ###
  VlnPlot(sobj, features = c('nCount_RNA', 'nCount_ATAC', 'percent.mt'), ncol = 3,
          group.by = 'orig.ident', log = F, pt.size = 0) + NoLegend()
  VlnPlot(sobj, features = c('nCount_RNA', 'nCount_ATAC', 'percent.mt'), ncol = 3,
          group.by = 'orig.ident', log = T, pt.size = 0) + NoLegend()
  
  sobj.filt = subset(
    x = sobj, 
    subset = nCount_RNA < 3e3 & nCount_RNA > 50 & 
      nCount_ATAC < 5e4 & nCount_ATAC > 500 & 
      percent.mt < 15
  )
  
  VlnPlot(sobj.filt, features = c('nCount_RNA', 'nCount_ATAC', 'percent.mt'), ncol = 3,
          group.by = 'orig.ident', log = T, pt.size = 0) + NoLegend()
}

{
  DefaultAssay(sobj.filt) = 'ATAC'
  
  sobj.filt = NucleosomeSignal(sobj.filt)
  sobj.filt = TSSEnrichment(sobj.filt)

  VlnPlot(
    object = sobj.filt, 
    features = c('TSS.enrichment', 'nucleosome_signal'), 
    group.by = 'orig.ident', 
    ncol = 2, 
    pt.size = 0
  )
}

##### Dimensional reduction by UMAP #####
### RNA analysis ###
DefaultAssay(sobj.filt) = 'RNA'
sobj.filt = SCTransform(sobj.filt, verbose = F) %>% RunPCA() %>% 
  RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

### ATAC analysis ###
DefaultAssay(sobj.filt) = 'ATAC'
sobj.filt = RunTFIDF(sobj.filt)
sobj.filt = FindTopFeatures(sobj.filt, min.cutoff = 'q0')
sobj.filt = RunSVD(sobj.filt)
sobj.filt = RunUMAP(sobj.filt, reduction = 'lsi', dims = 2:50, 
                    reduction.name = 'umap.atac', reduction.key = 'atacUMAP_')

### WNN graph ###
sobj.filt = FindMultiModalNeighbors(sobj.filt, reduction.list = list('pca', 'lsi'), dims.list = list(1:50, 2:50))
sobj.filt = RunUMAP(sobj.filt, nn.name = 'weighted.nn', reduction.name = 'wnn.umap', reduction.key = 'wnnUMAP_')

sobj.filt = FindClusters(sobj.filt, graph.name = 'wsnn', algorithm = 3, verbose = F)

VlnPlot(
  object = sobj.filt,
  features = c('nCount_RNA', 'percent.mt', 'nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
  group.by = 'seurat_clusters', 
  ncol = 5,
  pt.size = 0
)

{
  p1 = DimPlot(sobj.filt, reduction = 'umap.rna', label = T, label.size = 5, repel = T) + ggtitle('RNA')
  p2 = DimPlot(sobj.filt, reduction = 'umap.atac', label = T, label.size = 5, repel = T) + ggtitle('ATAC')
  p3 = DimPlot(sobj.filt, reduction = 'wnn.umap', label = T, label.size = 5, repel = T) + ggtitle('WNN')
  p1 + p2 + p3 & theme(plot.title = element_text(hjust = 0.5)) #& NoLegend()

  DimPlot(sobj.filt, reduction = 'wnn.umap', group.by = 'seurat_clusters', label = F) + NoLegend() + ggtitle('')
}

saveRDS(sobj.filt, file = OUTPUT.FILE)

