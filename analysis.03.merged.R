###################################################################################################

                                    ##### Merged Analyses #####

###################################################################################################

### !!! Before running the code below, 'analysis.01.wt.R' and 'analysis.01.mut.R' should be run first. !!! ###
{
  INPUT.FILENAME1 = 'ANL001a.wt.data.processed.rds' # from 'analysis.01.wt.R'
  INPUT.FILENAME2 = 'ANL001a.mut.data.processed.rds' # from ''analysis.01.mut.R'
}

{
  library(Seurat)
  library(Signac)
  library(BSgenome.Drerio.UCSC.danRer11)
  library(GenomeInfoDb)
  library(GenomicRanges)
  library(JASPAR2020)
  library(TFBSTools)
  library(motifmatchr)
  library(chromVAR)
  library(BiocParallel)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(patchwork)
  library(corrplot)
}


###################################################################################################

##### Basic comparisons of the WT and Mut data #####

###################################################################################################

(sobj.wt = readRDS(INPUT.FILENAME1))
(sobj.mut = readRDS(INPUT.FILENAME2))

### RNA ###
{
  rbind(WT = summary(sobj.wt@meta.data$nCount_RNA), 
        MUT = summary(sobj.mut@meta.data$nCount_RNA))

  t.test(sobj.wt@meta.data$nCount_RNA, sobj.mut@meta.data$nCount_RNA, alternative = 'l')$p.value
  wilcox.test(sobj.wt@meta.data$nCount_RNA, sobj.mut@meta.data$nCount_RNA, alternative = 'l')$p.value
  ks.test(sobj.wt@meta.data$nCount_RNA, sobj.mut@meta.data$nCount_RNA, alternative = 'g')$p.value

  par(omi = c(0, .2, 0, 0), cex = 1.5)
  hist(sobj.wt@meta.data$nCount_RNA, br = 20, freq = F, col = 'white', border = 'black', 
       xlab = 'Read counts', main = 'RNA')
  hist(sobj.mut@meta.data$nCount_RNA, br = 40, freq = F, col = 'white', border = 'red', add = T)
  
  par(omi = c(0, .2, 0, 0), cex = 1.5)
  plot(density(sobj.wt@meta.data$nCount_RNA), lwd = 2, col = 'black', ylim = c(0, 0.0023), 
       xlab = 'Read counts', main = 'RNA')
  lines(density(sobj.mut@meta.data$nCount_RNA), lwd = 2, col = 'red')
  legend('topright', bty = 'n', col = c('black', 'red'), lty = 1, lwd = 2, 
         legend = c('WT', 'MUT'))
  
  par(omi = c(0, .2, 0, 0), cex = 1.5)
  plot(sample(sobj.wt@meta.data$nCount_RNA, 6000), 
       sample(sobj.mut@meta.data$nCount_RNA, 6000), 
       xlim = c(50, 5000), ylim = c(50, 5000), log = 'xy', cex = .9, 
       xlab = 'Read counts, WT (log scale)', ylab = 'Read counts, Mut (log scale)', 
       main = 'RNA (6,000 random cells)')
  abline(0, 1, col = 'grey', lty = 5)
}

### ATAC ###
{
  rbind(WT = summary(sobj.wt@meta.data$nCount_ATAC), 
        MUT = summary(sobj.mut@meta.data$nCount_ATAC))

  t.test(sobj.wt@meta.data$nCount_ATAC, sobj.mut@meta.data$nCount_ATAC, alternative = 'g')$p.value
  wilcox.test(sobj.wt@meta.data$nCount_ATAC, sobj.mut@meta.data$nCount_ATAC, alternative = 'g')$p.value
  ks.test(sobj.wt@meta.data$nCount_ATAC, sobj.mut@meta.data$nCount_ATAC, alternative = 'l')$p.value

  par(omi = c(0, .2, 0, 0), cex = 1.5)
  hist(sobj.wt@meta.data$nCount_ATAC, br = 20, freq = F, col = 'white', border = 'black', 
       xlab = 'Fragment counts', main = 'ATAC')
  hist(sobj.mut@meta.data$nCount_ATAC, br = 40, freq = F, col = 'white', border = 'red', add = T)
  
  par(omi = c(0, .2, 0, 0), cex = 1.5)
  plot(density(sobj.wt@meta.data$nCount_ATAC), lwd = 2, col = 'black', ylim = c(0, 1.1e-4), 
       xlab = 'Fragment counts', main = 'ATAC')
  lines(density(sobj.mut@meta.data$nCount_ATAC), lwd = 2, col = 'red')
  legend('topright', bty = 'n', col = c('black', 'red'), lty = 1, lwd = 2, 
         legend = c('WT', 'MUT'))
  
  par(omi = c(0, .2, 0, 0), cex = 1.5)
  plot(sample(sobj.wt@meta.data$nCount_ATAC, 6000), 
       sample(sobj.mut@meta.data$nCount_ATAC, 6000), 
       xlim = c(500, 80000), ylim = c(500, 80000), log = 'xy', cex = .9, 
       xlab = 'Fragment counts, WT (log scale)', ylab = 'Fragment counts, Mut (log scale)', 
       main = 'ATAC (6,000 random cells)')
  abline(0, 1, col = 'grey', lty = 5)
}

### Nucleosome signal ###
{
  rbind(WT = summary(sobj.wt@meta.data$nucleosome_signal), 
        MUT = summary(sobj.mut@meta.data$nucleosome_signal))

  t.test(sobj.wt@meta.data$nucleosome_signal, sobj.mut@meta.data$nucleosome_signal, 
         alternative = 'l')$p.value
  wilcox.test(sobj.wt@meta.data$nucleosome_signal, sobj.mut@meta.data$nucleosome_signal, 
              alternative = 'l')$p.value
  ks.test(sobj.wt@meta.data$nucleosome_signal, sobj.mut@meta.data$nucleosome_signal, 
          alternative = 'g')$p.value

  par(omi = c(0, .2, 0, 0), cex = 1.5)
  hist(sobj.wt@meta.data$nucleosome_signal, br = 20, freq = F, col = 'white', border = 'black', 
       xlab = 'Signal values', main = 'Nucleosome signal')
  hist(sobj.mut@meta.data$nucleosome_signal, br = 40, freq = F, col = 'white', border = 'red', add = T)
  
  par(omi = c(0, .2, 0, 0), cex = 1.5)
  plot(density(sobj.wt@meta.data$nucleosome_signal), lwd = 2, col = 'black', ylim = c(0, 5.2), 
       xlab = 'Signal values', main = 'Nucleosome signal')
  lines(density(sobj.mut@meta.data$nucleosome_signal), lwd = 2, col = 'red')
  legend('topright', bty = 'n', col = c('black', 'red'), lty = 1, lwd = 2, 
         legend = c('WT', 'MUT'))
  
  par(omi = c(0, .2, 0, 0), cex = 1.5)
  plot(sample(sobj.wt@meta.data$nucleosome_signal, 6000), 
       sample(sobj.mut@meta.data$nucleosome_signal, 6000), 
       xlim = c(.25, 1.65), ylim = c(.25, 1.65), cex = .9, 
       xlab = 'Signal values, WT', ylab = 'Signal values, Mut', 
       main = 'Nucleosome signal (6,000 random cells)')
  abline(0, 1, col = 'grey', lty = 5)
}

### TSS enrichment ###
{
  rbind(WT = summary(sobj.wt@meta.data$TSS.enrichment), 
        MUT = summary(sobj.mut@meta.data$TSS.enrichment))

  t.test(sobj.wt@meta.data$TSS.enrichment, sobj.mut@meta.data$TSS.enrichment, 
         alternative = 'l')$p.value
  wilcox.test(sobj.wt@meta.data$TSS.enrichment, sobj.mut@meta.data$TSS.enrichment, 
              alternative = 'l')$p.value
  ks.test(sobj.wt@meta.data$TSS.enrichment, sobj.mut@meta.data$TSS.enrichment, 
          alternative = 'g')$p.value

  par(omi = c(0, .2, 0, 0), cex = 1.5)
  hist(sobj.wt@meta.data$TSS.enrichment, br = 20, freq = F, col = 'white', border = 'black', 
       xlab = 'Enrichment values', main = 'TSS enrichment')
  hist(sobj.mut@meta.data$TSS.enrichment, br = 40, freq = F, col = 'white', border = 'red', add = T)
  
  par(omi = c(0, .2, 0, 0), cex = 1.5)
  plot(density(sobj.wt@meta.data$TSS.enrichment), lwd = 2, col = 'black', ylim = c(0, 0.8), 
       xlab = 'Enrichment values', main = 'TSS enrichment')
  lines(density(sobj.mut@meta.data$TSS.enrichment), lwd = 2, col = 'red')
  legend('topright', bty = 'n', col = c('black', 'red'), lty = 1, lwd = 2, 
         legend = c('WT', 'MUT'))
  
  par(omi = c(0, .2, 0, 0), cex = 1.5)
  plot(sample(sobj.wt@meta.data$TSS.enrichment, 6000), 
       sample(sobj.mut@meta.data$TSS.enrichment, 6000), 
       xlim = c(1.7, 28), ylim = c(1.7, 28), log = 'xy', cex = .9, 
       xlab = 'Enrichment values, WT (log scale)', ylab = 'Enrichment values, Mut (log scale)', 
       main = 'TSS enrichment (6,000 random cells)')
  abline(0, 1, col = 'grey', lty = 5)
}



###################################################################################################

##### Create a merged dataset of WT and Mut #####

###################################################################################################

{
  #(sobj.wt = readRDS(INPUT.FILENAME1))
  #(sobj.mut = readRDS(INPUT.FILENAME2))

  DefaultAssay(sobj.wt) = "RNA"
  DefaultAssay(sobj.mut) = "RNA"
  
  system.time({
    sobj.filt = merge(sobj.wt, y = sobj.mut, add.cell.ids = c("WT", "MUT"), project = "Val")
  })
}

sobj.filt$treatment = sapply(strsplit(rownames(sobj.filt@meta.data), '_'), function(x) x[1])

##### Dimensional reduction by UMAP #####
### RNA analysis ###
system.time({
  DefaultAssay(sobj.filt) = "RNA"
  sobj.filt = SCTransform(sobj.filt, verbose = F) %>% RunPCA() %>% 
    RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
})

### ATAC analysis ###
system.time({
  DefaultAssay(sobj.filt) = "ATAC"
  sobj.filt = RunTFIDF(sobj.filt)
  sobj.filt = FindTopFeatures(sobj.filt, min.cutoff = 'q0')
  sobj.filt = RunSVD(sobj.filt)
  sobj.filt = RunUMAP(sobj.filt, reduction = 'lsi', dims = 2:50, 
                      reduction.name = "umap.atac", reduction.key = "atacUMAP_")
})

### WNN graph ###
system.time({
  sobj.filt = FindMultiModalNeighbors(sobj.filt, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
  sobj.filt = RunUMAP(sobj.filt, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
})

### Clustering ###
sobj.filt1 = FindClusters(sobj.filt, graph.name = "wsnn", algorithm = 3, verbose = F)

VlnPlot(
  object = sobj.filt,
  features = c("nCount_RNA", "percent.mt", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  group.by = 'treatment', 
  ncol = 5,
  pt.size = 0
)

VlnPlot(
  object = sobj.filt,
  features = c("nCount_RNA", "percent.mt", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  group.by = 'seurat_clusters', 
  ncol = 5,
  pt.size = 0
)

VlnPlot(
  object = sobj.filt,
  features = c("nCount_RNA", "percent.mt", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  group.by = 'wsnn_res.0.8', 
  ncol = 5,
  pt.size = 0
)

{
  p1 = DimPlot(sobj.filt, reduction = "wnn.umap", group.by = 'treatment', label = F, label.size = 5, repel = T) + ggtitle("WNN")
  p2 = DimPlot(sobj.filt, reduction = "umap.rna", group.by = 'treatment', label = F, label.size = 5, repel = T) + ggtitle("RNA")
  p3 = DimPlot(sobj.filt, reduction = "umap.atac", group.by = 'treatment', label = F, label.size = 5, repel = T) + ggtitle("ATAC")
  p1 + p2 + p3 & theme(plot.title = element_text(hjust = 0.5))
  ggsave('ANL001b.merged.umap.by.genotype.png', width = 24, height = 8)
}

{
  p1 = DimPlot(sobj.filt, reduction = "wnn.umap", label = T, label.size = 5, repel = T) + ggtitle("WNN")
  p2 = DimPlot(sobj.filt, reduction = "umap.rna", label = T, label.size = 5, repel = T) + ggtitle("RNA")
  p3 = DimPlot(sobj.filt, reduction = "umap.atac", label = T, label.size = 5, repel = T) + ggtitle("ATAC")
  p1 + p2 + p3 & theme(plot.title = element_text(hjust = 0.5)) & NoLegend()
  ggsave('ANL001b.merged.umap.by.cluster.png', width = 24, height = 8)
}

{
  DimPlot(sobj.filt, reduction = "wnn.umap", group.by = 'seurat_clusters', label = F) + NoLegend() + ggtitle('')
  ggsave('ANL001b.merged.wnn.umap.bycluster.v1.png', width = 8, height = 8)
}

{
  (clst.cell.tot = table(sobj.filt@meta.data$seurat_clusters))
  (genotype.cell.tot = table(sobj.filt@meta.data$treatment))

  clst.cell.each.tab = NULL
  for ( i in 1:length(clst.cell.tot) ) {
    tmp1 = sobj.filt@meta.data$treatment[sobj.filt@meta.data$seurat_clusters == names(clst.cell.tot)[i]]
    tmp2 = table(tmp1)
    tmp3 = c(names(clst.cell.tot)[i], clst.cell.tot[i], 
             tmp2, 
             round(tmp2 / as.numeric(clst.cell.tot[i]), 3), 
             round(tmp2 / genotype.cell.tot, 3))
    
    clst.cell.each.tab = rbind(clst.cell.each.tab, tmp3)
  }
  
  colnames(clst.cell.each.tab) = c('Cluster', 'Total', 'MUT', 'WT', 'MUT_lc', 'WT_lc', 'MUT_gb', 'WT_gb')
  rownames(clst.cell.each.tab) = NULL

  write.csv(clst.cell.each.tab, file = 'ANL001b.merged.clusters.cell.numbers.fractions.csv', row.names = F, quote = F)
}

### WNN weight distributions ###
{
  wnn.weights = sobj.filt@meta.data[, c('treatment', 'seurat_clusters', 'SCT.weight', 'ATAC.weight')]

  hist(wnn.weights$SCT.weight)
  hist(wnn.weights$ATAC.weight)

  hist(wnn.weights$SCT.weight, col = NULL, border = 'blue')
  hist(wnn.weights$ATAC.weight, col = NULL, border = 'red', add = T)
  
  wnn.weights.bycluster = split(wnn.weights[, 3:4], wnn.weights$seurat_clusters)

  sapply(1:length(wnn.weights.bycluster), 
         function(x) boxplot(wnn.weights.bycluster[[x]], 
                             main = paste0('Cluster ', x, ': ', nrow(wnn.weights.bycluster[[x]])))
  )

  sapply(1:length(wnn.weights.bycluster), function(x) {
    hist(wnn.weights.bycluster[[x]]$SCT.weight, col = NULL, border = 'blue', freq = F, 
         xlim = c(0, 1), 
         main = paste0('Cluster ', names(wnn.weights.bycluster)[x], ': ', nrow(wnn.weights.bycluster[[x]])))
    hist(wnn.weights.bycluster[[x]]$ATAC.weight, col = NULL, border = 'red', freq = F, add = T)
  })
  
  sapply(1:length(wnn.weights.bycluster), 
         function(x) plot(wnn.weights.bycluster[[x]]$SCT.weight, 
                          wnn.weights.bycluster[[x]]$ATAC.weight, 
                          xlim = c(0, 1), ylim = c(0, 1), 
                          main = paste0('Cluster ', names(wnn.weights.bycluster)[x], ': ', 
                                        nrow(wnn.weights.bycluster[[x]])))
  )
}

{
  write.csv(sobj.filt[['RNA']]@counts, file = 'ANL001b.merged.data.processed.rna.counts.d210913.csv', quote = F)
  write.csv(sobj.filt[['ATAC']]@counts, file = 'ANL001b.merged.data.processed.atac.counts.d210913.csv', quote = F)
  write.csv(sobj.filt@meta.data, file = 'ANL001b.merged.data.processed.metadata.d210913.csv', quote = F)
}


##### DEA #####
##### Finding differentially expressed features (cluster biomarkers) #####

DefaultAssay(sobj.filt) = 'SCT'

### find markers for every cluster compared to all remaining cells, report only the positive ones
system.time({
  dea.markers = FindAllMarkers(sobj.filt, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
})

write.csv(dea.markers, file = 'ANL001b.merged.dea.results.d210913.csv', quote = F, row.names = F)

### UMAP
p1 = FeaturePlot(sobj.filt, reduction = 'umap.rna', ncol = 9, 
                 features = (dea.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC))$gene)
ggsave(p1, filename = 'ANL001b.merged.dea.top.marker.umap.rna.png', width = 30, height = 15)

p1 = FeaturePlot(sobj.filt, reduction = 'umap.atac', ncol = 9, 
                 features = (dea.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC))$gene)
ggsave(p1, filename = 'ANL001b.merged.dea.top.marker.umap.atac.png', width = 30, height = 15)

p1 = FeaturePlot(sobj.filt, reduction = 'wnn.umap', ncol = 9, 
                 features = (dea.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC))$gene)
ggsave(p1, filename = 'ANL001b.merged.dea.top.marker.umap.wnn.png', width = 30, height = 15)


### Heatmap
dea.top10 = dea.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
p1 = DoHeatmap(sobj.filt, features = dea.top10$gene) + NoLegend()
ggsave(p1, filename = 'ANL001b.merged.dea.top10markers.heatmap.png', width = 20, height = 20)


saveRDS(sobj.filt, file = 'ANL001b.merged.data.processed.d210913.rds')



###################################################################################################

##### Confidence builders 1 #####

###################################################################################################

mygenes = c('si:ch73-1a9.3', 'si:ch73-281n10.2', 'hmgn2', 'hmgn6', 'hmgn7')

#(sobj.filt = readRDS('ANL001b.merged.data.processed.d210913.rds'))

DefaultAssay(sobj.filt) = 'SCT'

### Violin plots
VlnPlot(sobj.filt1, features = mygenes, slot = 'data', split.plot = T, split.by = 'treatment')
VlnPlot(sobj.filt1, features = mygenes, slot = 'counts', log = T) # raw counts

p1 = VlnPlot(sobj.filt1[, names(which(sobj.filt1$treatment == 'WT'))], 
             features = mygenes, slot = 'data', ncol = 5)
p2 = VlnPlot(sobj.filt1[, names(which(sobj.filt1$treatment == 'MUT'))], 
             features = mygenes, slot = 'data', ncol = 5)
p3 = plot_grid(plotlist = list(p1, p2), nrow = 2)
ggsave(p3, filename = 'ANL001b.merged.geneset1.violin.v1.png', width = 40, height = 10)

### UMAP
(p1 = FeaturePlot(sobj.filt1, reduction = 'wnn.umap', features = mygenes, order = T, split.by = 'treatment'))
ggsave(p1, filename = 'ANL001b.merged.geneset1.umap.wnn.v1.png', width = 10, height = 25)

(p1 = FeaturePlot(sobj.filt1, reduction = 'wnn.umap', ncol = 5, features = mygenes, order = T, 
                  cols = c('grey', 'blue'), 
                  cells = names(which(sobj.filt1$treatment == 'WT'))))
(p2 = FeaturePlot(sobj.filt1, reduction = 'wnn.umap', ncol = 5, features = mygenes, order = T, 
                  cols = c('grey', 'dark red'), 
                  cells = names(which(sobj.filt1$treatment == 'MUT'))))
p3 = plot_grid(plotlist = list(p1, p2), nrow = 2)
ggsave(p3, filename = 'ANL001b.merged.geneset1.umap.wnn.v2.png', width = 25, height = 10)

### Heatmap
(p1 = DoHeatmap(sobj.filt1, features = mygenes, slot = 'counts') + 
    scale_fill_gradientn(colors = c("black", "red")) + NoLegend())

p1 = DoHeatmap(sobj.filt1, features = mygenes, slot = 'counts', 
               cells = names(which(sobj.filt1$treatment == 'WT'))) + ggtitle('WT') + 
  scale_fill_gradientn(colors = c("black", "red")) + NoLegend()
p2 = DoHeatmap(sobj.filt1, features = mygenes, slot = 'counts', 
               cells = names(which(sobj.filt1$treatment == 'MUT'))) + ggtitle('Mutant') + 
  scale_fill_gradientn(colors = c("black", "red")) + NoLegend()
p1 + p2
ggsave(p1 + p2, filename = 'ANL001b.merged.geneset1.heatmap.png', width = 20, height = 10)

### Dot plot
DotPlot(sobj.filt1, features = mygenes, cols = c('blue', 'dark red'), split.by = 'treatment') + 
  xlab('') + ylab('Cluster by genotype')# + RotatedAxis()
ggsave(filename = 'ANL001b.merged.geneset1.dots.png', width = 10, height = 20)



###################################################################################################

##### Confidence builders 2 #####

###################################################################################################

mygenes = c('grin1a', 'grin1b')

### Violin plots
sobj.filt1@meta.data$treatment = factor(sobj.filt1@meta.data$treatment, 
                                        levels = c('WT', 'MUT'))

{
  VlnPlot(sobj.filt1, features = mygenes, slot = 'data', split.by = 'treatment', split.plot = T, 
          cols = c('orange', 'blue'))
  VlnPlot(sobj.filt1, features = mygenes, slot = 'data', split.by = 'treatment', split.plot = F, 
          cols = c('orange', 'blue'))
  VlnPlot(sobj.filt1, features = mygenes, slot = 'counts', log = T) # raw counts
}

{
  p.list = VlnPlot(sobj.filt1, features = mygenes, slot = 'data', split.by = 'treatment', split.plot = F, 
                   cols = c('orange', 'blue'), combine = F)
  ggsave(p.list[[1]] + xlab('Cluster'), 
         filename = 'ANL001b.merged.grin1a.violin.v1.d230119.png', width = 10, height = 6)
  ggsave(p.list[[2]] + xlab('Cluster'), 
         filename = 'ANL001b.merged.grin1b.violin.v1.d230119.png', width = 10, height = 6)
}

{
  p1 = VlnPlot(sobj.filt1[, names(which(sobj.filt1$treatment == 'WT'))], 
               features = mygenes, slot = 'data', ncol = length(mygenes))
  p2 = VlnPlot(sobj.filt1[, names(which(sobj.filt1$treatment == 'MUT'))], 
               features = mygenes, slot = 'data', ncol = length(mygenes))
  (p3 = plot_grid(plotlist = list(p1, p2), nrow = 2))
}

### UMAP
(p1 = FeaturePlot(sobj.filt1, reduction = 'wnn.umap', features = mygenes, order = T, split.by = 'treatment'))
ggsave(p1, filename = 'ANL001b.merged.grin1a.grin1b.umap.wnn.v1.d230119.png', width = 10, height = 10)

(p1 = FeaturePlot(sobj.filt1, reduction = 'wnn.umap', 
                  features = mygenes, ncol = length(mygenes), order = T, 
                  cols = c('grey', 'orange'), 
                  cells = names(which(sobj.filt1$treatment == 'WT'))))
(p2 = FeaturePlot(sobj.filt1, reduction = 'wnn.umap', 
                  features = mygenes, ncol = length(mygenes), order = T, 
                  cols = c('grey', 'blue'), 
                  cells = names(which(sobj.filt1$treatment == 'MUT'))))
(p3 = plot_grid(plotlist = list(p1, p2), nrow = 2))
ggsave(p3, filename = 'ANL001b.merged.grin1a.grin1b.umap.wnn.v2.d230119.png', width = 10, height = 10)

### Heatmap
(p1 = DoHeatmap(sobj.filt1, features = mygenes, slot = 'counts') + 
    scale_fill_gradientn(colors = c("black", "red")) + NoLegend())

p1 = DoHeatmap(sobj.filt1, features = mygenes, slot = 'counts', 
               cells = names(which(sobj.filt1$treatment == 'WT'))) + ggtitle('WT') + 
  scale_fill_gradientn(colors = c("black", "orange")) + NoLegend()
p2 = DoHeatmap(sobj.filt1, features = mygenes, slot = 'counts', 
               cells = names(which(sobj.filt1$treatment == 'MUT'))) + ggtitle('Mutant') + 
  scale_fill_gradientn(colors = c("black", "blue")) + NoLegend()
(p3 = plot_grid(plotlist = list(p1, p2), nrow = 2))
ggsave(p3, filename = 'ANL001b.merged.grin1a.grin1b.heatmap.d230119.png', width = 10, height = 10)

### Dot plot
DotPlot(sobj.filt1, features = mygenes, cols = c('orange', 'blue'), split.by = 'treatment') + 
  xlab('') + ylab('Cluster by genotype')# + RotatedAxis()
ggsave(filename = 'ANL001b.merged.grin1a.grin1b.dots.d230119.png', width = 8, height = 16)



###################################################################################################

##### Heatmap of the top 10 DEGs of each cluster for each condition #####

###################################################################################################

#(sobj.filt = readRDS('ANL001b.merged.data.processed.d210913.rds'))
#dea.markers = read.csv('ANL001b.merged.dea.results.d210913.csv')

dea.top10 = dea.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

p1 = DoHeatmap(sobj.filt, features = dea.top10$gene, 
               cells = names(which(sobj.filt$treatment == 'WT'))) + NoLegend()
p2 = DoHeatmap(sobj.filt, features = dea.top10$gene, 
               cells = names(which(sobj.filt$treatment == 'MUT'))) + NoLegend()

ggsave(p1 + p2, filename = 'ANL001b.merged.dea.top10markers.heatmap.each.png', width = 40, height = 20)



###################################################################################################

##### Heatmap of all DEGs of each cluster for each condition #####

###################################################################################################

#(sobj.filt = readRDS('ANL001b.merged.data.processed.d210913.rds'))
#dea.markers = read.csv('ANL001b.merged.dea.results.d210913.csv')

dea.markers = dea.markers %>% group_by(cluster)

p1 = DoHeatmap(sobj.filt, features = dea.markers$gene, 
               cells = names(which(sobj.filt$treatment == 'WT'))) + NoLegend()
p2 = DoHeatmap(sobj.filt, features = dea.markers$gene, 
               cells = names(which(sobj.filt$treatment == 'MUT'))) + NoLegend()

ggsave(p1 + p2, filename = 'ANL001b.merged.dea.allmarkers.heatmap.each.png', 
       width = 60, height = 30, limitsize = F)



###################################################################################################

##### Sub-clustering #####

###################################################################################################

(sobj.filt = readRDS('ANL001b.merged.data.processed.d210913.rds'))

DefaultAssay(sobj.filt) = 'RNA'

(clst.toobig = names(which(clst.cell.tot > 1000)))

system.time({
  sobj.subclst = lapply(clst.toobig, as.null)
  names(sobj.subclst) = clst.toobig
  
  for ( i in 1:length(clst.toobig) ) {
    tmp1 = subset(sobj.filt, subset = seurat_clusters == clst.toobig[i])
    
    ### RNA analysis ###
    DefaultAssay(tmp1) = "RNA"
    tmp1 = SCTransform(tmp1, verbose = F) %>% RunPCA() %>% 
      RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
    dim(tmp1@reductions$umap.rna@cell.embeddings)
    
    ### ATAC analysis ###
    DefaultAssay(tmp1) = "ATAC"
    tmp1 = RunTFIDF(tmp1)
    tmp1 = FindTopFeatures(tmp1, min.cutoff = 'q0')
    tmp1 = RunSVD(tmp1)
    tmp1 = RunUMAP(tmp1, reduction = 'lsi', dims = 2:50, 
                   reduction.name = "umap.atac", reduction.key = "atacUMAP_")
    
    ### WNN graph ###
    tmp1 = FindMultiModalNeighbors(tmp1, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
    tmp1 = RunUMAP(tmp1, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
    
    tmp1 = FindClusters(tmp1, graph.name = "wsnn", algorithm = 3, verbose = F)
    
    sobj.subclst[[i]] = tmp1
  }
})

for ( i in 1:length(sobj.subclst) ) {
  tmp.sobj = sobj.subclst[[i]]
  
  VlnPlot(
    object = tmp.sobj,
    features = c("nCount_RNA", "percent.mt", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
    group.by = 'treatment', 
    ncol = 5,
    pt.size = 0
  )
  
  VlnPlot(
    object = tmp.sobj,
    features = c("nCount_RNA", "percent.mt", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
    group.by = 'seurat_clusters', 
    ncol = 5,
    pt.size = 0
  )
  
  VlnPlot(
    object = tmp.sobj,
    features = c("nCount_RNA", "percent.mt", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
    group.by = 'wsnn_res.0.8', 
    ncol = 5,
    pt.size = 0
  )
  
  {
    p1 = DimPlot(tmp.sobj, reduction = "wnn.umap", group.by = 'treatment', label = F, label.size = 5, repel = T) + 
      ggtitle("WNN")
    p2 = DimPlot(tmp.sobj, reduction = "umap.rna", group.by = 'treatment', label = F, label.size = 5, repel = T) + 
      ggtitle("RNA")
    p3 = DimPlot(tmp.sobj, reduction = "umap.atac", group.by = 'treatment', label = F, label.size = 5, repel = T) + 
      ggtitle("ATAC")
    p1 + p2 + p3 & theme(plot.title = element_text(hjust = 0.5))
    ggsave(paste0('ANL001b.merged.subclust.', names(sobj.subclst)[i], '.umap.by.genotype.png'), 
           width = 24, height = 8)
  }
  
  {
    p1 = DimPlot(tmp.sobj, reduction = "wnn.umap", label = T, label.size = 10, repel = T) + ggtitle("WNN")
    p2 = DimPlot(tmp.sobj, reduction = "umap.rna", label = T, label.size = 10, repel = T) + ggtitle("RNA")
    p3 = DimPlot(tmp.sobj, reduction = "umap.atac", label = T, label.size = 10, repel = T) + ggtitle("ATAC")
    p1 + p2 + p3 & theme(plot.title = element_text(hjust = 0.5)) & NoLegend()
    ggsave(paste0('ANL001b.merged.subclust.', names(sobj.subclst)[i], '.umap.by.cluster.png'), 
           width = 24, height = 8)
  }
  
  {
    DimPlot(tmp.sobj, reduction = "wnn.umap", group.by = 'seurat_clusters', label = T, label.size = 10, repel = T) + 
      NoLegend() + ggtitle(paste0('Sub-clustering ', names(sobj.subclst)[i]))
    ggsave(paste0('ANL001b.merged.subclust.', names(sobj.subclst)[i], '.wnn.umap.by.cluster.v1.png'), 
           width = 8, height = 8)
  }
  
  {
    DimPlot(tmp.sobj, reduction = "wnn.umap", group.by = 'seurat_clusters', label = F) + 
      NoLegend() + ggtitle(paste0('Sub-clustering ', names(sobj.subclst)[i]))
    ggsave(paste0('ANL001b.merged.subclust.', names(sobj.subclst)[i], '.wnn.umap.by.cluster.v2.png'), 
           width = 8, height = 8)
  }
}

{
  sobj.filt.subclst = sobj.filt
  
  foo.subclst = as.character(sobj.filt.subclst@meta.data$seurat_clusters)
  names(foo.subclst) = rownames(sobj.filt.subclst@meta.data)
  
  for ( i in 1:length(sobj.subclst) ) {
    foo.subclst[rownames(sobj.subclst[[i]]@meta.data)] = 
      paste0(names(sobj.subclst)[i], '-', sobj.subclst[[i]]@meta.data$seurat_clusters)
  }
  
  foo.subclst = factor(foo.subclst)
  levels(foo.subclst)
  table(foo.subclst)
  
  hist(sort(table(foo.subclst)))
  plot(sort(table(foo.subclst)))
  
  sobj.filt.subclst@meta.data$seurat_clusters_v2 = foo.subclst
}

{
  (clst.cell.tot = table(sobj.filt.subclst@meta.data$seurat_clusters_v2))
  (genotype.cell.tot = table(sobj.filt.subclst@meta.data$treatment))
  
  clst.cell.each.tab = NULL
  for ( i in 1:length(clst.cell.tot) ) {
    tmp1 = sobj.filt.subclst@meta.data$treatment[sobj.filt.subclst@meta.data$seurat_clusters_v2 == 
                                                   names(clst.cell.tot)[i]]
    tmp2 = table(tmp1)
    tmp3 = c(names(clst.cell.tot)[i], clst.cell.tot[i], 
             tmp2, 
             round(tmp2 / as.numeric(clst.cell.tot[i]), 3), 
             round(tmp2 / genotype.cell.tot, 3))
    
    clst.cell.each.tab = rbind(clst.cell.each.tab, tmp3)
  }
  
  colnames(clst.cell.each.tab) = c('Cluster', 'Total', 'MUT', 'WT', 'MUT_lc', 'WT_lc', 'MUT_gb', 'WT_gb')
  rownames(clst.cell.each.tab) = NULL

  write.csv(clst.cell.each.tab, file = 'ANL001b.merged.subclust.clusters.cell.numbers.fractions.csv', 
            row.names = F, quote = F)
}

saveRDS(sobj.filt.subclst, file = 'ANL001b.merged.subclust.d211221.rds')


##### DEA #####
DefaultAssay(sobj.filt) = 'SCT'
(clusters.all = sort((unique(sobj.filt$seurat_clusters_v2))))
system.time({
  dea.markers = NULL
  for ( i in 1:length(clusters.all) ) {
    tmp = FindMarkers(sobj.filt, group.by = 'seurat_clusters_v2', ident.1 = clusters.all[i], 
                      only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
    dea.markers = rbind(dea.markers, cbind(cluster = clusters.all[i], gene = rownames(tmp), tmp))
  }
})

rownames(dea.markers) = NULL

write.csv(dea.markers, file = 'ANL001b.merged.subclust.dea.results.d211221.csv', quote = F, row.names = F)


### UMAP
(top.genes = data.frame(dea.markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)))
write.csv(top.genes, row.names = F, quote = F, 
          file = 'ANL001b.merged.subclust.dea.top.marker.by.log2fc.csv')

{
  p1 = lapply(1:nrow(top.genes), function(i) 
    FeaturePlot(sobj.filt, reduction = 'umap.rna', features = top.genes$gene[i]) + 
      ggtitle(paste0(top.genes$cluster[i], ': ', top.genes$gene[i]))
  )
  
  p2 = plot_grid(plotlist = p1, ncol = 9)
  
  ggsave(p2, filename = 'ANL001b.merged.subclust.dea.top.marker.umap.rna.png', width = 27, height = 21)
}

{
  p1 = lapply(1:nrow(top.genes), function(i) 
    FeaturePlot(sobj.filt, reduction = 'umap.atac', features = top.genes$gene[i]) + 
      ggtitle(paste0(top.genes$cluster[i], ': ', top.genes$gene[i]))
  )
  
  p2 = plot_grid(plotlist = p1, ncol = 9)
  
  ggsave(p2, filename = 'ANL001b.merged.subclust.dea.top.marker.umap.atac.png', width = 27, height = 21)
}

{
  p1 = lapply(1:nrow(top.genes), function(i) 
    FeaturePlot(sobj.filt, reduction = 'wnn.umap', features = top.genes$gene[i]) + 
      ggtitle(paste0(top.genes$cluster[i], ': ', top.genes$gene[i]))
  )
  
  p2 = plot_grid(plotlist = p1, ncol = 9)
  
  ggsave(p2, filename = 'ANL001b.merged.subclust.dea.top.marker.umap.wnn.png', width = 27, height = 21)
}

### Heatmap
dea.top10 = dea.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
p1 = DoHeatmap(sobj.filt, features = dea.top10$gene, group.by = 'seurat_clusters_v2') + NoLegend()
ggsave(p1, filename = 'ANL001b.merged.subclust.dea.top10markers.heatmap.png', width = 40, height = 30)

### Heatmap of the top 10 DEGs of each cluster for each condition
dea.top10 = dea.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

p1 = DoHeatmap(sobj.filt, features = dea.top10$gene, group.by = 'seurat_clusters_v2', 
               cells = names(which(sobj.filt$treatment == 'WT'))) + NoLegend()
p2 = DoHeatmap(sobj.filt, features = dea.top10$gene, group.by = 'seurat_clusters_v2', 
               cells = names(which(sobj.filt$treatment == 'MUT'))) + NoLegend()

ggsave(p1 + p2, filename = 'ANL001b.merged.subclust.dea.top10markers.heatmap.each.png', 
       width = 60, height = 30, limitsize = F)



###################################################################################################

##### DEA for WT/Mut for each cluster: original clusters #####

###################################################################################################

(sobj.filt = readRDS('ANL001b.merged.data.processed.d210913.rds'))

(clusters.all = sort((unique(sobj.filt$seurat_clusters))))

{ ### WT markers ###
  system.time({
    dea.markers.wt = NULL
    for ( i in 1:length(clusters.all) ) {
      cat(i, '\n')
      tmp = FindMarkers(sobj.filt, ident.1 = 'WT', group.by = 'treatment', subset.ident = clusters.all[i], 
                        only.pos = T, min.pct = 0.2, logfc.threshold = 0.2)
      dea.markers.wt = rbind(dea.markers.wt, cbind(cluster = clusters.all[i], gene = rownames(tmp), tmp))
    }
  })

  rownames(dea.markers.wt) = NULL
}

{ ### Mut markers ###
  system.time({
    dea.markers.mut = NULL
    for ( i in 1:length(clusters.all) ) {
      cat(i, '\n')
      tmp = FindMarkers(sobj.filt, ident.1 = 'MUT', group.by = 'treatment', subset.ident = clusters.all[i], 
                        only.pos = T, min.pct = 0.2, logfc.threshold = 0.2)
      dea.markers.mut = rbind(dea.markers.mut, cbind(cluster = clusters.all[i], gene = rownames(tmp), tmp))
    }
  })

  rownames(dea.markers.mut) = NULL
}

{
  write.csv(dea.markers.wt, file = 'ANL001b.merged.dea.clusters.wt.results.d220322.csv', 
            quote = F, row.names = F)
  write.csv(dea.markers.mut, file = 'ANL001b.merged.dea.clusters.mut.results.d220322.csv', 
            quote = F, row.names = F)
}

{ ### UMAP ###
  (top.genes.wt = data.frame(dea.markers.wt %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)))
  write.csv(top.genes.wt, row.names = F, quote = F, 
            file = 'ANL001b.merged.dea.clusters.wt.top.marker.by.log2fc.d220322.csv')
  
  (top.genes.mut = data.frame(dea.markers.mut %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)))
  write.csv(top.genes.mut, row.names = F, quote = F, 
            file = 'ANL001b.merged.dea.clusters.mut.top.marker.by.log2fc.d220322.csv')
}

{ ### UMAP-RNA ###
  ## WT ##
  p1 = lapply(1:nrow(top.genes.wt), function(i) 
    FeaturePlot(sobj.filt, reduction = 'umap.rna', features = top.genes.wt$gene[i]) + 
      ggtitle(paste0(top.genes.wt$cluster[i], ': ', top.genes.wt$gene[i]))
  )
  
  p2 = plot_grid(plotlist = p1, ncol = 9)
  
  ggsave(p2, filename = 'ANL001b.merged.dea.clusters.wt.top.marker.umap.rna.d220322.png', 
         width = 36, height = 21)

  ## Mut ##
  p1 = lapply(1:nrow(top.genes.mut), function(i) 
    FeaturePlot(sobj.filt, reduction = 'umap.rna', features = top.genes.mut$gene[i]) + 
      ggtitle(paste0(top.genes.mut$cluster[i], ': ', top.genes.mut$gene[i]))
  )
  
  p2 = plot_grid(plotlist = p1, ncol = 9)
  
  ggsave(p2, filename = 'ANL001b.merged.dea.clusters.mut.top.marker.umap.rna.d220322.png', 
         width = 36, height = 21)
}

{ ### UMAP-ATAC ###
  ## WT ##
  p1 = lapply(1:nrow(top.genes.wt), function(i) 
    FeaturePlot(sobj.filt, reduction = 'umap.atac', features = top.genes.wt$gene[i]) + 
      ggtitle(paste0(top.genes.wt$cluster[i], ': ', top.genes.wt$gene[i]))
  )
  
  p2 = plot_grid(plotlist = p1, ncol = 9)
  
  ggsave(p2, filename = 'ANL001b.merged.dea.clusters.wt.top.marker.umap.atac.d220322.png', 
         width = 36, height = 21)

  ## Mut ##
  p1 = lapply(1:nrow(top.genes.mut), function(i) 
    FeaturePlot(sobj.filt, reduction = 'umap.atac', features = top.genes.mut$gene[i]) + 
      ggtitle(paste0(top.genes.mut$cluster[i], ': ', top.genes.mut$gene[i]))
  )
  
  p2 = plot_grid(plotlist = p1, ncol = 9)
  
  ggsave(p2, filename = 'ANL001b.merged.dea.clusters.mut.top.marker.umap.atac.d220322.png', 
         width = 36, height = 21)
}

{ ### UMAP-WNN ###
  ## WT ##
  p1 = lapply(1:nrow(top.genes.wt), function(i) 
    FeaturePlot(sobj.filt, reduction = 'wnn.umap', features = top.genes.wt$gene[i]) + 
      ggtitle(paste0(top.genes.wt$cluster[i], ': ', top.genes.wt$gene[i]))
  )
  
  p2 = plot_grid(plotlist = p1, ncol = 9)
  
  ggsave(p2, filename = 'ANL001b.merged.dea.clusters.wt.top.marker.umap.wnn.d220222.png', 
         width = 36, height = 21)
  
  ## Mut ##
  p1 = lapply(1:nrow(top.genes.mut), function(i) 
    FeaturePlot(sobj.filt, reduction = 'wnn.umap', features = top.genes.mut$gene[i]) + 
      ggtitle(paste0(top.genes.mut$cluster[i], ': ', top.genes.mut$gene[i]))
  )
  
  p2 = plot_grid(plotlist = p1, ncol = 9)
  
  ggsave(p2, filename = 'ANL001b.merged.dea.clusters.mut.top.marker.umap.wnn.d220222.png', 
         width = 36, height = 21)
}

{ ### Heatmap of the top 10 DEGs ###
  dea.top10.wt = dea.markers.wt %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  p1 = DoHeatmap(sobj.filt, features = dea.top10.wt$gene, group.by = 'seurat_clusters') + NoLegend()
  ggsave(p1, filename = 'ANL001b.merged.dea.clusters.wt.top10markers.heatmap.d220322.png', 
         width = 40, height = 30)

  dea.top10.mut = dea.markers.mut %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  p1 = DoHeatmap(sobj.filt, features = dea.top10.mut$gene, group.by = 'seurat_clusters') + NoLegend()
  ggsave(p1, filename = 'ANL001b.merged.dea.clusters.mut.top10markers.heatmap.d220322.png', 
         width = 40, height = 30)
}

{ ### Heatmap of the top 10 DEGs of each cluster for each condition ###
  p1 = DoHeatmap(sobj.filt, features = dea.top10.wt$gene, group.by = 'seurat_clusters', 
                 cells = names(which(sobj.filt$treatment == 'WT'))) + NoLegend()
  p2 = DoHeatmap(sobj.filt, features = dea.top10.wt$gene, group.by = 'seurat_clusters', 
                 cells = names(which(sobj.filt$treatment == 'MUT'))) + NoLegend()

  ggsave(p1 + p2, filename = 'ANL001b.merged.dea.clusters.wt.top10markers.heatmap.each.d220322.png', 
         width = 60, height = 30, limitsize = F)

  p1 = DoHeatmap(sobj.filt, features = dea.top10.mut$gene, group.by = 'seurat_clusters', 
                 cells = names(which(sobj.filt$treatment == 'WT'))) + NoLegend()
  p2 = DoHeatmap(sobj.filt, features = dea.top10.mut$gene, group.by = 'seurat_clusters', 
                 cells = names(which(sobj.filt$treatment == 'MUT'))) + NoLegend()
  
  ggsave(p1 + p2, filename = 'ANL001b.merged.dea.clusters.mut.top10markers.heatmap.each.d220322.png', 
         width = 60, height = 30, limitsize = F)
}



###################################################################################################

##### DEA for WT/Mut for each cluster: sub-clustering #####

###################################################################################################

(sobj.filt = readRDS('ANL001b.merged.subclust.d211221.rds'))

DefaultAssay(sobj.filt) = 'SCT'

(clusters.all = sort((unique(sobj.filt$seurat_clusters_v2))))

{ ### WT markers ###
  system.time({
    dea.markers.wt = NULL
    for ( i in 1:length(clusters.all) ) {
      cat(i, '\n')
      tmp = FindMarkers(sobj.filt, ident.1 = 'WT', group.by = 'treatment', subset.ident = clusters.all[i], 
                        only.pos = T, min.pct = 0.2, logfc.threshold = 0.2)
      dea.markers.wt = rbind(dea.markers.wt, cbind(cluster = clusters.all[i], gene = rownames(tmp), tmp))
    }
  })
  
  rownames(dea.markers.wt) = NULL
}

{ ### Mut markers ###
  system.time({
    dea.markers.mut = NULL
    for ( i in 1:length(clusters.all) ) {
      cat(i, '\n')
      tmp = FindMarkers(sobj.filt, ident.1 = 'MUT', group.by = 'treatment', subset.ident = clusters.all[i], 
                        only.pos = T, min.pct = 0.2, logfc.threshold = 0.2)
      dea.markers.mut = rbind(dea.markers.mut, cbind(cluster = clusters.all[i], gene = rownames(tmp), tmp))
    }
  })

  rownames(dea.markers.mut) = NULL
}

{
  write.csv(dea.markers.wt, file = 'ANL001b.merged.dea.subclusters.wt.results.d220316.csv', quote = F, row.names = F)
  write.csv(dea.markers.mut, file = 'ANL001b.merged.dea.subclusters.mut.results.d220316.csv', quote = F, row.names = F)
}

{ ### UMAP ###
  (top.genes.wt = data.frame(dea.markers.wt %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)))
  write.csv(top.genes.wt, row.names = F, quote = F, 
            file = 'ANL001b.merged.dea.subclusters.wt.top.marker.by.log2fc.d220316.csv')

  (top.genes.mut = data.frame(dea.markers.mut %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)))
  write.csv(top.genes.mut, row.names = F, quote = F, 
            file = 'ANL001b.merged.dea.subclusters.mut.top.marker.by.log2fc.d220316.csv')
}

{ ### UMAP-RNA ###
  ## WT ##
  p1 = lapply(1:nrow(top.genes.wt), function(i) 
    FeaturePlot(sobj.filt, reduction = 'umap.rna', features = top.genes.wt$gene[i]) + 
      ggtitle(paste0(top.genes.wt$cluster[i], ': ', top.genes.wt$gene[i]))
  )
  
  p2 = plot_grid(plotlist = p1, ncol = 9)
  
  ggsave(p2, filename = 'ANL001b.merged.dea.subclusters.wt.top.marker.umap.rna.d220316.png', 
         width = 27, height = 21)

  ## Mut ##
  p1 = lapply(1:nrow(top.genes.mut), function(i) 
    FeaturePlot(sobj.filt, reduction = 'umap.rna', features = top.genes.mut$gene[i]) + 
      ggtitle(paste0(top.genes.mut$cluster[i], ': ', top.genes.mut$gene[i]))
  )
  
  p2 = plot_grid(plotlist = p1, ncol = 9)
  
  ggsave(p2, filename = 'ANL001b.merged.dea.subclusters.mut.top.marker.umap.rna.d220316.png', 
         width = 27, height = 21)
}

{ ### UMAP-ATAC ###
  ## WT ##
  p1 = lapply(1:nrow(top.genes.wt), function(i) 
    FeaturePlot(sobj.filt, reduction = 'umap.atac', features = top.genes.wt$gene[i]) + 
      ggtitle(paste0(top.genes.wt$cluster[i], ': ', top.genes.wt$gene[i]))
  )
  
  p2 = plot_grid(plotlist = p1, ncol = 9)
  
  ggsave(p2, filename = 'ANL001b.merged.dea.subclusters.wt.top.marker.umap.atac.d220316.png', 
         width = 27, height = 21)
  
  ## Mut ##
  p1 = lapply(1:nrow(top.genes.mut), function(i) 
    FeaturePlot(sobj.filt, reduction = 'umap.atac', features = top.genes.mut$gene[i]) + 
      ggtitle(paste0(top.genes.mut$cluster[i], ': ', top.genes.mut$gene[i]))
  )
  
  p2 = plot_grid(plotlist = p1, ncol = 9)
  
  ggsave(p2, filename = 'ANL001b.merged.dea.subclusters.mut.top.marker.umap.atac.d220316.png', 
         width = 27, height = 21)
}

{ ### UMAP-WNN ###
  ## WT ##
  p1 = lapply(1:nrow(top.genes.wt), function(i) 
    FeaturePlot(sobj.filt, reduction = 'wnn.umap', features = top.genes.wt$gene[i]) + 
      ggtitle(paste0(top.genes.wt$cluster[i], ': ', top.genes.wt$gene[i]))
  )
  
  p2 = plot_grid(plotlist = p1, ncol = 9)
  
  ggsave(p2, filename = 'ANL001b.merged.dea.subclusters.wt.top.marker.umap.wnn.d220216.png', 
         width = 27, height = 21)
  
  ## Mut ##
  p1 = lapply(1:nrow(top.genes.mut), function(i) 
    FeaturePlot(sobj.filt, reduction = 'wnn.umap', features = top.genes.mut$gene[i]) + 
      ggtitle(paste0(top.genes.mut$cluster[i], ': ', top.genes.mut$gene[i]))
  )
  
  p2 = plot_grid(plotlist = p1, ncol = 9)
  
  ggsave(p2, filename = 'ANL001b.merged.dea.subclusters.mut.top.marker.umap.wnn.d220216.png', 
         width = 27, height = 21)
}

{ ### Heatmap of the top 10 DEGs ###
  dea.top10.wt = dea.markers.wt %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  p1 = DoHeatmap(sobj.filt, features = dea.top10.wt$gene, group.by = 'seurat_clusters_v2') + NoLegend()
  ggsave(p1, filename = 'ANL001b.merged.dea.subclusters.wt.top10markers.heatmap.d220316.png', 
         width = 40, height = 30)
  
  dea.top10.mut = dea.markers.mut %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  p1 = DoHeatmap(sobj.filt, features = dea.top10.mut$gene, group.by = 'seurat_clusters_v2') + NoLegend()
  ggsave(p1, filename = 'ANL001b.merged.dea.subclusters.mut.top10markers.heatmap.d220316.png', 
         width = 40, height = 30)
}

{ ### Heatmap of the top 10 DEGs of each cluster for each condition ###
  p1 = DoHeatmap(sobj.filt, features = dea.top10.wt$gene, group.by = 'seurat_clusters_v2', 
                 cells = names(which(sobj.filt$treatment == 'WT'))) + NoLegend()
  p2 = DoHeatmap(sobj.filt, features = dea.top10.wt$gene, group.by = 'seurat_clusters_v2', 
                 cells = names(which(sobj.filt$treatment == 'MUT'))) + NoLegend()
  
  ggsave(p1 + p2, filename = 'ANL001b.merged.dea.subclusters.wt.top10markers.heatmap.each.d220316.png', 
         width = 60, height = 30, limitsize = F)
  
  p1 = DoHeatmap(sobj.filt, features = dea.top10.mut$gene, group.by = 'seurat_clusters_v2', 
                 cells = names(which(sobj.filt$treatment == 'WT'))) + NoLegend()
  p2 = DoHeatmap(sobj.filt, features = dea.top10.mut$gene, group.by = 'seurat_clusters_v2', 
                 cells = names(which(sobj.filt$treatment == 'MUT'))) + NoLegend()
  
  ggsave(p1 + p2, filename = 'ANL001b.merged.dea.subclusters.mut.top10markers.heatmap.each.d220316.png', 
         width = 60, height = 30, limitsize = F)
}



###############################################################################################

##### Enriched motifs by cluster #####

###############################################################################################

(sobj.filt = readRDS('ANL001b.merged.data.processed.d210913.rds'))

DefaultAssay(sobj.filt) = 'ATAC'

(pfm.set = getMatrixSet(x = JASPAR2020, opts = list(all = T)))

{
  gr = granges(sobj.filt)
  levels(gr@seqnames) = paste0('chr', levels(gr@seqnames))
  seqlevels(gr) = levels(gr@seqnames)
  seqlevels(seqinfo(gr))
}

{
  peaks.1based = GRangesToString(
    grange = StringToGRanges(
      regions = rownames(sobj.filt),
      sep = c(":", "-"),
      starts.in.df.are.0based = T
    ), sep = c(":", "-")
  )
  
  rownames(sobj.filt@assays$ATAC@counts) = peaks.1based
  rownames(sobj.filt@assays$ATAC@data) = peaks.1based
  
  gr[seqnames(gr) == 'chr12'][1]
  ranges(gr[seqnames(gr) == 'chr12'][1]) = IRanges(1, 1488)
  
  gr[seqnames(gr) == 'chr18'][1]
  ranges(gr[seqnames(gr) == 'chr18'][1]) = IRanges(1, 1079)
  
  gr[seqnames(gr) == 'chr25'][1]
  ranges(gr[seqnames(gr) == 'chr25'][1]) = IRanges(1, 876)
  
  motif.matrix = CreateMotifMatrix(features = gr, pwm = pfm.set, 
                                   genome = BSgenome.Drerio.UCSC.danRer11)
}

(motif.object = CreateMotifObject(data = motif.matrix, pwm = pfm.set))
rownames(motif.object@data) = sub('^chr', '', rownames(motif.object@data))

rownames(sobj.filt@assays$ATAC@counts) = sub(':', '-', rownames(sobj.filt))
rownames(sobj.filt@assays$ATAC@data) = sub(':', '-', rownames(sobj.filt))

sobj.filt = SetAssayData(sobj.filt, assay = 'ATAC', slot = 'motifs', new.data = motif.object)

### chromVar ###
system.time({
  sobj.filt@assays$ATAC@ranges = gr
  
  sobj.filt = RunChromVAR(
    object = sobj.filt, 
    motif.matrix = motif.matrix, 
    genome = BSgenome.Drerio.UCSC.danRer11
  )
  
  saveRDS(sobj.filt, 'ANL001b.merged.data.processed.motif.rds')
})

motif.names = markers_motifs$feature
colnames(markers_rna) = paste0("RNA.", colnames(markers_rna))
colnames(markers_motifs) = paste0("motif.", colnames(markers_motifs))
markers_rna$gene = markers_rna$RNA.feature
markers_motifs$gene = ConvertMotifID(sobj.filt, id = motif.names)

topTFs = function(celltype, padj.cutoff = 1e-2) {
  ctmarkers_rna = dplyr::filter(
    markers_rna, RNA.group == celltype, RNA.padj < padj.cutoff, RNA.logFC > 0) %>% 
    arrange(-RNA.auc)
  ctmarkers_motif = dplyr::filter(
    markers_motifs, motif.group == celltype, motif.padj < padj.cutoff, motif.logFC > 0) %>% 
    arrange(-motif.auc)
  ctmarkers_motif$gene = tolower(ctmarkers_motif$gene)
  top_tfs = inner_join(
    x = ctmarkers_rna[, c('RNA.group', 'gene', 'RNA.auc', 'RNA.pval')], 
    y = ctmarkers_motif[, c('motif.group', 'motif.feature', 'gene', 'motif.auc', 'motif.pval')], 
    by = "gene"
  )
  top_tfs$avg_auc = (top_tfs$RNA.auc + top_tfs$motif.auc) / 2
  top_tfs = arrange(top_tfs, -avg_auc)
  return(top_tfs)
}

### identify top markers in all clusters ###
top.tfs.all = NULL
for ( i in sort(unique(sobj.filt@meta.data$seurat_clusters)) ) {
  top.tfs.all = rbind(top.tfs.all, topTFs(i))
}

write.csv(top.tfs.all, file = 'ANL001b.merged.top.motifs.all.csv', quote = F, row.names = F)


##### WNN-UMAP plots for top motifs/TFs in each cluster #####
### Multiple png files ###
ftn.wnn.umap.toptfs = function(FILENAME, sobj) {
  top.tfs.all = read.csv(FILENAME)
  all.clusters = unique(top.tfs.all$RNA.group)
  for ( tmp.cluster in all.clusters ) {
    cat(tmp.cluster, '\n')
    (top.tfs.sub = top.tfs.all[top.tfs.all$RNA.group == tmp.cluster, ])
    
    tmp.title = paste0('Cluster ', tmp.cluster)
    for ( i in 1:nrow(top.tfs.sub) ) {
      tmp.gene = top.tfs.sub$gene[i]
      tmp.motif = top.tfs.sub$motif.feature[i]
      
      DefaultAssay(sobj) = 'SCT'
      gene_plot = FeaturePlot(sobj, features = paste0('sct_', tmp.gene), 
                              reduction = 'wnn.umap', pt.size = .1) + ggtitle(tmp.gene)
      
      DefaultAssay(sobj) = 'chromvar'
      motif_plot = FeaturePlot(sobj, features = tmp.motif, min.cutoff = 0, pt.size = .1, 
                               cols = c("lightgrey", "darkred"), reduction = 'wnn.umap')
      
      LabelClusters(gene_plot, id = 'ident', clusters = as.character(tmp.cluster), size = 10, repel = F) + 
        LabelClusters(motif_plot, id = 'ident', clusters = as.character(tmp.cluster), size = 10, repel = F) + 
        plot_layout(ncol = 2) + 
        plot_annotation(title = tmp.title, 
                        theme = theme(plot.title = element_text(size = 25, hjust = .5, vjust = .5)))
      
      ggsave(paste0('wnn.umap.top.motifs.tfs/merged/wnnumap.merged.cluster', tmp.cluster, '.', 
                    tmp.gene, '.', tmp.motif, '.png'), width = 12, height = 6)
    }
  }
}

system.time({
  ftn.wnn.umap.toptfs('ANL001b.merged.top.motifs.all.csv', sobj.filt)
})


### A single pdf file ###
ftn.wnn.umap.toptfs.plots.list = function(FILENAME, sobj) {
  top.tfs.all = read.csv(FILENAME)
  all.clusters = unique(top.tfs.all$RNA.group)
  plots.list = list()
  for ( tmp.cluster in all.clusters ) {
    cat(tmp.cluster, '\n')
    (top.tfs.sub = top.tfs.all[top.tfs.all$RNA.group == tmp.cluster, ])
    
    tmp.title = paste0('Cluster ', tmp.cluster)
    tmp.list = list()
    for ( i in 1:nrow(top.tfs.sub) ) {
      tmp.gene = top.tfs.sub$gene[i]
      tmp.motif = top.tfs.sub$motif.feature[i]
      
      DefaultAssay(sobj) = 'SCT'
      gene_plot = FeaturePlot(sobj, features = paste0('sct_', tmp.gene), 
                              reduction = 'wnn.umap', pt.size = .1) + ggtitle(tmp.gene)
      
      DefaultAssay(sobj) = 'chromvar'
      motif_plot = FeaturePlot(sobj, features = tmp.motif, min.cutoff = 0, pt.size = .1, 
                               cols = c("lightgrey", "darkred"), reduction = 'wnn.umap')
      
      tmp.p = 
        LabelClusters(gene_plot, id = 'ident', clusters = as.character(tmp.cluster), size = 10, repel = F) + 
        LabelClusters(motif_plot, id = 'ident', clusters = as.character(tmp.cluster), size = 10, repel = F) + 
        plot_layout(ncol = 2) + 
        plot_annotation(title = tmp.title, 
                        theme = theme(plot.title = element_text(size = 25, hjust = .5, vjust = .5)))
      
      tmp.list[[i]] = tmp.p
    }
    
    plots.list = c(plots.list, tmp.list)
  }
  
  return(plots.list)
}

system.time({
  pdf('wnn.umap.top.motifs.tfs/wnnumap.merged.all.pdf', width = 12, height = 6)
  invisible(lapply(ftn.wnn.umap.toptfs.plots.list('ANL001b.merged.top.motifs.all.csv', sobj.filt), print))
  dev.off()
})

all.genes = sort(unique(top.tfs.all[, 'gene']))

### List of individual plots ###
system.time({
  plots.list = list()
  for ( i in 1:length(all.genes) ) {
    cat(i, '\n')
    
    if ( i == 98 ) next
    
    g = all.genes[i]
    tmp.anno = Annotation(sobj.filt[['ATAC']])
    tmp.anno = tmp.anno[tmp.anno$gene_name == g]
    
    p = CoveragePlot(sobj.filt, region = g, features = g, peaks = F, 
                     assay = 'ATAC', expression.assay = 'SCT', 
                     extend.upstream = 5000, extend.downstream = 5000) +
      plot_annotation(title = g, 
                      theme = theme(plot.title = element_text(size = 25, hjust = .5, vjust = .5)))

    plots.list[[g]] = p
  }
})

### Individual files ###
system.time({
  for ( i in 1:length(plots.list) ) {
    cat(i, '\n')
    ggsave(plot = plots.list[[i]], width = 10, height = 12, 
           paste0('peak.coverage/merged/peak.coverage.merged.', names(plots.list)[i], '.png'))
  }
})

### A single pdf file ###
system.time({
  pdf('peak.coverage/peak.coverage.merged.all.pdf', width = 10, height = 12)
  invisible(lapply(plots.list, print))
  dev.off()
})



###############################################################################################

##### Enriched motifs by condition #####

###############################################################################################

#(sobj.filt = readRDS('ANL001b.merged.data.processed.motif.rds'))

DefaultAssay(sobj.filt) = 'ATAC'
Idents(sobj.filt) = 'treatment'

data.chva = sobj.filt@assays$chromvar@data

markers_motifs = presto:::wilcoxauc.Seurat(X = sobj.filt, group_by = 'treatment', 
                                           assay = 'data', seurat_assay = 'chromvar')
markers_motifs$gene = ConvertMotifID(sobj.filt, id = markers_motifs$feature)
markers_motifs = markers_motifs[order(markers_motifs$padj), ]

top.motifs = markers_motifs[markers_motifs$padj < 1e-10, ]

{
  ks.tests = NULL
  for ( i in 1:nrow(data.chva) ) {
    tmp1 = data.chva[i, grep('WT', colnames(data.chva))]
    tmp2 = data.chva[i, grep('MUT', colnames(data.chva))]
    
    tmp3 = ks.test(tmp1, tmp2, alternative = 'l')$p.value
    tmp4 = ks.test(tmp1, tmp2, alternative = 'g')$p.value
    
    tmp5 = mean(tmp1)
    tmp6 = mean(tmp2)
    
    ks.tests = rbind(ks.tests, c(tmp3, tmp4, tmp5, tmp6, if (tmp5 > tmp6) 'WT' else 'Mut'))
  }
}

### identify top markers in all clusters ###
top.tfs.all = NULL
for ( i in sort(unique(sobj.filt@meta.data$treatment)) ) {
  top.tfs.all = rbind(top.tfs.all, topTFs(i))
}

write.csv(top.tfs.all, file = 'ANL001b.merged.top.motifs.by.cond.all.d2211.csv', quote = F, row.names = F)


##### WNN-UMAP plots for top motifs/TFs in each cluster #####
### Multiple png files ###
ftn.wnn.umap.toptfs = function(FILENAME, sobj) {
  top.tfs.all = read.csv(FILENAME)
  all.clusters = unique(top.tfs.all$RNA.group)
  for ( tmp.cluster in all.clusters ) {
    cat(tmp.cluster, '\n')
    (top.tfs.sub = top.tfs.all[top.tfs.all$RNA.group == tmp.cluster, ])
    
    tmp.title = paste0('Cluster ', tmp.cluster)
    for ( i in 1:nrow(top.tfs.sub) ) {
      tmp.gene = top.tfs.sub$gene[i]
      tmp.motif = top.tfs.sub$motif.feature[i]
      
      DefaultAssay(sobj) = 'SCT'
      gene_plot = FeaturePlot(sobj, features = paste0('sct_', tmp.gene), 
                              reduction = 'wnn.umap', pt.size = .1) + ggtitle(tmp.gene)
      
      DefaultAssay(sobj) = 'chromvar'
      motif_plot = FeaturePlot(sobj, features = tmp.motif, min.cutoff = 0, pt.size = .1, 
                               cols = c("lightgrey", "darkred"), reduction = 'wnn.umap')
      
      LabelClusters(gene_plot, id = 'ident', clusters = as.character(tmp.cluster), size = 10, repel = F) + 
        LabelClusters(motif_plot, id = 'ident', clusters = as.character(tmp.cluster), size = 10, repel = F) + 
        plot_layout(ncol = 2) + 
        plot_annotation(title = tmp.title, 
                        theme = theme(plot.title = element_text(size = 25, hjust = .5, vjust = .5)))
      
      ggsave(paste0('wnn.umap.top.motifs.tfs/merged/wnnumap.merged.cluster', tmp.cluster, '.', 
                    tmp.gene, '.', tmp.motif, '.png'), width = 12, height = 6)
    }
  }
}

system.time({
  ftn.wnn.umap.toptfs('ANL001b.merged.top.motifs.all.csv', sobj.filt)
})


### A single pdf file ###
ftn.wnn.umap.toptfs.plots.list = function(FILENAME, sobj) {
  top.tfs.all = read.csv(FILENAME)
  all.clusters = unique(top.tfs.all$RNA.group)
  plots.list = list()
  for ( tmp.cluster in all.clusters ) {
    cat(tmp.cluster, '\n')
    (top.tfs.sub = top.tfs.all[top.tfs.all$RNA.group == tmp.cluster, ])
    
    tmp.title = paste0('Cluster ', tmp.cluster)
    tmp.list = list()
    for ( i in 1:nrow(top.tfs.sub) ) {
      tmp.gene = top.tfs.sub$gene[i]
      tmp.motif = top.tfs.sub$motif.feature[i]
      
      DefaultAssay(sobj) = 'SCT'
      gene_plot = FeaturePlot(sobj, features = paste0('sct_', tmp.gene), 
                              reduction = 'wnn.umap', pt.size = .1) + ggtitle(tmp.gene)
      
      DefaultAssay(sobj) = 'chromvar'
      motif_plot = FeaturePlot(sobj, features = tmp.motif, min.cutoff = 0, pt.size = .1, 
                               cols = c("lightgrey", "darkred"), reduction = 'wnn.umap')
      
      tmp.p = 
        LabelClusters(gene_plot, id = 'ident', clusters = as.character(tmp.cluster), size = 10, repel = F) + 
        LabelClusters(motif_plot, id = 'ident', clusters = as.character(tmp.cluster), size = 10, repel = F) + 
        plot_layout(ncol = 2) + 
        plot_annotation(title = tmp.title, 
                        theme = theme(plot.title = element_text(size = 25, hjust = .5, vjust = .5)))
      
      tmp.list[[i]] = tmp.p
    }
    
    plots.list = c(plots.list, tmp.list)
  }
  
  return(plots.list)
}

system.time({
  pdf('wnn.umap.top.motifs.tfs/wnnumap.merged.all.pdf', width = 12, height = 6)
  invisible(lapply(ftn.wnn.umap.toptfs.plots.list('ANL001b.merged.top.motifs.all.csv', sobj.filt), print))
  dev.off()
})

all.genes = sort(unique(top.tfs.all[, 'gene']))

### List of individual plots ###
system.time({
  plots.list = list()
  for ( i in 1:length(all.genes) ) {
    cat(i, '\n')
    
    if ( i == 98 ) next
    
    g = all.genes[i]
    tmp.anno = Annotation(sobj.filt[['ATAC']])
    tmp.anno = tmp.anno[tmp.anno$gene_name == g]
    
    p = CoveragePlot(sobj.filt, region = g, features = g, peaks = F, show.bulk = T, 
                     assay = 'ATAC', expression.assay = 'SCT', 
                     extend.upstream = 5000, extend.downstream = 5000) +
      plot_annotation(title = g, 
                      theme = theme(plot.title = element_text(size = 25, hjust = .5, vjust = .5)))

    plots.list[[g]] = p
  }
})

### Individual files ###
system.time({
  for ( i in 1:length(plots.list) ) {
    cat(i, '\n')
    ggsave(plot = plots.list[[i]], width = 10, height = 12, 
           paste0('peak.coverage/merged/peak.coverage.merged.', names(plots.list)[i], '.png'))
  }
})

### A single pdf file ###
system.time({
  pdf('peak.coverage/peak.coverage.merged.all.pdf', width = 10, height = 12)
  invisible(lapply(plots.list, print))
  dev.off()
})



###################################################################################################

##### Venn diagrams of ATAC data #####

###################################################################################################

#(sobj.filt = readRDS('ANL001b.merged.data.processed.motif.rds'))

DefaultAssay(sobj.filt) = "ATAC"

atac.cnt = data.matrix(sobj.filt[['ATAC']]@counts)
atac.norm = data.matrix(sobj.filt[['ATAC']]@data)

cnt.vals = c(atac.cnt)
norm.vals = c(atac.norm)

{
  cnt1 = atac.cnt[, grep('WT', colnames(atac.cnt))]
  cnt2 = atac.cnt[, grep('MUT', colnames(atac.cnt))]
  
  cnt1.vals = c(cnt1)
  cnt2.vals = c(cnt2)
}

{
  norm1 = atac.norm[, grep('WT', colnames(atac.norm))]
  norm2 = atac.norm[, grep('MUT', colnames(atac.norm))]
  
  norm1.vals = c(norm1)
  norm2.vals = c(norm2)
}

{
  cnt.means = cbind(WT = rowMeans(cnt1), MUT = rowMeans(cnt2))
  norm.means = cbind(WT = rowMeans(norm1), MUT = rowMeans(norm2))
  norm.means1 = cbind(WT = apply(norm1, 1, function(x) if ( any(x > 0) ) mean(x[x > 0]) else 0), 
                      MUT = apply(norm2, 1, function(x) if ( any(x > 0) ) mean(x[x > 0]) else 0))
}

{
  (QUANTILE.THRES = seq(0.1, 0.9, 0.1))
  (THRES.PEAKS = sapply(QUANTILE.THRES, function(x) quantile(norm.means[norm.means > 0], probs = x)))

  norm.means.bin.list = list()
  for ( i in THRES.PEAKS ) {
    norm.means.bin = norm.means
    norm.means.bin[norm.means.bin > i] = 1
    norm.means.bin[norm.means.bin <= i] = 0
    norm.means.bin.list[[names(which(THRES.PEAKS == i))]] = norm.means.bin
  }

  (norm.means.bin.corr.thres = sapply(norm.means.bin.list, function(x) round(cor(x)[1, 2], 2)))

  (norm.means.bin.manhat.thres = sapply(norm.means.bin.list, function(x) 
    round(dist(t(x), method = 'manhattan') / nrow(x), 2)))

  (norm.means.corr.thres = sapply(norm.means.bin.list, function(x) {
    y = norm.means[rowSums(x) > 0, ]; round(cor(y)[1, 2], 2)
  }))

  (peak.patt.thres = sapply(norm.means.bin.list, function(y) 
    table(apply(y, 1, function(x) paste(x, collapse = '-')))))

  (peak.patt.thres.pct = round(peak.patt.thres / nrow(norm.means), 3) * 100)

  save(QUANTILE.THRES, THRES.PEAKS, norm.means, norm.means.bin.list, 
       norm.means.bin.corr.thres, norm.means.bin.manhat.thres, norm.means.corr.thres, peak.patt.thres, 
       file = 'ANL001b.merged.peaks.binary.signals.all.thres.rdata')
}

{
  write.csv(peak.patt.thres.pct, file = 'ANL001b.merged.atac.peaks.patt.all.thres.pct.csv')
  
  png(file = 'ANL001b.merged.atac.peaks.patt.all.thres.pct.heatmap.png', width = 800, height = 800)
  par(omi = c(.2, 0, 0, .2), cex.main = 2)
  heatmap.2(peak.patt.thres.pct, trace = 'none', 
            keysize = .9, key.title = '', key.xlab = 'Percentage', 
            col = colorpanel(10, 'light grey', 'blue'), 
            Colv = NA, dendrogram = 'row', 
            cexRow = 2, cexCol = 2, 
            cellnote = peak.patt.thres.pct, notecex = 2, notecol = 'black', 
            main = 'ATAC peak patterns')
  dev.off()
}

{
  tmp1 = norm.means.bin.list$`50%`
  
  tmp2 = names(which(apply(tmp1, 1, function(x) all(x == c(1, 0)))))
  tmp3 = names(which(apply(tmp1, 1, function(x) all(x == c(0, 1)))))
  
  peak.sig1 = norm.means[tmp2, ]
  peak.sig2 = norm.means[tmp3, ]
  
  peak.sig1 = cbind(t(matrix(as.numeric(unlist(strsplit(rownames(peak.sig1), '[-]'))), nrow = 3, byrow = F)), 
                    peak.sig1)
  colnames(peak.sig1)[1:3] = c('Chr', 'Start', 'End')
  
  peak.sig2 = cbind(t(matrix(as.numeric(unlist(strsplit(rownames(peak.sig2), '[-]'))), nrow = 3, byrow = F)), 
                    peak.sig2)
  colnames(peak.sig2)[1:3] = c('Chr', 'Start', 'End')
  
  peak.sig1 = peak.sig1[order(peak.sig1[, 'WT'], decreasing = T), ]

  peak.sig2 = peak.sig2[order(peak.sig2[, 'MUT'], decreasing = T), ]

  write.csv(peak.sig1, file = 'ANL001b.merged.peaks.signal.top50pct.allcells.mean.intensities.wt.specific.csv', quote = F)
  write.csv(peak.sig2, file = 'ANL001b.merged.peaks.signal.top50pct.allcells.mean.intensities.mut.specific.csv', quote = F)
  
  par(omi = c(.6, 0, 0, 0), cex.main = 1.6)
  heatmap.2(data.matrix(peak.sig1[1:1000, c('WT', 'MUT')]), trace = 'none', 
            keysize = .9, key.xlab = 'Peak intensity', key.par = list(cex.main = 1), 
            col = colorpanel(10, 'grey', 'blue'), 
            labRow = '', 
            Colv = F, Rowv = F, dendrogram = 'none', 
            main = 'WT-specific top 1,000 peaks')
  heatmap.2(data.matrix(peak.sig2[1:1000, c('WT', 'MUT')]), trace = 'none', 
            keysize = .9, key.xlab = 'Peak intensity', key.par = list(cex.main = 1), 
            col = colorpanel(10, 'grey', 'blue'), 
            labRow = '', 
            Colv = F, Rowv = F, dendrogram = 'none', 
            main = 'Mutant-specific top 1,000 peaks')
  heatmap.2(data.matrix(rbind(peak.sig1[1:1000, c('WT', 'MUT')], peak.sig2[1:1000, c('WT', 'MUT')])), 
            trace = 'none', keysize = .9, key.xlab = 'Peak intensity', key.par = list(cex.main = 1), 
            col = colorpanel(10, 'grey', 'blue'), 
            labRow = '', 
            Colv = F, Rowv = F, dendrogram = 'none', 
            main = 'WT/Mut-specific top 1,000 peaks')
}

{
  tmp.medi = median(norm.means[norm.means > 0])
  norm.means.bin = norm.means
  norm.means.bin[norm.means.bin > tmp.medi] = 1
  norm.means.bin[norm.means.bin <= tmp.medi] = 0

  peak.patt = apply(norm.means.bin, 1, function(x) paste(x, collapse = '-'))

  my.patterns = c('1-0', '0-1', '1-1')
  
  sigs.all = lapply(my.patterns, function(x) names(which(peak.patt == x)))
  names(sigs.all) = my.patterns
}

{
  tmp1 = apply(norm.means.bin, 2, function(x) names(which(x == 1)))
  tmp2 = lapply(tmp1, function(x) matrix(unlist(strsplit(x, '-')), ncol = 3, byrow = T))
  tmp2$WT[, 1] = paste0('chr', tmp2$WT[, 1])
  tmp2$MUT[, 1] = paste0('chr', tmp2$MUT[, 1])
  
  write.table(tmp2$WT, file = 'ANL001b.merged.atac.peaks.thres50.coordis.wt.all.tsv', 
              quote = F, row.names = F, col.names = F, sep = '\t')
  write.table(tmp2$MUT, file = 'ANL001b.merged.atac.peaks.thres50.coordis.mut.all.tsv', 
              quote = F, row.names = F, col.names = F, sep = '\t')
  
  sigs.mtx.all = lapply(sigs.all, function(x) {
    tmp = matrix(unlist(strsplit(x, '-')), ncol = 3, byrow = T)
    tmp[, 1] = paste0('chr', tmp[, 1])
    tmp
  })

  for ( i in 1:length(sigs.mtx.all) ) {
    tmp1 = paste0('ANL001b.merged.atac.peaks.thres50.coordis.patt', gsub('-', '', names(sigs.mtx.all)[i]), '.tsv')
    write.table(sigs.mtx.all[[i]], file = tmp1, 
                quote = F, row.names = F, col.names = F, sep = '\t')
  }
}

{
  (peaks.by.chr = sapply(sigs.mtx.all, function(x) table(x[, 1])))
  (peaks.by.chr = peaks.by.chr[order(as.numeric(sub('chr', '', rownames(peaks.by.chr)))), ])
  
  (peaks.by.chr.pct = round(peaks.by.chr / rowSums(peaks.by.chr) * 100, 1))
  
  write.table(peaks.by.chr.pct, file = 'ANL001b.merged.atac.peaks.thres50.patt.by.chr.pct.tsv', 
              quote = F, row.names = T, col.names = T, sep = '\t')

  png(file = 'ANL001b.merged.atac.peaks.thres50.patt.by.chr.pct.heatmap.png', width = 800, height = 1000)
  par(omi = c(.2, 0, 0, .2), cex.main = 2)
  heatmap.2(peaks.by.chr.pct, trace = 'none', 
            keysize = .9, key.title = '', key.xlab = 'Percentage', 
            col = colorpanel(10, 'light grey', 'blue'), 
            cexRow = 2, cexCol = 2, 
            cellnote = peaks.by.chr.pct, notecex = 1.5, notecol = 'black', 
            main = 'ATAC peak patterns (top 50%)')
  dev.off()
}

{
  annot = data.frame(Annotation(sobj.filt))
  annot$ID1 = apply(annot[, 1:3], 1, function(x) gsub(' +', '', paste(x, collapse = '-')))
  
  annot.unq = annot[!duplicated(annot$ID1), ]

  sig1 = sigs.all$`1-0`
  sig2 = sigs.all$`0-1`
  sig3 = sigs.all$`1-1`
}

{ ### Peak signals for WT: 1-0 ###
  sig1.tab = t(sapply(strsplit(sig1, '-'), function(x) as.numeric(x)))

  system.time({
    FLANKING = 5000
    tmp.cols = NULL
    for ( i in 1:nrow(sig1.tab) ) {
      cat(i, '\n')
      
      y1 = rep(0, nrow(annot.unq))
      y2 = which(annot.unq$seqnames == sig1.tab[i, 1])
      y1[y2] = as.numeric(sig1.tab[i, 2] >= annot.unq$start[y2] - FLANKING & sig1.tab[i, 3] <= annot.unq$end[y2] + FLANKING)
      
      if ( sum(y1) == 0 ) next
      
      tmp.cols = c(tmp.cols, sig1[i])
      
      cat(y1, file = 'tmp.values', sep = '\n', append = T)
    }
    
    cat(tmp.cols, file = 'tmp.cols', sep = '\n') # For backup
  })

  system.time({
    tmp.values = scan('tmp.values')
  })

  system.time({
    peak.annot.mtx1 = matrix(tmp.values, nr = nrow(annot.unq), nc = length(tmp.cols), byrow = F)
    rownames(peak.annot.mtx1) = annot.unq$ID1
    colnames(peak.annot.mtx1) = tmp.cols
  })

  system.time({
    signalmap1.list = apply(peak.annot.mtx1, 2, function(x) annot.unq[x == 1, ])
  })

  (signal.genes.freq1 = sort(table(unlist(sapply(signalmap1.list, function(x) unique(x$gene_name)))), decreasing = T))[1:5]

  par(omi = c(0, .2, 0, 0), cex = 1.5)
  plot(sort(as.numeric(table(unlist(sapply(signalmap1.list, function(x) unique(x$gene_name))))), decreasing = T), 
       ylab = 'Number of overlapping signal peaks', xlab = 'Genes (ranked)', main = 'WT-specific top 50% signal peaks', 
       log = 'xy')
  legend('topright', bty = 'n', cex = 1, 
         legend = paste0(1:5, '. ', names(signal.genes.freq1)[1:5], ': ', signal.genes.freq1[1:5]))
  
  signalmap1.genes = sapply(signalmap1.list, function(x) paste(unique(x$gene_name), collapse = '|'))

  system.time({
    save(peak.annot.mtx1, signalmap1.list, signalmap1.genes, 
         file = 'ANL001b.merged.atac.peak.signals.top50pct.wt.annot.rdata')
  })

  file.remove('tmp.values', 'tmp.cols')

  peak.sig.int1 = read.csv('ANL001b.merged.peaks.signal.top50pct.allcells.mean.intensities.wt.specific.csv', row.names = 1)

  peak.sig.int1.genes = cbind(peak.sig.int1, Genes = signalmap1.genes[rownames(peak.sig.int1)])

  write.csv(peak.sig.int1.genes, quote = F, 
            file = 'ANL001b.merged.peaks.signal.top50pct.allcells.mean.intensities.wt.specific.genes.csv')
}

{ ### Peak signals for Mut: 0-1 ###
  sig2.tab = t(sapply(strsplit(sig2, '-'), function(x) as.numeric(x)))

  system.time({
    FLANKING = 5000
    tmp.cols = NULL
    for ( i in 1:nrow(sig2.tab) ) {
      cat(i, '\n')
      
      y1 = rep(0, nrow(annot.unq))
      y2 = which(annot.unq$seqnames == sig2.tab[i, 1])
      y1[y2] = as.numeric(sig2.tab[i, 2] >= annot.unq$start[y2] - FLANKING & sig2.tab[i, 3] <= annot.unq$end[y2] + FLANKING)

      if ( sum(y1) == 0 ) next
      
      tmp.cols = c(tmp.cols, sig2[i])
      
      cat(y1, file = 'tmp.values', sep = '\n', append = T)
    }
    
    cat(tmp.cols, file = 'tmp.cols', sep = '\n') # For backup
  })

  system.time({
    tmp.values = scan('tmp.values')
  })

  system.time({
    peak.annot.mtx2 = matrix(tmp.values, nr = nrow(annot.unq), nc = length(tmp.cols), byrow = F)
    rownames(peak.annot.mtx2) = annot.unq$ID1
    colnames(peak.annot.mtx2) = tmp.cols
  })

  system.time({
    signalmap2.list = apply(peak.annot.mtx2, 2, function(x) annot.unq[x == 1, ])
  })

  (signal.genes.freq2 = sort(table(unlist(sapply(signalmap2.list, function(x) unique(x$gene_name)))), decreasing = T))[1:5]

  par(omi = c(0, .2, 0, 0), cex = 1.5)
  plot(sort(as.numeric(table(unlist(sapply(signalmap2.list, function(x) unique(x$gene_name))))), decreasing = T), 
       ylab = 'Number of overlapping signal peaks', xlab = 'Genes (ranked)', main = 'Mutant-specific top 50% signal peaks', 
       log = 'xy')
  legend('topright', bty = 'n', cex = 1, 
         legend = paste0(1:5, '. ', names(signal.genes.freq2)[1:5], ': ', signal.genes.freq2[1:5]))
  
  signalmap2.genes = sapply(signalmap2.list, function(x) paste(unique(x$gene_name), collapse = '|'))

  system.time({
    save(peak.annot.mtx2, signalmap2.list, signalmap2.genes, 
         file = 'ANL001b.merged.atac.peak.signals.top50pct.mut.annot.rdata')
  })

  file.remove('tmp.values', 'tmp.cols')

  peak.sig.int2 = read.csv('ANL001b.merged.peaks.signal.top50pct.allcells.mean.intensities.mut.specific.csv', row.names = 1)

  peak.sig.int2.genes = cbind(peak.sig.int2, Genes = signalmap2.genes[rownames(peak.sig.int2)])

  write.csv(peak.sig.int2.genes, quote = F, 
            file = 'ANL001b.merged.peaks.signal.top50pct.allcells.mean.intensities.mut.specific.genes.csv')
}

{ ### Comparison between WT and Mut signal genes ###
  signal.genes.cmm = intersect(names(signal.genes.freq1), names(signal.genes.freq2))
  deg = read.csv('ANL001b.merged.dea.results.d210913.csv')

  length(intersect(names(signal.genes.freq1), deg$gene))
  length(intersect(names(signal.genes.freq2), deg$gene))
}



###################################################################################################

##### Condition-specific peaks #####

###################################################################################################

{
  #(sobj.filt = readRDS('ANL001b.merged.data.processed.motif.rds'))

  annot = data.frame(Annotation(sobj.filt))
  annot$ID1 = apply(annot[, 1:3], 1, function(x) gsub(' +', '', paste(x, collapse = '-')))
  
  annot.unq = annot[!duplicated(annot$ID1), ]
}

#(load('ANL001b.merged.peaks.binary.signals.all.thres.rdata'))

{
  peaks.unq1 = names(which(apply(norm.means, 1, function(x) x['WT'] > 0 & x['MUT'] == 0)))
  peaks.unq2 = names(which(apply(norm.means, 1, function(x) x['WT'] == 0 & x['MUT'] > 0)))

  peaks.each = apply(norm.means, 2, function(x) names(which(x > 0)))

  VennDiagram::venn.diagram(peaks.each, filename = 'peaks.norm.means.venn.png', 
                            width = 600, height = 600, 
                            resolution = 100, compression = 'jpeg')
}

{ ### WT-unique peaks ###
  peaks.unq.tab1 = t(sapply(strsplit(peaks.unq1, '-'), function(x) as.numeric(x)))
  rownames(peaks.unq.tab1) = peaks.unq1
  colnames(peaks.unq.tab1) = c('Chr', 'Start', 'End')

  system.time({
    FLANKING = 5000
    tmp.cols = NULL
    for ( i in 1:nrow(peaks.unq.tab1) ) {
      cat(i, '\n')
      
      y1 = rep(0, nrow(annot.unq))
      y2 = which(annot.unq$seqnames == peaks.unq.tab1[i, 'Chr'])
      y1[y2] = as.numeric(peaks.unq.tab1[i, 'Start'] >= annot.unq$start[y2] - FLANKING & 
                            peaks.unq.tab1[i, 'End'] <= annot.unq$end[y2] + FLANKING)
      
      if ( sum(y1) == 0 ) next
      
      tmp.cols = c(tmp.cols, peaks.unq1[i]) # peaks with gene annotations
      
      cat(y1, file = 'tmp.values', sep = '\n', append = T)
    }
    
    cat(tmp.cols, file = 'tmp.cols', sep = '\n') # For backup
  })

  system.time({
    tmp.values = scan('tmp.values')
  })

  system.time({
    peak.annot.mtx1 = matrix(tmp.values, nr = nrow(annot.unq), nc = length(tmp.cols), byrow = F)
    rownames(peak.annot.mtx1) = annot.unq$ID1
    colnames(peak.annot.mtx1) = tmp.cols
  })

  system.time({
    signalmap1.list = apply(peak.annot.mtx1, 2, function(x) annot.unq[x == 1, ])
  })

  (signal.genes.freq1 = sort(table(unlist(sapply(signalmap1.list, 
                                                 function(x) unique(x$gene_name)))), decreasing = T))[1:5]

  par(omi = c(0, .2, 0, 0), cex = 1.5)
  plot(sort(as.numeric(table(unlist(sapply(signalmap1.list, function(x) unique(x$gene_name))))), decreasing = T), 
       ylab = 'Number of overlapping peaks', xlab = 'Genes (ranked)', 
       main = 'WT-specific peaks', 
       log = 'xy')
  legend('topright', bty = 'n', cex = 1, 
         legend = paste0(1:5, '. ', names(signal.genes.freq1)[1:5], ': ', signal.genes.freq1[1:5]))
  
  signalmap1.genes = sapply(signalmap1.list, function(x) paste(unique(x$gene_name), collapse = '|'))

  system.time({
    save(peak.annot.mtx1, signalmap1.list, signalmap1.genes, 
         file = 'ANL001b.wt.mut.peaks.unique.wt.annot.rdata')
  })

  file.remove('tmp.values', 'tmp.cols')

  peaks.unq1.int.genes = cbind(peaks.unq.tab1, 
                               norm.means[rownames(peaks.unq.tab1), ], 
                               Genes = signalmap1.genes[rownames(peaks.unq.tab1)])

  write.csv(peaks.unq1.int.genes, quote = F, 
            file = 'ANL001b.wt.mut.peaks.unqiue.wt.allcells.mean.intensities.genes.csv')
}

{ ### Mutant-unique peaks ###
  peaks.unq.tab2 = t(sapply(strsplit(peaks.unq2, '-'), function(x) as.numeric(x)))
  rownames(peaks.unq.tab2) = peaks.unq2
  colnames(peaks.unq.tab2) = c('Chr', 'Start', 'End')

  system.time({
    FLANKING = 5000
    tmp.cols = NULL
    for ( i in 1:nrow(peaks.unq.tab2) ) {
      cat(i, '\n')
      
      y1 = rep(0, nrow(annot.unq))
      y2 = which(annot.unq$seqnames == peaks.unq.tab2[i, 'Chr'])
      y1[y2] = as.numeric(peaks.unq.tab2[i, 'Start'] >= annot.unq$start[y2] - FLANKING & 
                            peaks.unq.tab2[i, 'End'] <= annot.unq$end[y2] + FLANKING)

      if ( sum(y1) == 0 ) next
      
      tmp.cols = c(tmp.cols, peaks.unq2[i]) # peaks with gene annotations
      
      cat(y1, file = 'tmp.values', sep = '\n', append = T)
    }
    
    cat(tmp.cols, file = 'tmp.cols', sep = '\n') # For backup
  })

  system.time({
    tmp.values = scan('tmp.values')
  })

  system.time({
    peak.annot.mtx2 = matrix(tmp.values, nr = nrow(annot.unq), nc = length(tmp.cols), byrow = F)
    rownames(peak.annot.mtx2) = annot.unq$ID1
    colnames(peak.annot.mtx2) = tmp.cols
  })
  
  system.time({
    signalmap2.list = apply(peak.annot.mtx2, 2, function(x) annot.unq[x == 1, ])
  })

  (signal.genes.freq2 = sort(table(unlist(sapply(signalmap2.list, 
                                                 function(x) unique(x$gene_name)))), decreasing = T))[1:5]

  par(omi = c(0, .2, 0, 0), cex = 1.5)
  plot(sort(as.numeric(table(unlist(sapply(signalmap2.list, function(x) unique(x$gene_name))))), decreasing = T), 
       ylab = 'Number of overlapping peaks', xlab = 'Genes (ranked)', 
       main = 'Mutant-specific peaks', 
       log = 'xy')
  legend('topright', bty = 'n', cex = 1, 
         legend = paste0(1:5, '. ', names(signal.genes.freq2)[1:5], ': ', signal.genes.freq2[1:5]))
  
  signalmap2.genes = sapply(signalmap2.list, function(x) paste(unique(x$gene_name), collapse = '|'))

  system.time({
    save(peak.annot.mtx2, signalmap2.list, signalmap2.genes, 
         file = 'ANL001b.wt.mut.peaks.unique.mut.annot.rdata')
  })

  file.remove('tmp.values', 'tmp.cols')

  peaks.unq2.int.genes = cbind(peaks.unq.tab2, 
                               norm.means[rownames(peaks.unq.tab2), ], 
                               Genes = signalmap2.genes[rownames(peaks.unq.tab2)])

  write.csv(peaks.unq2.int.genes, quote = F, 
            file = 'ANL001b.wt.mut.peaks.unqiue.mut.allcells.mean.intensities.genes.csv')
}



###################################################################################################

##### Cluster-specific peaks : heatmaps #####

###################################################################################################

#(sobj.filt = readRDS('ANL001b.merged.data.processed.motif.rds'))

cells.by.cluster = split(rownames(sobj.filt@meta.data), sobj.filt@meta.data$seurat_clusters)

atac.norm = data.matrix(sobj.filt[['ATAC']]@data)

system.time({
  for ( i in 1:length(cells.by.cluster) ) {
    cat(i, '\n')
    
    CLUSTER = names(cells.by.cluster)[i]
    
    norm1 = data.matrix(atac.norm[, grep('WT', cells.by.cluster[[CLUSTER]], value = T)])
    norm2 = data.matrix(atac.norm[, grep('MUT', cells.by.cluster[[CLUSTER]], value = T)])
    
    norm.means = cbind(WT = rowMeans(norm1, na.rm = T), MUT = rowMeans(norm2, na.rm = T))
    
    peaks.each = apply(norm.means, 2, function(x) names(which(x > 0)))
    
    VennDiagram::venn.diagram(peaks.each, 
                              filename = paste0('ANL001b.merged.wt.mut.peaks.norm.means.cluster', 
                                                CLUSTER, '.venn.png'), 
                              main = paste0('Cluster ', CLUSTER), main.cex = 1.8, 
                              width = 600, height = 600, 
                              resolution = 100, compression = 'jpeg')
  }
})

system.time({
  KS.THRES = 0.01
  WC.THRES = 0.001
  
  results.all = NULL
  stats.all = NULL
  for ( i in 1:length(cells.by.cluster) ) {
    cat(i, '\n')
    
    CLUSTER = names(cells.by.cluster)[i]
    
    norm1 = atac.norm[, grep('WT', cells.by.cluster[[CLUSTER]], value = T)]
    norm2 = atac.norm[, grep('MUT', cells.by.cluster[[CLUSTER]], value = T)]
    
    tmp.stats = c(Cluster = as.integer(CLUSTER), 
                  N_Total = length(cells.by.cluster[[CLUSTER]]), 
                  N_WT = ncol(norm1), 
                  N_MUT = ncol(norm2))
    
    {
      norm1.means.each = rowMeans(norm1)
      norm2.means.each = rowMeans(norm2)
      
      norm1.medians.each = apply(norm1, 1, median)
      norm2.medians.each = apply(norm2, 1, median)
    }
    
    {
      tmp.test = cbind(Cluster = as.integer(CLUSTER), 
                       Mean_WT = norm1.means.each, Mean_MUT = norm2.means.each, 
                       Median_WT = norm1.medians.each, Median_MUT = norm2.medians.each)
      
      tmp.test = cbind(tmp.test, log2FC = apply(tmp.test, 1, function(x) log2(x['Mean_MUT'] / x['Mean_WT'])))
      
      pvals.ks = sapply(1:nrow(norm1), function(x) ks.test(norm1[x, ], norm2[x, ])$p.value)
      pvals.wc = sapply(1:nrow(norm1), function(x) wilcox.test(norm1[x, ], norm2[x, ])$p.value)
      
      names(pvals.ks) = rownames(norm1)
      names(pvals.wc) = rownames(norm1)
      
      results.all[[CLUSTER]]$all = cbind(tmp.test, P_KS = pvals.ks, P_WC = pvals.wc)
    }
    
    {
      tmp.signi = sort(unique(c(names(which(pvals.ks < KS.THRES)), names(which(pvals.wc < WC.THRES)))))
      results.all[[CLUSTER]]$signi = cbind(norm1[tmp.signi, ], norm2[tmp.signi, ])
    }
    
    {
      tmp.stats = c(tmp.stats, 
                    sapply(c(0.05, 0.01, 0.005, 0.001), function(x) sum(pvals.ks < x, na.rm = T)), 
                    sapply(c(0.05, 0.01, 0.005, 0.001), function(x) sum(pvals.wc < x, na.rm = T)))
      names(tmp.stats)[-c(1:4)] = c('KS_p05', 'KS_p01', 'KS_p005', 'KS_p001', 
                                    'WC_p05', 'WC_p01', 'WC_p005', 'WC_p001')
    }
    
    {
      norm.means = cbind(WT = rowMeans(norm1), MUT = rowMeans(norm2))
      
      peaks.each = apply(norm.means, 2, function(x) names(which(x > 0)))
      
      tmp.each = sapply(peaks.each, length)
      names(tmp.each) = c('N_peaks_WT', 'N_peaks_MUT')
      tmp.cmm = sum(apply(norm.means, 1, function(x) all(x > 0)))
    }
    
    stats.all = rbind(stats.all, c(tmp.stats, tmp.each, N_peaks_cmm = tmp.cmm))
  }

  rownames(stats.all) = NULL
  
  save(results.all, stats.all, file = 'ANL001b.merged.wt.mut.peaks.by.cluster.rdata')
  write.csv(stats.all, file = 'ANL001b.merged.wt.mut.peaks.by.cluster.stats.csv', row.names = F, quote = F)
})

{
  tmp.res = lapply(results.all, function(x) x$all[rownames(x$signi), ])
  
  res.signi = NULL
  for ( i in 1:length(tmp.res) ) {
    tmp1 = cbind(Cluster = tmp.res[[i]][, 'Cluster'], Peak = rownames(tmp.res[[i]]), tmp.res[[i]][, -1])
    res.signi = rbind(res.signi, tmp1)
  }
  
  rownames(res.signi) = NULL

  write.csv(res.signi, file = 'ANL001b.merged.wt.mut.peaks.by.cluster.signi.peaks.csv', row.names = F, quote = F)
}

system.time({
  results.corr = list()
  for ( i in 1:length(results.all) ) {
    cat(i, '\n')
    
    tmp.cluster = names(results.all)[i]

    norm.peaks.sig = results.all[[i]]$signi
    
    results.corr[[tmp.cluster]]$samples = cor(norm.peaks.sig, use = 'pairwise.complete.obs')
    results.corr[[tmp.cluster]]$features = cor(t(norm.peaks.sig), use = 'pairwise.complete.obs')
  }
  
  save(results.corr, file = 'ANL001b.merged.wt.mut.peaks.by.cluster.corr.rdata')
})

system.time({
  for ( i in 1:length(results.all) ) {
  #for ( i in c(35, 38, 39, 11) + 1 ) {
    cat(i, '\n')
    
    tmp.cluster = names(results.all)[i]
    tmp.res = results.all[[i]]$all
    
    tmp.p.ks = tmp.res[, 'P_KS']
    tmp.p.wc = tmp.res[, 'P_WC']
    
    tmp.p.ks = tmp.p.ks[!is.na(tmp.p.ks)]
    tmp.p.wc = tmp.p.wc[!is.na(tmp.p.wc)]
    
    { ### Histogram ###
      png(paste0('ANL001b.merged.wt.mut.peaks.norm.means.cluster', tmp.cluster, '.histo.png'), 
          width = 1600, height = 800)
      par(mfrow = c(1, 2), omi = c(0, .2, 0, 0), cex = 2)
      hist(-log10(tmp.p.ks), freq = F, col = 'white', border = 'blue', xlab = '-log10(p-value)', 
           main = paste0('Cluster ', tmp.cluster))
      hist(-log10(tmp.p.wc), freq = F, col = 'white', border = 'red', add = T)
      
      plot(density(-log10(tmp.p.ks)), col = 'blue', xlab = '-log10(p-value)', 
           main = paste0('Cluster ', tmp.cluster))
      lines(density(-log10(tmp.p.wc)), col = 'red')
      dev.off()
    }
    
    { ### Waterfall ###
      png(paste0('ANL001b.merged.wt.mut.peaks.norm.means.cluster', tmp.cluster, '.waterfall.png'), 
          width = 1000, height = 1000)
      par(omi = c(0, .2, 0, 0), cex = 2)
      plot(sort(tmp.p.ks), col = 'blue', log = 'xy', 
           ylim = range(c(tmp.p.ks, tmp.p.wc)), 
           xlab = 'Peaks (ordered)', 
           ylab = 'P value (log scale)', 
           main = paste0('Cluster ', tmp.cluster))
      points(sort(tmp.p.wc), col = 'red')
      abline(h = c(0.01), col = 'blue', lty = 2, lwd = 1)
      abline(h = c(0.001), col = 'red', lty = 2, lwd = 1)
      dev.off()
    }
    
    { ### Heatmap ###
      norm.peaks.sig = results.all[[i]]$signi
      
      png(paste0('ANL001b.merged.wt.mut.peaks.norm.means.cluster', tmp.cluster, '.heatmap.intensity.png'), 
          width = 1000, height = 1000)
      par(omi = c(1, 0, .1, 1), cex.main = 2)
      heatmap.2(norm.peaks.sig, trace = 'none', 
                keysize = .9, key.xlab = 'Peak intensity', key.par = list(cex.main = 1), 
                col = colorpanel(10, 'black', 'yellow'), 
                ColSideColors = sapply(colnames(norm.peaks.sig), 
                                       function(x) if ( substr(x, 1, 3) == 'WT_' ) 'blue' else 'red'), 
                main = paste0('Cluster ', tmp.cluster))
      dev.off()
      
      png(paste0('ANL001b.merged.wt.mut.peaks.norm.means.cluster', tmp.cluster, '.heatmap.corr.peaks.png'), 
          width = 1000, height = 1000)
      par(omi = c(1, 0, .1, 1), cex.main = 2)
      heatmap.2(results.corr[[tmp.cluster]]$features, trace = 'none', symm = T, 
                keysize = .9, key.xlab = 'Correlation', key.par = list(cex.main = 1), 
                col = colorpanel(10, 'orange', 'black', 'cyan'), 
                main = paste0('Cluster ', tmp.cluster))
      dev.off()
      
      tmp.corr.samp = results.corr[[tmp.cluster]]$samples
      tmp.sd = apply(tmp.corr.samp, 1, function(x) sd(x, na.rm = T))
      tmp.corr.samp = tmp.corr.samp[!is.na(tmp.sd), !is.na(tmp.sd)]
      
      png(paste0('ANL001b.merged.wt.mut.peaks.norm.means.cluster', tmp.cluster, '.heatmap.corr.cells.png'), 
          width = 1000, height = 1000)
      par(omi = c(1, 0, .1, 1), cex.main = 2)
      heatmap.2(tmp.corr.samp, trace = 'none', symm = T, 
                keysize = .9, key.xlab = 'Correlation', key.par = list(cex.main = 1), 
                col = colorpanel(10, 'orange', 'black', 'cyan'), 
                ColSideColors = sapply(colnames(tmp.corr.samp), 
                                       function(x) if ( substr(x, 1, 3) == 'WT_' ) 'blue' else 'red'), 
                RowSideColors = sapply(colnames(tmp.corr.samp), 
                                       function(x) if ( substr(x, 1, 3) == 'WT_' ) 'blue' else 'red'), 
                main = paste0('Cluster ', tmp.cluster))
      dev.off()
    }
  }
})



###################################################################################################

##### Cluster-specific peaks for WT and MUT each : heatmaps #####

###################################################################################################

KS.THRES = 0.01
WC.THRES = 0.001

#(sobj.filt = readRDS('ANL001b.merged.data.processed.motif.rds'))

cells.by.cluster = split(rownames(sobj.filt@meta.data), sobj.filt@meta.data$seurat_clusters)

atac.norm = data.matrix(sobj.filt[['ATAC']]@data)

{ ### Monte Carlo simulations ###
  system.time({
    N.SIMUL = 1000
    rnd.signi.counts.all = NULL
    rnd.signi.peaks.all = NULL
    for ( i in 1:length(cells.by.cluster) ) {
      CLUSTER = names(cells.by.cluster)[i]
      
      rnd.signi.counts = NULL
      rnd.signi.peaks = NULL
      for ( j in 1:N.SIMUL ) {
        cat(i, ':', CLUSTER, ':', j, '\n')
        
        { ### Random data ###
          data.rnd = atac.norm[, cells.by.cluster[[CLUSTER]]]
          
          sample1 = grep('WT', colnames(data.rnd))
          sample2 = grep('MUT', colnames(data.rnd))
          
          for ( k in 1:ncol(data.rnd) ) {
            tmp1a = sample(sample1, 1)
            tmp1b = sample(sample2, 1)
            
            tmp2a = data.rnd[, tmp1a]
            tmp2b = data.rnd[, tmp1b]
            
            data.rnd[, tmp1a] = tmp2b
            data.rnd[, tmp1b] = tmp2a
          }
          
          norm1 = data.rnd[, sample1]
          norm2 = data.rnd[, sample2]
        }
        
        {
          norm1.means.each = rowMeans(norm1)
          norm2.means.each = rowMeans(norm2)
          
          norm1.medians.each = apply(norm1, 1, median)
          norm2.medians.each = apply(norm2, 1, median)
        }
        
        {
          tmp.test = cbind(Cluster = as.integer(CLUSTER), 
                           Mean_WT = norm1.means.each, Mean_MUT = norm2.means.each, 
                           Median_WT = norm1.medians.each, Median_MUT = norm2.medians.each)
          
          tmp.test = cbind(tmp.test, 
                           Delta = apply(tmp.test, 1, function(x) x['Mean_MUT'] - x['Mean_WT']), 
                           log2FC = apply(tmp.test, 1, function(x) log2(x['Mean_MUT'] / x['Mean_WT'])))
          
          pvals.ks.mut = sapply(1:nrow(norm1), function(x) ks.test(norm1[x, ], norm2[x, ], alternative = 'g')$p.value)
          pvals.wc.mut = sapply(1:nrow(norm1), function(x) wilcox.test(norm1[x, ], norm2[x, ], alternative = 'l')$p.value)
          names(pvals.ks.mut) = rownames(norm1)
          names(pvals.wc.mut) = rownames(norm1)
          
          pvals.ks.wt = sapply(1:nrow(norm1), function(x) ks.test(norm1[x, ], norm2[x, ], alternative = 'l')$p.value)
          pvals.wc.wt = sapply(1:nrow(norm1), function(x) wilcox.test(norm1[x, ], norm2[x, ], alternative = 'g')$p.value)
          names(pvals.ks.wt) = rownames(norm1)
          names(pvals.wc.wt) = rownames(norm1)
          
          tmp.signi.mut = unique(c(names(which(pvals.ks.mut < KS.THRES)), names(which(pvals.wc.mut < WC.THRES))))
          tmp.signi.wt = unique(c(names(which(pvals.ks.wt < KS.THRES)), names(which(pvals.wc.wt < WC.THRES))))
          
          tmp.signi.counts = c(N_SigPeaks_MUT = length(tmp.signi.mut), N_SigPeaks_WT = length(tmp.signi.wt))
        }
        
        rnd.signi.counts = rbind(rnd.signi.counts, tmp.signi.counts)
        rnd.signi.peaks$M[[j]] = tmp.signi.mut
        rnd.signi.peaks$W[[j]] = tmp.signi.wt
      }
      
      rnd.signi.counts.all[[CLUSTER]] = rnd.signi.counts
      rnd.signi.peaks.all[[CLUSTER]] = rnd.signi.peaks
      
      write.csv(rnd.signi.counts, quote = F, row.names = F, 
                file = paste0('ANL001b.merged.wt.mut.peaks.by.cluster.by.cond.signi.simul.cluster', CLUSTER, '.csv'))
      save(rnd.signi.peaks, 
           file = paste0('ANL001b.merged.wt.mut.peaks.by.cluster.by.cond.signi.simul.cluster', CLUSTER, '.rdata'))
    }
    
    save(rnd.signi.counts.all, rnd.signi.peaks.all, 
         file = 'ANL001b.merged.wt.mut.peaks.by.cluster.by.cond.signi.simul.d220324.rdata')
  })
}

### Real data ###
system.time({
  results.all = NULL
  stats.all = NULL
  for ( i in 1:length(cells.by.cluster) ) {
    cat(i, '\n')
    
    CLUSTER = names(cells.by.cluster)[i]
    
    norm1 = atac.norm[, grep('WT', cells.by.cluster[[CLUSTER]], value = T)]
    norm2 = atac.norm[, grep('MUT', cells.by.cluster[[CLUSTER]], value = T)]
    
    tmp.stats = c(Cluster = as.integer(CLUSTER), 
                  N_cells_Total = length(cells.by.cluster[[CLUSTER]]), 
                  N_cells_WT = ncol(norm1), 
                  N_cells_MUT = ncol(norm2))
    
    {
      norm1.means.each = rowMeans(norm1)
      norm2.means.each = rowMeans(norm2)
      
      norm1.medians.each = apply(norm1, 1, median)
      norm2.medians.each = apply(norm2, 1, median)
    }
    
    {
      tmp.test = cbind(Cluster = as.integer(CLUSTER), 
                       Mean_WT = norm1.means.each, Mean_MUT = norm2.means.each, 
                       Median_WT = norm1.medians.each, Median_MUT = norm2.medians.each)
      
      tmp.test = cbind(tmp.test, 
                       Delta = apply(tmp.test, 1, function(x) x['Mean_MUT'] - x['Mean_WT']), 
                       log2FC = apply(tmp.test, 1, function(x) log2(x['Mean_MUT'] / x['Mean_WT'])))
      
      pvals.ks.mut = sapply(1:nrow(norm1), function(x) ks.test(norm1[x, ], norm2[x, ], alternative = 'g')$p.value)
      pvals.wc.mut = sapply(1:nrow(norm1), function(x) wilcox.test(norm1[x, ], norm2[x, ], alternative = 'l')$p.value)
      names(pvals.ks.mut) = rownames(norm1)
      names(pvals.wc.mut) = rownames(norm1)
      
      pvals.ks.wt = sapply(1:nrow(norm1), function(x) ks.test(norm1[x, ], norm2[x, ], alternative = 'l')$p.value)
      pvals.wc.wt = sapply(1:nrow(norm1), function(x) wilcox.test(norm1[x, ], norm2[x, ], alternative = 'g')$p.value)
      names(pvals.ks.wt) = rownames(norm1)
      names(pvals.wc.wt) = rownames(norm1)
      
      results.all[[CLUSTER]]$all = data.frame(tmp.test, 
                                              P_MUT_KS = pvals.ks.mut, P_MUT_WC = pvals.wc.mut, 
                                              P_WT_KS = pvals.ks.wt, P_WT_WC = pvals.wc.wt)
    }
    
    {
      tmp.signi.mut = unique(c(names(which(pvals.ks.mut < KS.THRES)), names(which(pvals.wc.mut < WC.THRES))))
      tmp.signi.mut = tmp.signi.mut[order(tmp.test[tmp.signi.mut, 'Delta'], decreasing = T)]
      tmp.signi.mut1 = data.frame(norm1)[tmp.signi.mut, ]
      tmp.signi.mut1 = tmp.signi.mut1[order(colMeans(tmp.signi.mut1), decreasing = T)]
      tmp.signi.mut2 = data.frame(norm2)[tmp.signi.mut, ]
      tmp.signi.mut2 = tmp.signi.mut2[order(colMeans(tmp.signi.mut2), decreasing = F)]
      results.all[[CLUSTER]]$signi_mut = cbind(tmp.signi.mut1, tmp.signi.mut2)

      tmp.signi.wt = unique(c(names(which(pvals.ks.wt < KS.THRES)), names(which(pvals.wc.wt < WC.THRES))))
      tmp.signi.wt = tmp.signi.wt[order(tmp.test[tmp.signi.wt, 'Delta'], decreasing = F)]
      tmp.signi.wt1 = data.frame(norm1)[tmp.signi.wt, ]
      tmp.signi.wt1 = tmp.signi.wt1[order(colMeans(tmp.signi.wt1), decreasing = T)]
      tmp.signi.wt2 = data.frame(norm2)[tmp.signi.wt, ]
      tmp.signi.wt2 = tmp.signi.wt2[order(colMeans(tmp.signi.wt2), decreasing = F)]
      results.all[[CLUSTER]]$signi_wt = cbind(tmp.signi.wt1, tmp.signi.wt2)
      
      tmp.stats = c(tmp.stats, N_SigPeaks_MUT = length(tmp.signi.mut), N_SigPeaks_WT = length(tmp.signi.wt))
    }
    
    {
      tmp.stats = c(tmp.stats, 
                    sapply(c(0.05, 0.01, 0.005, 0.001), function(x) sum(pvals.ks.mut < x, na.rm = T)), 
                    sapply(c(0.05, 0.01, 0.005, 0.001), function(x) sum(pvals.wc.mut < x, na.rm = T)), 
                    sapply(c(0.05, 0.01, 0.005, 0.001), function(x) sum(pvals.ks.wt < x, na.rm = T)), 
                    sapply(c(0.05, 0.01, 0.005, 0.001), function(x) sum(pvals.wc.wt < x, na.rm = T)))
      names(tmp.stats)[-c(1:6)] = c('MUT_KS_p05', 'MUT_KS_p01', 'MUT_KS_p005', 'MUT_KS_p001', 
                                    'MUT_WC_p05', 'MUT_WC_p01', 'MUT_WC_p005', 'MUT_WC_p001', 
                                    'WT_KS_p05', 'WT_KS_p01', 'WT_KS_p005', 'WT_KS_p001', 
                                    'WT_WC_p05', 'WT_WC_p01', 'WT_WC_p005', 'WT_WC_p001')
    }
    
    {
      norm.means = cbind(WT = rowMeans(norm1), MUT = rowMeans(norm2))
      
      peaks.each = apply(norm.means, 2, function(x) names(which(x > 0)))
      
      tmp.each = sapply(peaks.each, length)
      names(tmp.each) = c('N_peaks_WT', 'N_peaks_MUT')
      tmp.cmm = sum(apply(norm.means, 1, function(x) all(x > 0)))
    }
    
    stats.all = rbind(stats.all, c(tmp.stats, tmp.each, N_peaks_cmm = tmp.cmm))
  }
  
  rownames(stats.all) = NULL
  
  save(results.all, stats.all, file = 'ANL001b.merged.wt.mut.peaks.by.cluster.by.cond.d220301.rdata')
  write.csv(stats.all, file = 'ANL001b.merged.wt.mut.peaks.by.cluster.by.cond.stats.d220301.csv', 
            row.names = F, quote = F)
})

{ ### P-value for the number of significant peaks ###
  p.all = NULL
  for ( i in c('3', '8', '35', '38', '39', '42') ) {
    real.signi.cnt = stats.all[stats.all[, 'Cluster'] == i, c('N_SigPeaks_MUT', 'N_SigPeaks_WT')]
    
    rnd.signi.counts.all[[i]] = read.csv(paste0('ANL001b.merged.wt.mut.peaks.by.cluster.by.cond.signi.simul.cluster', i, '.csv'))
    rnd.signi.cnt = rnd.signi.counts.all[[i]]
    rnd.signi.cnt.mean = round(colMeans(rnd.signi.cnt), 2)
    rnd.signi.cnt.sd = round(apply(rnd.signi.cnt, 2, sd), 2)
    
    p.mut = sum(rnd.signi.cnt[, 'N_SigPeaks_MUT'] >= as.numeric(real.signi.cnt['N_SigPeaks_MUT'])) / nrow(rnd.signi.cnt)
    p.wt = sum(rnd.signi.cnt[, 'N_SigPeaks_WT'] >= as.numeric(real.signi.cnt['N_SigPeaks_WT'])) / nrow(rnd.signi.cnt)
    
    p.all = rbind(p.all, c(i, real.signi.cnt, rnd.signi.cnt.mean, rnd.signi.cnt.sd, p.mut, p.wt))
    
    {
      par(mfrow = c(1, 2), omi = c(0, .2, 0, 0))
      hist(rnd.signi.cnt[, 'N_SigPeaks_MUT'], 
           xlab = 'Number of significant peaks', main = paste0('Cluster', i, ': MUT'), 
           xlim = range(c(rnd.signi.cnt[, 'N_SigPeaks_MUT'], real.signi.cnt['N_SigPeaks_MUT'])))
      abline(v = real.signi.cnt['N_SigPeaks_MUT'], col = 'red')
      hist(rnd.signi.cnt[, 'N_SigPeaks_WT'], 
           xlab = 'Number of significant peaks', main = paste0('Cluster', i, ': WT'), 
           xlim = range(c(rnd.signi.cnt[, 'N_SigPeaks_WT'], real.signi.cnt['N_SigPeaks_WT'])))
      abline(v = real.signi.cnt['N_SigPeaks_WT'], col = 'red')
    }
  }
  
  colnames(p.all) = c('Cluster', 'SigPeaks_M', 'SigPeaks_W', 'Mean_M', 'Mean_W', 'SD_M', 'SD_W', 'p_M', 'p_W')
}

{
  tmp.res = lapply(results.all, function(x) x$all[c(rownames(x$signi_mut), rownames(x$signi_wt)), ])
  
  res.signi = NULL
  for ( i in 1:length(tmp.res) ) {
    tmp1 = cbind(Cluster = tmp.res[[i]][, 'Cluster'], Peak = rownames(tmp.res[[i]]), tmp.res[[i]][, -1])
    tmp1 = tmp1[order(tmp1$Delta, decreasing = T), ]
    res.signi = rbind(res.signi, tmp1)
  }
  
  rownames(res.signi) = NULL
  
  write.csv(res.signi, file = 'ANL001b.merged.wt.mut.peaks.by.cluster.by.cond.signi.peaks.d220301.csv', 
            row.names = F, quote = F)
}

system.time({
  results.corr = list()
  for ( i in 1:length(results.all) ) {
    cat(i, '\n')
    
    tmp.cluster = names(results.all)[i]
    norm.peaks.sig.mut = results.all[[i]]$signi_mut
    norm.peaks.sig.wt = results.all[[i]]$signi_wt
    
    results.corr[[tmp.cluster]]$samples_mut = cor(norm.peaks.sig.mut, use = 'pairwise.complete.obs')
    results.corr[[tmp.cluster]]$samples_wt = cor(norm.peaks.sig.wt, use = 'pairwise.complete.obs')
    
    results.corr[[tmp.cluster]]$features_mut = cor(t(norm.peaks.sig.mut), use = 'pairwise.complete.obs')
    results.corr[[tmp.cluster]]$features_wt = cor(t(norm.peaks.sig.wt), use = 'pairwise.complete.obs')
  }
  
  save(results.corr, file = 'ANL001b.merged.wt.mut.peaks.by.cluster.by.cond.corr.d220301.rdata')
})

system.time({
  for ( i in 1:length(results.all) ) {
    cat(i, '\n')
    
    tmp.cluster = names(results.all)[i]
    tmp.res = results.all[[i]]$all
    
    { ### Mutant ###
      tmp.p.ks.mut = tmp.res[, 'P_MUT_KS']
      tmp.p.wc.mut = tmp.res[, 'P_MUT_WC']
      
      tmp.p.ks.mut = tmp.p.ks.mut[!is.na(tmp.p.ks.mut)]
      tmp.p.wc.mut = tmp.p.wc.mut[!is.na(tmp.p.wc.mut)]
      
      { ### Histogram ###
        png(paste0('ANL001b.merged.wt.mut.peaks.norm.means.cluster', tmp.cluster, '.mut.histo.png'), 
            width = 1600, height = 800)
        par(mfrow = c(1, 2), omi = c(0, .2, 0, 0), cex = 2)
        hist(-log10(tmp.p.ks.mut), freq = F, col = 'white', border = 'blue', xlab = '-log10(p-value)', 
             main = paste0('Cluster ', tmp.cluster, ' : mutant'))
        hist(-log10(tmp.p.wc.mut), freq = F, col = 'white', border = 'red', add = T)
        
        plot(density(-log10(tmp.p.ks.mut)), col = 'blue', xlab = '-log10(p-value)', 
             main = paste0('Cluster ', tmp.cluster, ' : mutant'))
        lines(density(-log10(tmp.p.wc.mut)), col = 'red')
        dev.off()
      }
      
      { ### Waterfall ###
        png(paste0('ANL001b.merged.wt.mut.peaks.norm.means.cluster', tmp.cluster, '.mut.waterfall.png'), 
            width = 1000, height = 1000)
        par(omi = c(0, .2, 0, 0), cex = 2)
        plot(sort(tmp.p.ks.mut), col = 'blue', log = 'xy', 
             ylim = range(c(tmp.p.ks.mut, tmp.p.wc.mut)), 
             xlab = 'Peaks (ordered)', 
             ylab = 'P value (log scale)', 
             main = paste0('Cluster ', tmp.cluster, ' : mutant'))
        points(sort(tmp.p.wc.mut), col = 'red')
        abline(h = c(0.01), col = 'blue', lty = 2, lwd = 1)
        abline(h = c(0.001), col = 'red', lty = 2, lwd = 1)
        dev.off()
      }
    }
    
    { ### WT ###
      tmp.p.ks.wt = tmp.res[, 'P_WT_KS']
      tmp.p.wc.wt = tmp.res[, 'P_WT_WC']
      
      tmp.p.ks.wt = tmp.p.ks.wt[!is.na(tmp.p.ks.wt)]
      tmp.p.wc.wt = tmp.p.wc.wt[!is.na(tmp.p.wc.wt)]
      
      { ### Histogram ###
        png(paste0('ANL001b.merged.wt.mut.peaks.norm.means.cluster', tmp.cluster, '.wt.histo.png'), 
            width = 1600, height = 800)
        par(mfrow = c(1, 2), omi = c(0, .2, 0, 0), cex = 2)
        hist(-log10(tmp.p.ks.wt), freq = F, col = 'white', border = 'blue', xlab = '-log10(p-value)', 
             main = paste0('Cluster ', tmp.cluster, ' : WT'))
        hist(-log10(tmp.p.wc.wt), freq = F, col = 'white', border = 'red', add = T)
        
        plot(density(-log10(tmp.p.ks.wt)), col = 'blue', xlab = '-log10(p-value)', 
             main = paste0('Cluster ', tmp.cluster, ' : WT'))
        lines(density(-log10(tmp.p.wc.wt)), col = 'red')
        dev.off()
      }
      
      { ### Waterfall ###
        png(paste0('ANL001b.merged.wt.mut.peaks.norm.means.cluster', tmp.cluster, '.wt.waterfall.png'), 
            width = 1000, height = 1000)
        par(omi = c(0, .2, 0, 0), cex = 2)
        plot(sort(tmp.p.ks.wt), col = 'blue', log = 'xy', 
             ylim = range(c(tmp.p.ks.wt, tmp.p.wc.wt)), 
             xlab = 'Peaks (ordered)', 
             ylab = 'P value (log scale)', 
             main = paste0('Cluster ', tmp.cluster, ' : WT'))
        points(sort(tmp.p.wc.wt), col = 'red')
        abline(h = c(0.01), col = 'blue', lty = 2, lwd = 1)
        abline(h = c(0.001), col = 'red', lty = 2, lwd = 1)
        dev.off()
      }
    }
    
    { ### Heatmap: mutant ###
      norm.peaks.sig = data.matrix(results.all[[i]]$signi_mut)
      
      if ( nrow(norm.peaks.sig) >= 2 ) {
        png(paste0('ANL001b.merged.wt.mut.peaks.norm.means.cluster', tmp.cluster, '.mut.heatmap.intensity.png'), 
            width = 1000, height = 1000)
        par(omi = c(1, 0, .1, 1), cex.main = 2)
        heatmap.2(norm.peaks.sig, trace = 'none', 
                  keysize = .6, key.xlab = 'Peak intensity', key.par = list(cex.main = 1), 
                  col = colorpanel(10, 'black', 'yellow'), 
                  ColSideColors = sapply(colnames(norm.peaks.sig), 
                                         function(x) if ( substr(x, 1, 3) == 'WT_' ) 'blue' else 'red'), 
                  dendrogram = 'none', Colv = F, Rowv = F, 
                  main = paste0('Cluster ', tmp.cluster, ' : mutant'))
        dev.off()
        
        png(paste0('ANL001b.merged.wt.mut.peaks.norm.means.cluster', tmp.cluster, '.mut.heatmap.corr.peaks.png'), 
            width = 1000, height = 1000)
        par(omi = c(1, 0, .1, 1), cex.main = 2)
        heatmap.2(data.matrix(results.corr[[tmp.cluster]]$features_mut), trace = 'none', symm = T, 
                  keysize = .9, key.xlab = 'Correlation', key.par = list(cex.main = 1), 
                  col = colorpanel(10, 'orange', 'black', 'cyan'), 
                  main = paste0('Cluster ', tmp.cluster, ' : mutant'))
        dev.off()
        
        tmp.corr.samp = results.corr[[tmp.cluster]]$samples_mut
        tmp.sd = apply(tmp.corr.samp, 1, function(x) sd(x, na.rm = T))
        tmp.corr.samp = tmp.corr.samp[!is.na(tmp.sd), !is.na(tmp.sd)]
        
        png(paste0('ANL001b.merged.wt.mut.peaks.norm.means.cluster', tmp.cluster, '.mut.heatmap.corr.cells.png'), 
            width = 1000, height = 1000)
        par(omi = c(1, 0, .1, 1), cex.main = 2)
        heatmap.2(tmp.corr.samp, trace = 'none', symm = T, 
                  keysize = .9, key.xlab = 'Correlation', key.par = list(cex.main = 1), 
                  col = colorpanel(10, 'orange', 'black', 'cyan'), 
                  ColSideColors = sapply(colnames(tmp.corr.samp), 
                                         function(x) if ( substr(x, 1, 3) == 'WT_' ) 'blue' else 'red'), 
                  RowSideColors = sapply(colnames(tmp.corr.samp), 
                                         function(x) if ( substr(x, 1, 3) == 'WT_' ) 'blue' else 'red'), 
                  main = paste0('Cluster ', tmp.cluster, ' : mutant'))
        dev.off()
      }
    }
    
    { ### Heatmap: WT ###
      norm.peaks.sig = data.matrix(results.all[[i]]$signi_wt)
      
      if ( nrow(norm.peaks.sig) >= 2 ) {
        png(paste0('ANL001b.merged.wt.mut.peaks.norm.means.cluster', tmp.cluster, '.wt.heatmap.intensity.png'), 
            width = 1000, height = 1000)
        par(omi = c(1, 0, .1, 1), cex.main = 2)
        heatmap.2(norm.peaks.sig, trace = 'none', 
                  keysize = .6, key.xlab = 'Peak intensity', key.par = list(cex.main = 1), 
                  col = colorpanel(10, 'black', 'yellow'), 
                  ColSideColors = sapply(colnames(norm.peaks.sig), 
                                         function(x) if ( substr(x, 1, 3) == 'WT_' ) 'blue' else 'red'), 
                  dendrogram = 'none', Colv = F, Rowv = F, 
                  main = paste0('Cluster ', tmp.cluster, ' : WT'))
        dev.off()
        
        png(paste0('ANL001b.merged.wt.mut.peaks.norm.means.cluster', tmp.cluster, '.wt.heatmap.corr.peaks.png'), 
            width = 1000, height = 1000)
        par(omi = c(1, 0, .1, 1), cex.main = 2)
        heatmap.2(data.matrix(results.corr[[tmp.cluster]]$features_wt), trace = 'none', symm = T, 
                  keysize = .9, key.xlab = 'Correlation', key.par = list(cex.main = 1), 
                  col = colorpanel(10, 'orange', 'black', 'cyan'), 
                  main = paste0('Cluster ', tmp.cluster, ' : WT'))
        dev.off()
        
        tmp.corr.samp = results.corr[[tmp.cluster]]$samples_wt
        tmp.sd = apply(tmp.corr.samp, 1, function(x) sd(x, na.rm = T))
        tmp.corr.samp = tmp.corr.samp[!is.na(tmp.sd), !is.na(tmp.sd)]
        
        png(paste0('ANL001b.merged.wt.mut.peaks.norm.means.cluster', tmp.cluster, '.wt.heatmap.corr.cells.png'), 
            width = 1000, height = 1000)
        par(omi = c(1, 0, .1, 1), cex.main = 2)
        heatmap.2(tmp.corr.samp, trace = 'none', symm = T, 
                  keysize = .9, key.xlab = 'Correlation', key.par = list(cex.main = 1), 
                  col = colorpanel(10, 'orange', 'black', 'cyan'), 
                  ColSideColors = sapply(colnames(tmp.corr.samp), 
                                         function(x) if ( substr(x, 1, 3) == 'WT_' ) 'blue' else 'red'), 
                  RowSideColors = sapply(colnames(tmp.corr.samp), 
                                         function(x) if ( substr(x, 1, 3) == 'WT_' ) 'blue' else 'red'), 
                  main = paste0('Cluster ', tmp.cluster, ' : WT'))
        dev.off()
      }
    }
    
    { ### Heatmap: mutant + WT ###
      norm.peaks.sig = data.matrix(rbind(results.all[[i]]$signi_mut, results.all[[i]]$signi_wt))
      
      png(paste0('ANL001b.merged.wt.mut.peaks.norm.means.cluster', tmp.cluster, '.mut.wt.heatmap.intensity.png'), 
          width = 1000, height = 1000)
      par(omi = c(1, 0, .1, 1), cex.main = 2)
      heatmap.2(norm.peaks.sig, trace = 'none', 
                keysize = .6, key.xlab = 'Peak intensity', key.par = list(cex.main = 1), 
                col = colorpanel(10, 'black', 'yellow'), 
                ColSideColors = sapply(colnames(norm.peaks.sig), 
                                       function(x) if ( substr(x, 1, 3) == 'WT_' ) 'blue' else 'red'), 
                dendrogram = 'none', Colv = F, Rowv = F, 
                main = paste0('Cluster ', tmp.cluster, ' : mutant + WT'))
      dev.off()
      
      tmp.corr.feat = cor(t(norm.peaks.sig), use = 'pairwise.complete.obs')
      tmp.sd = apply(tmp.corr.feat, 1, function(x) sd(x, na.rm = T))
      tmp.corr.feat = tmp.corr.feat[!is.na(tmp.sd), !is.na(tmp.sd)]
      
      png(paste0('ANL001b.merged.wt.mut.peaks.norm.means.cluster', tmp.cluster, '.mut.wt.heatmap.corr.peaks.png'), 
          width = 1000, height = 1000)
      par(omi = c(1, 0, .1, 1), cex.main = 2)
      heatmap.2(tmp.corr.feat, trace = 'none', symm = T, 
                keysize = .9, key.xlab = 'Correlation', key.par = list(cex.main = 1), 
                col = colorpanel(10, 'orange', 'black', 'cyan'), 
                main = paste0('Cluster ', tmp.cluster, ' : mutant + WT'))
      dev.off()
      
      tmp.corr.samp = cor(norm.peaks.sig, use = 'pairwise.complete.obs')
      tmp.sd = apply(tmp.corr.samp, 1, function(x) sd(x, na.rm = T))
      tmp.corr.samp = tmp.corr.samp[!is.na(tmp.sd), !is.na(tmp.sd)]
      
      png(paste0('ANL001b.merged.wt.mut.peaks.norm.means.cluster', tmp.cluster, '.mut.wt.heatmap.corr.cells.png'), 
          width = 1000, height = 1000)
      par(omi = c(1, 0, .1, 1), cex.main = 2)
      heatmap.2(tmp.corr.samp, trace = 'none', symm = T, 
                keysize = .9, key.xlab = 'Correlation', key.par = list(cex.main = 1), 
                col = colorpanel(10, 'orange', 'black', 'cyan'), 
                ColSideColors = sapply(colnames(tmp.corr.samp), 
                                       function(x) if ( substr(x, 1, 3) == 'WT_' ) 'blue' else 'red'), 
                RowSideColors = sapply(colnames(tmp.corr.samp), 
                                       function(x) if ( substr(x, 1, 3) == 'WT_' ) 'blue' else 'red'), 
                main = paste0('Cluster ', tmp.cluster, ' : mutant + WT'))
      dev.off()
    }
  }
})



###################################################################################################

##### Significant peaks by cluster : gene annotation within +/-5k #####

###################################################################################################

#(sobj.filt = readRDS('ANL001b.merged.data.processed.motif.rds'))

#res.signi = read.csv('ANL001b.merged.wt.mut.peaks.by.cluster.signi.peaks.csv')
peaks.signi = unique(res.signi$Peak)

{
  annot = data.frame(Annotation(sobj.filt))
  annot$ID1 = apply(annot[, 1:3], 1, function(x) gsub(' +', '', paste(x, collapse = '-')))
  
  annot.unq = annot[!duplicated(annot$ID1), ]
}

{
  peaks.unq.tab1 = t(sapply(strsplit(peaks.signi, '-'), function(x) as.numeric(x)))
  rownames(peaks.unq.tab1) = peaks.signi
  colnames(peaks.unq.tab1) = c('Chr', 'Start', 'End')
  
  system.time({
    FLANKING = 5000
    tmp.cols = NULL
    for ( i in 1:nrow(peaks.unq.tab1) ) {
      cat(i, '\n')
      
      y1 = rep(0, nrow(annot.unq))
      y2 = which(annot.unq$seqnames == peaks.unq.tab1[i, 'Chr'])
      y1[y2] = as.numeric(peaks.unq.tab1[i, 'Start'] >= annot.unq$start[y2] - FLANKING & 
                            peaks.unq.tab1[i, 'End'] <= annot.unq$end[y2] + FLANKING)
      
      if ( sum(y1) == 0 ) next
      
      tmp.cols = c(tmp.cols, peaks.signi[i]) # peaks with gene annotations
      
      cat(y1, file = 'tmp.values', sep = '\n', append = T)
    }
    
    cat(tmp.cols, file = 'tmp.cols', sep = '\n') # For backup
  })

  system.time({
    tmp.values = scan('tmp.values')
  })

  system.time({
    peak.annot.mtx1 = matrix(tmp.values, nr = nrow(annot.unq), nc = length(tmp.cols), byrow = F)
    rownames(peak.annot.mtx1) = annot.unq$ID1
    colnames(peak.annot.mtx1) = tmp.cols
  })

  save(peak.annot.mtx1, annot.unq, peaks.unq.tab1, file = 'tmp.rdata')

  rm(tmp.values)
  
  system.time({
    signalmap1.list = apply(peak.annot.mtx1, 2, function(x) annot.unq[x == 1, ])
  })

  signal.genes.freq1 = sort(table(unlist(sapply(signalmap1.list, 
                                                function(x) unique(x$gene_name)))), decreasing = T)
  
  signalmap1.genes = sapply(signalmap1.list, function(x) paste(unique(x$gene_name), collapse = '|'))
  
  system.time({
    save(peak.annot.mtx1, signalmap1.list, signalmap1.genes, 
         file = 'ANL001b.merged.wt.mut.peaks.by.cluster.signi.annot.rdata')
  })

  file.remove('tmp.values', 'tmp.cols', 'tmp.rdata')

  peaks.signi.int.genes = cbind(peaks.unq.tab1, 
                                Genes = signalmap1.genes[rownames(peaks.unq.tab1)])
}

{
  res.signi = cbind(res.signi, peaks.signi.int.genes[res.signi$Peak, ])

  write.csv(res.signi, quote = F, row.names = F, 
            file = 'ANL001b.merged.wt.mut.peaks.by.cluster.signi.peaks.csv')
}



###################################################################################################

##### Significant peaks by cluster and by condition: gene annotation within +/-5k #####

###################################################################################################

#(sobj.filt = readRDS('ANL001b.merged.data.processed.motif.rds'))

#res.signi = read.csv('ANL001b.merged.wt.mut.peaks.by.cluster.by.cond.signi.peaks.d220301.csv')
peaks.signi = unique(res.signi$Peak)

{
  annot = data.frame(Annotation(sobj.filt))
  annot$ID1 = apply(annot[, 1:3], 1, function(x) gsub(' +', '', paste(x, collapse = '-')))
  
  annot.unq = annot[!duplicated(annot$ID1), ]
}

{
  peaks.unq.tab1 = t(sapply(strsplit(peaks.signi, '-'), function(x) as.numeric(x)))
  rownames(peaks.unq.tab1) = peaks.signi
  colnames(peaks.unq.tab1) = c('Chr', 'Start', 'End')
  
  system.time({
    FLANKING = 5000
    tmp.cols = NULL
    for ( i in 1:nrow(peaks.unq.tab1) ) {
      cat(i, '\n')
      
      y1 = rep(0, nrow(annot.unq))
      y2 = which(annot.unq$seqnames == peaks.unq.tab1[i, 'Chr'])
      y1[y2] = as.numeric(peaks.unq.tab1[i, 'Start'] >= annot.unq$start[y2] - FLANKING & 
                            peaks.unq.tab1[i, 'End'] <= annot.unq$end[y2] + FLANKING)
      
      if ( sum(y1) == 0 ) next
      
      tmp.cols = c(tmp.cols, peaks.signi[i]) # peaks with gene annotations
      
      cat(y1, file = 'tmp.values', sep = '\n', append = T)
    }
    
    cat(tmp.cols, file = 'tmp.cols', sep = '\n') # For backup
  })

  system.time({
    tmp.values = scan('tmp.values')
  })

  system.time({
    peak.annot.mtx1 = matrix(tmp.values, nr = nrow(annot.unq), nc = length(tmp.cols), byrow = F)
    rownames(peak.annot.mtx1) = annot.unq$ID1
    colnames(peak.annot.mtx1) = tmp.cols
  })

  save(peak.annot.mtx1, annot.unq, peaks.unq.tab1, file = 'tmp.rdata')

  rm(tmp.values)
  
  system.time({
    signalmap1.list = apply(peak.annot.mtx1, 2, function(x) annot.unq[x == 1, ])
  })

  signal.genes.freq1 = sort(table(unlist(sapply(signalmap1.list, 
                                                function(x) unique(x$gene_name)))), decreasing = T)
  
  signalmap1.genes = sapply(signalmap1.list, function(x) paste(unique(x$gene_name), collapse = '|'))
  
  system.time({
    save(peak.annot.mtx1, signalmap1.list, signalmap1.genes, 
         file = 'ANL001b.merged.wt.mut.peaks.by.cluster.by.cond.signi.annot.d220309.rdata')
  })

  file.remove('tmp.values', 'tmp.cols', 'tmp.rdata')

  peaks.signi.int.genes = cbind(peaks.unq.tab1, 
                                Genes = signalmap1.genes[rownames(peaks.unq.tab1)])
}

{
  res.signi = cbind(res.signi, peaks.signi.int.genes[res.signi$Peak, ])

  write.csv(res.signi, quote = F, row.names = F, 
            file = 'ANL001b.merged.wt.mut.peaks.by.cluster.by.cond.signi.peaks.d220309.csv')
}



###################################################################################################

##### Significant peaks by cluster and by condition: enriched motif analysis #####

###################################################################################################

#(sobj.filt = readRDS('ANL001b.merged.data.processed.motif.rds'))

DefaultAssay(sobj.filt) = 'ATAC'

#res.signi = read.csv('ANL001b.merged.wt.mut.peaks.by.cluster.by.cond.signi.peaks.d220309.csv')

{
  (gr.wt = GRanges(sub('-', ':', res.signi$Peak[res.signi$Delta < 0])))
  (gr.mut = GRanges(sub('-', ':', res.signi$Peak[res.signi$Delta > 0])))
}

{
  (atac.peaks.wt = CreateChromatinAssay(
    data = sobj.filt[['ATAC']]@data[res.signi$Peak[res.signi$Delta < 0], ], 
    ranges = gr.wt, 
    genome = 'danRer11'
  ))

  (atac.peaks.mut = CreateChromatinAssay(
    data = sobj.filt[['ATAC']]@data[res.signi$Peak[res.signi$Delta > 0], ], 
    ranges = gr.mut, 
    genome = 'danRer11'
  ))

  atac.peaks.wt@counts = sobj.filt[['ATAC']]@counts[res.signi$Peak[res.signi$Delta < 0], ]
  atac.peaks.mut@counts = sobj.filt[['ATAC']]@counts[res.signi$Peak[res.signi$Delta > 0], ]
}

{
  levels(gr.wt@seqnames) = paste0('chr', levels(gr.wt@seqnames)) # to match to BSgenome.Drerio.UCSC.danRer11
  seqlevels(gr.wt) = levels(gr.wt@seqnames)

  levels(gr.mut@seqnames) = paste0('chr', levels(gr.mut@seqnames)) # to match to BSgenome.Drerio.UCSC.danRer11
  seqlevels(gr.mut) = levels(gr.mut@seqnames)
}

# Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object
(pfm.set = getMatrixSet(x = JASPAR2020, opts = list(all = T)))

{
  system.time({
    motif.matrix.wt = CreateMotifMatrix(features = gr.wt, pwm = pfm.set, 
                                        genome = BSgenome.Drerio.UCSC.danRer11)
  })

  system.time({
    motif.matrix.mut = CreateMotifMatrix(features = gr.mut, pwm = pfm.set, 
                                         genome = BSgenome.Drerio.UCSC.danRer11)
  })
}

{
  rownames(motif.matrix.wt) = sub('^chr', '', rownames(motif.matrix.wt))
  rownames(motif.matrix.mut) = sub('^chr', '', rownames(motif.matrix.mut))
  
  # Create a new Mofif object to store the results
  (motif.object.wt = CreateMotifObject(data = motif.matrix.wt, pwm = pfm.set))
  (motif.object.mut = CreateMotifObject(data = motif.matrix.mut, pwm = pfm.set))

  atac.peaks.wt = SetAssayData(atac.peaks.wt, slot = 'motifs', new.data = motif.object.wt)
  atac.peaks.mut = SetAssayData(atac.peaks.mut, slot = 'motifs', new.data = motif.object.mut)
}

{
  sobj.filt.wt = sobj.filt
  sobj.filt.wt[['ATAC']] = atac.peaks.wt
  
  sobj.filt.mut = sobj.filt
  sobj.filt.mut[['ATAC']] = atac.peaks.mut
  
  sobj.filt.wt[['ACTIVITY']] = NULL
  sobj.filt.mut[['ACTIVITY']] = NULL
  
  sobj.filt.wt[['chromvar']] = NULL
  sobj.filt.mut[['chromvar']] = NULL
}

{ ### chromVar: Note that this step takes a long time ###
  system.time({
    sobj.filt.wt[['ATAC']]@ranges = gr.wt=
    
    register(MulticoreParam(10))
    
    sobj.filt.wt = RunChromVAR(
      object = sobj.filt.wt, # use @count
      motif.matrix = motif.matrix.wt, 
      genome = BSgenome.Drerio.UCSC.danRer11
    )

    saveRDS(sobj.filt.wt, 'ANL001b.merged.wt.mut.peaks.by.cluster.by.cond.signi.peaks.wt.motif.rds')
  })

  system.time({
    sobj.filt.mut[['ATAC']]@ranges = gr.mut
    
    register(MulticoreParam(10))
    
    sobj.filt.mut = RunChromVAR(
      object = sobj.filt.mut, # use @count
      motif.matrix = motif.matrix.mut, 
      genome = BSgenome.Drerio.UCSC.danRer11
    )

    saveRDS(sobj.filt.mut, 'ANL001b.merged.wt.mut.peaks.by.cluster.by.cond.signi.peaks.mut.motif.rds')
  })
}

{ ### presto: WT ###
  markers.rna.wt = presto:::wilcoxauc(X = sobj.filt.wt, group_by = 'seurat_clusters', 
                                      assay = 'data', seurat_assay = 'RNA')
  markers.motifs.wt = presto:::wilcoxauc(X = sobj.filt.wt, group_by = 'seurat_clusters', 
                                         assay = 'data', seurat_assay = 'chromvar')
  
  colnames(markers.rna.wt) = paste0("RNA.", colnames(markers.rna.wt))
  markers.rna.wt$gene = markers.rna.wt$RNA.feature
  
  motif.names.wt = markers.motifs.wt$feature
  colnames(markers.motifs.wt) = paste0("motif.", colnames(markers.motifs.wt))
  markers.motifs.wt$gene = ConvertMotifID(sobj.filt.wt, id = motif.names.wt)
}

{ ### presto: Mut ###
  markers.rna.mut = presto:::wilcoxauc(X = sobj.filt.mut, group_by = 'seurat_clusters', 
                                       assay = 'data', seurat_assay = 'RNA')
  markers.motifs.mut = presto:::wilcoxauc(X = sobj.filt.mut, group_by = 'seurat_clusters', 
                                          assay = 'data', seurat_assay = 'chromvar')
  
  colnames(markers.rna.mut) = paste0("RNA.", colnames(markers.rna.mut))
  markers.rna.mut$gene = markers.rna.mut$RNA.feature
  
  motif.names.mut = markers.motifs.mut$feature
  colnames(markers.motifs.mut) = paste0("motif.", colnames(markers.motifs.mut))
  markers.motifs.mut$gene = ConvertMotifID(sobj.filt.mut, id = motif.names.mut)
}

topTFs = function(markers_rna, markers_motifs, cluster, padj.cutoff = 0.2) {
  ctmarkers_rna = dplyr::filter(
    markers_rna, RNA.group == cluster, RNA.padj < padj.cutoff, RNA.logFC > 0) %>% 
    arrange(-RNA.auc)
  
  ctmarkers_motif = dplyr::filter(
    markers_motifs, motif.group == cluster, motif.padj < padj.cutoff, motif.logFC > 0) %>% 
    arrange(-motif.auc)
  ctmarkers_motif$gene = tolower(ctmarkers_motif$gene)
  
  top_tfs = inner_join(
    x = ctmarkers_rna[, c('RNA.group', 'gene', 'RNA.avgExpr', 'RNA.logFC', 
                          'RNA.auc', 'RNA.pval', 'RNA.padj', 'RNA.pct_in', 'RNA.pct_out')], 
    y = ctmarkers_motif[, c('motif.feature', 'gene', 'motif.avgExpr', 'motif.logFC', 
                            'motif.auc', 'motif.pval', 'motif.padj', 'motif.pct_in', 'motif.pct_out')], 
    by = "gene"
  )
  
  colnames(top_tfs)[1] = 'cluster'
  
  top_tfs$avg_auc = (top_tfs$RNA.auc + top_tfs$motif.auc) / 2
  top_tfs = arrange(top_tfs, -avg_auc)
  
  return(top_tfs)
}

{ ### identify top markers in all clusters ###
  top.tfs.all.wt = NULL
  for ( i in sort(unique(sobj.filt.wt@meta.data$seurat_clusters)) ) {
    top.tfs.all.wt = rbind(top.tfs.all.wt, topTFs(markers.rna.wt, markers.motifs.wt, 
                                                  cluster = i, padj.cutoff = 0.2))
  }
  
  top.tfs.all.mut = NULL
  for ( i in sort(unique(sobj.filt.mut@meta.data$seurat_clusters)) ) {
    top.tfs.all.mut = rbind(top.tfs.all.mut, topTFs(markers.rna.mut, markers.motifs.mut, 
                                                    cluster = i, padj.cutoff = 0.2))
  }
  
  write.csv(top.tfs.all.wt, file = 'ANL001b.merged.wt.mut.peaks.by.cluster.by.cond.signi.peaks.wt.motif.top.csv', 
            quote = F, row.names = F)
  write.csv(top.tfs.all.mut, file = 'ANL001b.merged.wt.mut.peaks.by.cluster.by.cond.signi.peaks.mut.motif.top.csv', 
            quote = F, row.names = F)
}

{
  top.tfs.unq.wt = top.tfs.all.wt[!top.tfs.all.wt$gene %in% top.tfs.all.mut$gene, ]
  top.tfs.unq.mut = top.tfs.all.mut[!top.tfs.all.mut$gene %in% top.tfs.all.wt$gene, ]
  
  top.tfs.unq = rbind(cbind(Cond = 'WT', top.tfs.unq.wt), cbind(Cond = 'Mut', top.tfs.unq.mut))

  write.csv(top.tfs.unq, file = 'ANL001b.merged.wt.mut.peaks.by.cluster.by.cond.signi.peaks.motif.top.unq.csv', 
            quote = F, row.names = F)
}

{
  motifs.top.peaks.wt = motif.matrix.wt[, top.tfs.unq.wt$motif.feature]
  motifs.top.peaks.wt = motifs.top.peaks.wt[rowSums(motifs.top.peaks.wt) > 0, ]

  motifs.top.peaks.mut = motif.matrix.mut[, top.tfs.unq.mut$motif.feature]
  motifs.top.peaks.mut = motifs.top.peaks.mut[rowSums(motifs.top.peaks.mut) > 0, ]
}


##### WNN-UMAP plots for top motifs/TFs in each cluster #####
if ( !dir.exists('signi.peaks.by.cluster.by.cond') ) dir.create('signi.peaks.by.cluster.by.cond')
if ( !dir.exists('signi.peaks.by.cluster.by.cond/wnn.umap.top.motifs.tfs') ) 
  dir.create('signi.peaks.by.cluster.by.cond/wnn.umap.top.motifs.tfs')

{ ### Multiple png files ###
  ftn.wnn.umap.toptfs = function(FILENAME, sobj, cond) {
    top.tfs.all = read.csv(FILENAME)
    all.clusters = unique(top.tfs.all$cluster)
    for ( tmp.cluster in all.clusters ) {
      cat(tmp.cluster, '\n')
      top.tfs.sub = top.tfs.all[top.tfs.all$cluster == tmp.cluster, ]
      
      tmp.title = paste0(toupper(cond), ': Cluster ', tmp.cluster)
      for ( i in 1:nrow(top.tfs.sub) ) {
        tmp.gene = top.tfs.sub$gene[i]
        tmp.motif = top.tfs.sub$motif.feature[i]
        
        if ( tmp.gene %in% rownames(sobj[['SCT']]) ) {
          DefaultAssay(sobj) = 'SCT'
          gene_plot = FeaturePlot(sobj, features = paste0('sct_', tmp.gene), 
                                  reduction = 'wnn.umap', pt.size = .1) + ggtitle(tmp.gene)
        } else {
          DefaultAssay(sobj) = 'RNA'
          gene_plot = FeaturePlot(sobj, features = tmp.gene, 
                                  reduction = 'wnn.umap', pt.size = .1) + ggtitle(tmp.gene)
        }
        
        DefaultAssay(sobj) = 'chromvar'
        motif_plot = FeaturePlot(sobj, features = tmp.motif, min.cutoff = 0, pt.size = .1, 
                                 cols = c("lightgrey", "darkred"), reduction = 'wnn.umap')
        
        p = LabelClusters(gene_plot, id = 'ident', clusters = as.character(tmp.cluster), size = 10, repel = F) + 
          LabelClusters(motif_plot, id = 'ident', clusters = as.character(tmp.cluster), size = 10, repel = F) + 
          plot_layout(ncol = 2) + 
          plot_annotation(title = tmp.title, 
                          theme = theme(plot.title = element_text(size = 25, hjust = .5, vjust = .5)))
        
        ggsave(plot = p, width = 12, height = 6, 
               paste0('signi.peaks.by.cluster.by.cond/wnn.umap.top.motifs.tfs/wnnumap.merged.', cond, 
                      '.cluster', tmp.cluster, '.', tmp.gene, '.', tmp.motif, '.png'))
      }
    }
  }
  
  system.time({
    ftn.wnn.umap.toptfs('ANL001b.merged.wt.mut.peaks.by.cluster.by.cond.signi.peaks.wt.motif.top.csv', 
                        sobj.filt.wt, filename.anno = 'wt')
  })

  system.time({
    ftn.wnn.umap.toptfs('ANL001b.merged.wt.mut.peaks.by.cluster.by.cond.signi.peaks.mut.motif.top.csv', 
                        sobj.filt.mut, filename.anno = 'mut')
  })
}

{ ### A single pdf file ###
  ftn.wnn.umap.toptfs.plots.list = function(FILENAME, sobj, cond) {
    top.tfs.all = read.csv(FILENAME)
    all.clusters = unique(top.tfs.all$cluster)
    plots.list = list()
    for ( tmp.cluster in all.clusters ) {
      cat(tmp.cluster, '\n')
      (top.tfs.sub = top.tfs.all[top.tfs.all$cluster == tmp.cluster, ])
      
      tmp.title = paste0(toupper(cond), ': Cluster ', tmp.cluster)
      tmp.list = list()
      for ( i in 1:nrow(top.tfs.sub) ) {
        tmp.gene = top.tfs.sub$gene[i]
        tmp.motif = top.tfs.sub$motif.feature[i]
        
        if ( tmp.gene %in% rownames(sobj[['SCT']]) ) {
          DefaultAssay(sobj) = 'SCT'
          gene_plot = FeaturePlot(sobj, features = paste0('sct_', tmp.gene), 
                                  reduction = 'wnn.umap', pt.size = .1) + ggtitle(tmp.gene)
        } else {
          DefaultAssay(sobj) = 'RNA'
          gene_plot = FeaturePlot(sobj, features = tmp.gene, 
                                  reduction = 'wnn.umap', pt.size = .1) + ggtitle(tmp.gene)
        }
        
        DefaultAssay(sobj) = 'chromvar'
        motif_plot = FeaturePlot(sobj, features = tmp.motif, min.cutoff = 0, pt.size = .1, 
                                 cols = c("lightgrey", "darkred"), reduction = 'wnn.umap')
        
        tmp.p = 
          LabelClusters(gene_plot, id = 'ident', clusters = as.character(tmp.cluster), size = 10, repel = F) + 
          LabelClusters(motif_plot, id = 'ident', clusters = as.character(tmp.cluster), size = 10, repel = F) + 
          plot_layout(ncol = 2) + 
          plot_annotation(title = tmp.title, 
                          theme = theme(plot.title = element_text(size = 25, hjust = .5, vjust = .5)))
        
        tmp.list[[i]] = tmp.p
      }
      
      plots.list = c(plots.list, tmp.list)
    }
    
    return(plots.list)
  }
  
  system.time({
    pdf('signi.peaks.by.cluster.by.cond/wnn.umap.top.motifs.tfs/umap.wnn.top.motifs.tfs.wt.all.pdf', 
        width = 12, height = 6)
    invisible(lapply(
      ftn.wnn.umap.toptfs.plots.list('ANL001b.merged.wt.mut.peaks.by.cluster.by.cond.signi.peaks.wt.motif.top.csv', 
                                     sobj.filt.wt, cond = 'wt'), print))
    dev.off()
  })

  system.time({
    pdf('signi.peaks.by.cluster.by.cond/wnn.umap.top.motifs.tfs/umap.wnn.top.motifs.tfs.mut.all.pdf', 
        width = 12, height = 6)
    invisible(lapply(
      ftn.wnn.umap.toptfs.plots.list('ANL001b.merged.wt.mut.peaks.by.cluster.by.cond.signi.peaks.mut.motif.top.csv', 
                                     sobj.filt.mut, cond = 'mut'), print))
    dev.off()
  })
}


### ATAC peak coverage plots ###
system.time({
  sobj.sub.list = list()
  for ( i in levels(sobj.filt) ) {
    cat(i, '\n')
    
    cells.sub = WhichCells(sobj.filt, idents = i)
    
    (sobj.sub = CreateSeuratObject(counts = sobj.filt[['RNA']]@counts[, cells.sub]))
    (sobj.sub[['SCT']] = CreateSCTAssayObject(counts = sobj.filt[['SCT']]@counts[, cells.sub], 
                                              scale.data = sobj.filt[['SCT']]@scale.data[, cells.sub], 
                                              SCTModel.list = sobj.filt[['SCT']]@SCTModel.list))
    (sobj.sub[['ATAC']] = CreateChromatinAssay(data = sobj.filt[['ATAC']]@data[, cells.sub], 
                                               ranges = sobj.filt[['ATAC']]@ranges, 
                                               annotation = sobj.filt[['ATAC']]@annotation, 
                                               fragments = sobj.filt[['ATAC']]@fragments))
    sobj.sub[['ATAC']]@seqinfo = sobj.filt[['ATAC']]@seqinfo
    sobj.sub[['ATAC']]@meta.features = sobj.filt[['ATAC']]@meta.features
    sobj.sub[['ATAC']]@var.features = sobj.filt[['ATAC']]@var.features
    sobj.sub[['ATAC']]@motifs = sobj.filt[['ATAC']]@motifs
    
    sobj.sub.list[[i]] = sobj.sub
  }
})

all.genes.wt = sort(unique(top.tfs.all.wt[, 'gene']))

### List of individual plots ###
system.time({
  plots.list = list()
  for ( i in 1:length(all.genes) ) {
    cat(i, '\n')
    
    if ( i == 98 ) next
    
    g = all.genes[i]
    tmp.anno = Annotation(sobj.filt[['ATAC']])
    tmp.anno = tmp.anno[tmp.anno$gene_name == g]
    
    p = CoveragePlot(sobj.filt, region = g, features = g, peaks = F, 
                     assay = 'ATAC', expression.assay = 'SCT', 
                     extend.upstream = 5000, extend.downstream = 5000) +
      plot_annotation(title = g, 
                      theme = theme(plot.title = element_text(size = 25, hjust = .5, vjust = .5)))

    plots.list[[g]] = p
  }
})

if ( !dir.exists('signi.peaks.by.cluster.by.cond/tracks') ) 
  dir.create('signi.peaks.by.cluster.by.cond/tracks')

### Individual files ###
system.time({
  for ( i in 1:length(plots.list) ) {
    cat(i, '\n')
    ggsave(plot = plots.list[[i]], width = 10, height = 12, 
           paste0('tracks/merged/peak.coverage.merged.', names(plots.list)[i], '.png'))
  }
})

### A single pdf file ###
system.time({
  pdf('tracks/peak.coverage.merged.all.pdf', width = 10, height = 12)
  invisible(lapply(plots.list, print))
  dev.off()
})



####################################################################################################

##### A dot plot of p-values for enriched motifs/TFs unique for each condition in each cluster #####

####################################################################################################

#top.tfs.unq = read.csv('ANL001b.merged.wt.mut.peaks.by.cluster.by.cond.signi.peaks.motif.top.unq.csv')
top.tfs.unq$gene = paste0(top.tfs.unq$gene, '|', top.tfs.unq$motif.feature)

pval.mtx = matrix(NA, 
                 nrow = length(unique(top.tfs.unq$gene)), 
                 ncol = length(unique(top.tfs.unq$cluster)), 
                 dimnames = list(sort(unique(top.tfs.unq$gene)), 
                                 sort(unique(top.tfs.unq$cluster))))

for ( i in 1:nrow(top.tfs.unq) ) {
  tmp1 = top.tfs.unq$gene[i]
  tmp2 = as.character(top.tfs.unq$cluster[i])
  tmp3 = top.tfs.unq$Cond[i]
  if ( tmp3 == 'Mut' ) {
    tmp4 = -log10(top.tfs.unq$motif.pval[i])
  } else tmp4 = log10(top.tfs.unq$motif.pval[i])
  
  pval.mtx[tmp1, tmp2] = tmp4
}

padj.mtx = pval.mtx
for ( i in 1:nrow(top.tfs.unq) ) {
  tmp1 = top.tfs.unq$gene[i]
  tmp2 = as.character(top.tfs.unq$cluster[i])
  tmp3 = top.tfs.unq$Cond[i]
  if ( tmp3 == 'Mut' ) {
    tmp4 = -log10(top.tfs.unq$motif.padj[i])
  } else tmp4 = log10(top.tfs.unq$motif.padj[i])
  
  padj.mtx[tmp1, tmp2] = tmp4
}

auc.mtx = pval.mtx
for ( i in 1:nrow(top.tfs.unq) ) {
  tmp1 = top.tfs.unq$gene[i]
  tmp2 = as.character(top.tfs.unq$cluster[i])
  tmp3 = top.tfs.unq$Cond[i]
  if ( tmp3 == 'Mut' ) {
    tmp4 = top.tfs.unq$motif.auc[i]
  } else tmp4 = -top.tfs.unq$motif.auc[i]
  
  auc.mtx[tmp1, tmp2] = tmp4
}

{
  corrplot(pval.mtx, is.corr = F, method = 'circle', col = c('orange', 'blue'), 
           na.label = 'square', na.label.col = 'white', 
           tl.col = 'black', tl.cex = .5, 
           addgrid.col = 'beige', 
           cl.pos = 'n')

  corrplot(padj.mtx, is.corr = F, method = 'circle', col = c('orange', 'blue'), 
           na.label = 'square', na.label.col = 'white', 
           tl.col = 'black', tl.cex = .5, 
           addgrid.col = 'beige', 
           cl.pos = 'n')

  corrplot(auc.mtx, is.corr = F, method = 'circle', col = c('orange', 'blue'), 
           na.label = 'square', na.label.col = 'white', 
           tl.col = 'black', tl.cex = .5, 
           addgrid.col = 'beige', 
           cl.pos = 'n')
}

{
  pval.data = NULL
  for ( i in 1:nrow(top.tfs.unq) ) {
    tmp1 = top.tfs.unq$gene[i]
    tmp2 = as.character(top.tfs.unq$cluster[i])
    tmp3 = top.tfs.unq$Cond[i]
    if ( tmp3 == 'Mut' ) {
      tmp4 = -log10(top.tfs.unq$motif.pval[i])
    } else tmp4 = log10(top.tfs.unq$motif.pval[i])
    
    pval.data = rbind(pval.data, c('gene' = tmp1, 'cluster' = tmp2, 'pval' = tmp4))
  }
  
  padj.data = NULL
  for ( i in 1:nrow(top.tfs.unq) ) {
    tmp1 = top.tfs.unq$gene[i]
    tmp2 = as.character(top.tfs.unq$cluster[i])
    tmp3 = top.tfs.unq$Cond[i]
    if ( tmp3 == 'Mut' ) {
      tmp4 = -log(top.tfs.unq$motif.padj[i], 10)
    } else tmp4 = log(top.tfs.unq$motif.padj[i], 10)
    
    padj.data = rbind(padj.data, c('gene' = tmp1, 'cluster' = tmp2, 'pval' = tmp4))
  }

  pval.data = data.frame(pval.data)
  padj.data = data.frame(padj.data)
}

{
  padj.data1 = padj.data
  padj.data1$pval = abs(as.numeric(padj.data$pval))
  ggplot(padj.data1) + 
    geom_point(aes(cluster, gene, size = pval)) + 
    scale_size(range = c(1, 5)) +
    theme_bw()
}

{
  png('ANL001b.merged.enriched.motifs2tfs.unq.for.cond.by.cluster.padj.dots.png', width = 1000, height = 3000)
  par(omi = c(0, 1.7, 0, 0))
  plot(c(1, ncol(padj.mtx)), c(1, nrow(padj.mtx)), type = 'n', 
       xlab = '', ylab = '', axes = F, frame.plot = T, 
       ylim = c(0, nrow(padj.mtx)+1), yaxs = "i")
  axis(1, at = 1:ncol(padj.mtx), labels = colnames(padj.mtx), tick = F, las = 2, cex.axis = 1.8)
  axis(2, at = 1:nrow(padj.mtx), labels = rev(rownames(padj.mtx)), tick = F, las = 2, cex.axis = 1.5)
  abline(h = 1:nrow(padj.mtx), v = 1:ncol(padj.mtx), col = 'grey', lwd = .5)
  points(match(padj.data$cluster, colnames(padj.mtx)), match(padj.data$gene, rev(rownames(padj.mtx))), 
         pch = 19, 
         cex = log(abs(as.numeric(padj.data$pval)) / min(abs(as.numeric(padj.data$pval))), 5) + 1,
         col = sapply(as.numeric(padj.data$pval), function(x) if ( x < 0 ) 'orange' else 'blue')
  )
  dev.off()
}



####################################################################################################

##### A dot plot of p-values for enriched motifs/TFs unique for each condition in each cluster #####

####################################################################################################

#top.tfs.all.wt = read.csv('ANL001b.merged.wt.mut.peaks.by.cluster.by.cond.signi.peaks.wt.motif.top.csv')
#top.tfs.all.mut = read.csv('ANL001b.merged.wt.mut.peaks.by.cluster.by.cond.signi.peaks.mut.motif.top.csv')

top.tfs.all = rbind(cbind(Cond = 'WT', top.tfs.all.wt), cbind(Cond = 'Mut', top.tfs.all.mut))

{
  motif.class = sapply(top.tfs.all$motif.feature, 
                       function(x) {tmp = getMatrixByID(JASPAR2020, ID = x)@matrixClass; 
                       if ( is.null(tmp) ) NA else tmp})
  motif.famil = sapply(top.tfs.all$motif.feature, 
                       function(x) {tmp = getMatrixByID(JASPAR2020, ID = x)@tags$family; 
                       if ( is.null(tmp) ) NA else tmp})
  
  top.tfs.all = cbind(top.tfs.all, motif.class = motif.class, motif.family = motif.famil)

  write.csv(top.tfs.all, 'ANL001b.merged.wt.mut.peaks.by.cluster.by.cond.signi.peaks.motif.top.anno.csv', 
            quote = F, row.names = F)
}

{
  top.tfs.all$gene = paste0(top.tfs.all$gene, '|', top.tfs.all$motif.feature)
  top.tfs.all = top.tfs.all[order(top.tfs.all$gene), ]
}

pval.mtx = matrix(NA, 
                  nrow = length(unique(top.tfs.all$gene)), 
                  ncol = length(unique(top.tfs.all$cluster)), 
                  dimnames = list(sort(unique(top.tfs.all$gene)), 
                                  sort(unique(top.tfs.all$cluster))))

for ( i in 1:nrow(top.tfs.all) ) {
  tmp1 = top.tfs.all$gene[i]
  tmp2 = as.character(top.tfs.all$cluster[i])
  tmp3 = top.tfs.all$Cond[i]
  if ( tmp3 == 'Mut' ) {
    tmp4 = -log10(top.tfs.all$motif.pval[i])
  } else tmp4 = log10(top.tfs.all$motif.pval[i])
  
  pval.mtx[tmp1, tmp2] = tmp4
}

padj.mtx = pval.mtx
for ( i in 1:nrow(top.tfs.all) ) {
  tmp1 = top.tfs.all$gene[i]
  tmp2 = as.character(top.tfs.all$cluster[i])
  tmp3 = top.tfs.all$Cond[i]
  if ( tmp3 == 'Mut' ) {
    tmp4 = -log10(top.tfs.all$motif.padj[i])
  } else tmp4 = log10(top.tfs.all$motif.padj[i])
  
  padj.mtx[tmp1, tmp2] = tmp4
}

auc.mtx = pval.mtx
for ( i in 1:nrow(top.tfs.all) ) {
  tmp1 = top.tfs.all$gene[i]
  tmp2 = as.character(top.tfs.all$cluster[i])
  tmp3 = top.tfs.all$Cond[i]
  if ( tmp3 == 'Mut' ) {
    tmp4 = top.tfs.all$motif.auc[i]
  } else tmp4 = -top.tfs.all$motif.auc[i]
  
  auc.mtx[tmp1, tmp2] = tmp4
}

{
  corrplot(pval.mtx, is.corr = F, method = 'circle', col = c('orange', 'blue'), 
           na.label = 'square', na.label.col = 'white', 
           tl.col = 'black', tl.cex = .5, 
           addgrid.col = 'beige', 
           cl.pos = 'n')
  
  corrplot(padj.mtx, is.corr = F, method = 'circle', col = c('orange', 'blue'), 
           na.label = 'square', na.label.col = 'white', 
           tl.col = 'black', tl.cex = .5, 
           addgrid.col = 'beige', 
           cl.pos = 'n')
  
  corrplot(auc.mtx, is.corr = F, method = 'circle', col = c('orange', 'blue'), 
           na.label = 'square', na.label.col = 'white', 
           tl.col = 'black', tl.cex = .5, 
           addgrid.col = 'beige', 
           cl.pos = 'n')
}

{
  pval.data = NULL
  for ( i in 1:nrow(top.tfs.all) ) {
    tmp1 = top.tfs.all$gene[i]
    tmp2 = as.character(top.tfs.all$cluster[i])
    tmp3 = top.tfs.all$Cond[i]
    if ( tmp3 == 'Mut' ) {
      tmp4 = -log10(top.tfs.all$motif.pval[i])
    } else tmp4 = log10(top.tfs.all$motif.pval[i])
    
    pval.data = rbind(pval.data, c('gene' = tmp1, 'cluster' = tmp2, 'pval' = tmp4))
  }
  
  padj.data = NULL
  for ( i in 1:nrow(top.tfs.all) ) {
    tmp1 = top.tfs.all$gene[i]
    tmp2 = as.character(top.tfs.all$cluster[i])
    tmp3 = top.tfs.all$Cond[i]
    if ( tmp3 == 'Mut' ) {
      tmp4 = -log(top.tfs.all$motif.padj[i], 10)
    } else tmp4 = log(top.tfs.all$motif.padj[i], 10)
    
    padj.data = rbind(padj.data, c('gene' = tmp1, 'cluster' = tmp2, 'pval' = tmp4))
  }
  
  pval.data = data.frame(pval.data)
  padj.data = data.frame(padj.data)
  
  pval.data = pval.data[order(abs(as.numeric(pval.data$pval)), decreasing = T), ]
  padj.data = padj.data[order(abs(as.numeric(padj.data$pval)), decreasing = T), ]
  
  pval.data$id = paste0(pval.data$gene, '-', pval.data$cluster)
  padj.data$id = paste0(padj.data$gene, '-', padj.data$cluster)
}

{
  padj.data1 = padj.data
  padj.data1$pval = abs(as.numeric(padj.data$pval))
  ggplot(padj.data1) + 
    geom_point(aes(cluster, gene, size = pval)) + 
    scale_size(range = c(1, 5)) +
    theme_bw()
}

{
  dup.id = padj.data$id[duplicated(padj.data$id)]
  data.unq = padj.data[!padj.data$id %in% dup.id, ]

  data.dup = padj.data[padj.data$id %in% dup.id, ]
  data.dup = data.dup[order(data.dup$id), ]
}

{
  png('ANL001b.merged.enriched.motifs2tfs.all.by.cluster.padj.dots.v1.png', width = 1500, height = 6000)
  par(omi = c(0, 1.7, 0, 0))
  plot(c(1, ncol(padj.mtx)), c(1, nrow(padj.mtx)), type = 'n', 
       xlab = '', ylab = '', axes = F, frame.plot = T, 
       ylim = c(0, nrow(padj.mtx)+1), yaxs = "i")
  axis(1, at = 1:ncol(padj.mtx), labels = colnames(padj.mtx), tick = F, las = 2, cex.axis = 1.8)
  axis(2, at = 1:nrow(padj.mtx), labels = rev(rownames(padj.mtx)), tick = F, las = 2, cex.axis = 1.5)
  abline(h = 1:nrow(padj.mtx), v = 1:ncol(padj.mtx), col = 'grey', lwd = .5)
  points(match(padj.data$cluster, colnames(padj.mtx)), match(padj.data$gene, rev(rownames(padj.mtx))), 
         pch = 19, 
         cex = log(abs(as.numeric(padj.data$pval)) / min(abs(as.numeric(padj.data$pval))), 5) + 1,
         col = sapply(as.numeric(padj.data$pval), function(x) if ( x < 0 ) 'orange' else 'blue'))
  dev.off()
}

{
  png('ANL001b.merged.enriched.motifs2tfs.all.by.cluster.padj.dots.v2.png', width = 1500, height = 6000)
  par(omi = c(0, 1.7, 0, 0))
  plot(c(1, ncol(padj.mtx)), c(1, nrow(padj.mtx)), type = 'n', 
       xlab = '', ylab = '', axes = F, frame.plot = T, 
       ylim = c(0, nrow(padj.mtx)+1), yaxs = "i")
  axis(1, at = 1:ncol(padj.mtx), labels = colnames(padj.mtx), tick = F, las = 2, cex.axis = 1.8)
  axis(2, at = 1:nrow(padj.mtx), labels = rev(rownames(padj.mtx)), tick = F, las = 2, cex.axis = 1.5)
  abline(h = 1:nrow(padj.mtx), v = 1:ncol(padj.mtx), col = 'grey', lwd = .5)
  
  points(match(data.unq$cluster, colnames(padj.mtx)), match(data.unq$gene, rev(rownames(padj.mtx))), 
         pch = 19, 
         cex = log(abs(as.numeric(data.unq$pval)) / min(abs(as.numeric(data.unq$pval))), 5) + 1,
         col = sapply(as.numeric(data.unq$pval), function(x) if ( x < 0 ) 'orange' else 'blue'))
  
  data.dup.1 = data.dup[data.dup$pval < 0, ]
  head(data.dup.1)
  points(match(data.dup.1$cluster, colnames(padj.mtx))-0.3, match(data.dup.1$gene, rev(rownames(padj.mtx))), 
         pch = 19, 
         cex = log(abs(as.numeric(data.dup.1$pval)) / min(abs(as.numeric(data.dup.1$pval))), 5) + 1,
         col = 'orange')
  
  data.dup.2 = data.dup[data.dup$pval > 0, ]
  head(data.dup.2)
  points(match(data.dup.2$cluster, colnames(padj.mtx))+0.3, match(data.dup.2$gene, rev(rownames(padj.mtx))), 
         pch = 19, 
         cex = log(abs(as.numeric(data.dup.2$pval)) / min(abs(as.numeric(data.dup.2$pval))), 5) + 1,
         col = 'blue')
  dev.off()
}

{ ### Filtering v1 ###
  dt.splits = split(top.tfs.all, top.tfs.all$cluster)
  dt.splits = lapply(dt.splits, function(x) split(x, x$Cond))

  for ( i in 1:length(dt.splits) ) {
    tmp1 = dt.splits[[i]]
    for ( j in 1:length(tmp1) ) {
      tmp2 = tmp1[[j]]
      if ( nrow(tmp2) == 1 ) next
      tmp3 = table(sub('[.][0-9]$', '', tmp2$gene))
      tmp4 = names(which(tmp3 > 1))
      if ( length(tmp4) == 0 ) next
      tmp5 = names(which(tmp3 == 1))
      tmp.unq = tmp2[match(tmp5, sub('[.][0-9]$', '', tmp2$gene)), ]
      for ( k in 1:length(tmp4) ) {
        tmp6 = tmp2[grep(tmp4[k], tmp2$gene), ]
        tmp7 = tmp6[which.min(tmp6$motif.pval), ]
        tmp.unq = rbind(tmp.unq, tmp7)
      }
      
      dt.splits[[i]][[j]] = tmp.unq
    }
  }
  
  top.tfs.filt = NULL
  for ( i in 1:length(dt.splits) ) {
    tmp1 = dt.splits[[i]]
    tmp2 = NULL
    for ( j in 1:length(tmp1) ) {
      tmp2 = rbind(tmp2, tmp1[[j]])
    }
    top.tfs.filt = rbind(top.tfs.filt, tmp2)
  }
}

{ ### Filtering v2 ###
  dt.splits = split(top.tfs.all, sub('[.][0-9]$', '', top.tfs.all$gene))

  top.tfs.filt = NULL
  for ( i in 1:length(dt.splits) ) {
    tmp1 = dt.splits[[i]]
    if ( length(unique(tmp1$gene)) == 1 ) {
      top.tfs.filt = rbind(top.tfs.filt, tmp1)
      next
    }
    
    tmp2 = split(tmp1, tmp1$gene)
    tmp3 = sapply(tmp2, function(x) min(x$motif.pval))
    tmp4 = tmp1[tmp1$gene %in% names(which.min(tmp3)), ]
    
    top.tfs.filt = rbind(top.tfs.filt, tmp4)
  }
}

{ ### Filtering v3: top motif per motif class ###
  top.tfs.filt = top.tfs.filt[top.tfs.filt$motif.class != '', ]

  dt.splits = split(top.tfs.filt, top.tfs.filt$motif.class)

  par(omi = c(0, 2.6, 0, 0))
  barplot(sort(sapply(dt.splits, nrow), decreasing = F), las = 2, horiz = T, 
          xlab = 'Number of motifs')
  
  top.tfs.filt = NULL
  for ( i in 1:length(dt.splits) ) {
    tmp1 = dt.splits[[i]]
    if ( length(unique(tmp1$gene)) == 1 ) {
      top.tfs.filt = rbind(top.tfs.filt, tmp1)
      next
    }
    
    tmp2 = split(tmp1, tmp1$gene)
    tmp3 = sapply(tmp2, function(x) min(x$motif.pval))
    tmp4 = tmp1[tmp1$gene %in% names(which.min(tmp3)), ]
    
    top.tfs.filt = rbind(top.tfs.filt, tmp4)
  }
}

{ ### Filtering v4: top motif per cluster ###
  dt.splits = split(top.tfs.all, top.tfs.all$cluster)

  par(omi = c(0, .2, 0, 0))
  barplot(sort(sapply(dt.splits, nrow), decreasing = F), las = 2, horiz = T, 
          xlab = 'Number of motifs', ylab = 'Cluster')
  
  top.tfs.filt = NULL
  for ( i in 1:length(dt.splits) ) {
    tmp1 = dt.splits[[i]]
    if ( length(unique(tmp1$gene)) == 1 ) {
      top.tfs.filt = rbind(top.tfs.filt, tmp1)
      next
    }
    
    tmp2 = split(tmp1, tmp1$gene)
    tmp3 = sapply(tmp2, function(x) min(x$motif.pval))
    tmp4 = tmp1[tmp1$gene %in% names(which.min(tmp3)), ]
    
    top.tfs.filt = rbind(top.tfs.filt, tmp4)
  }
}

{
  padj.mtx.filt = matrix(NA, 
                         nrow = length(unique(top.tfs.filt$gene)), 
                         ncol = length(unique(top.tfs.filt$cluster)), 
                         dimnames = list(sort(unique(top.tfs.filt$gene)), 
                                         sort(unique(top.tfs.filt$cluster))))

  for ( i in 1:nrow(top.tfs.filt) ) {
    tmp1 = top.tfs.filt$gene[i]
    tmp2 = as.character(top.tfs.filt$cluster[i])
    tmp3 = top.tfs.filt$Cond[i]
    if ( tmp3 == 'Mut' ) {
      tmp4 = -log10(top.tfs.filt$motif.padj[i])
    } else tmp4 = log10(top.tfs.filt$motif.padj[i])
    
    padj.mtx.filt[tmp1, tmp2] = tmp4
  }
}

{
  padj.data.filt = NULL
  for ( i in 1:nrow(top.tfs.filt) ) {
    tmp1 = top.tfs.filt$gene[i]
    tmp2 = as.character(top.tfs.filt$cluster[i])
    tmp3 = top.tfs.filt$Cond[i]
    if ( tmp3 == 'Mut' ) {
      tmp4 = -log(top.tfs.filt$motif.padj[i], 10)
    } else tmp4 = log(top.tfs.filt$motif.padj[i], 10)
    
    padj.data.filt = rbind(padj.data.filt, c('gene' = tmp1, 'cluster' = tmp2, 'pval' = tmp4))
  }

  padj.data.filt = data.frame(padj.data.filt)
  padj.data.filt = padj.data.filt[order(abs(as.numeric(padj.data.filt$pval)), decreasing = T), ]
  padj.data.filt$id = paste0(padj.data.filt$gene, '-', padj.data.filt$cluster)
}

{
  dup.id = padj.data.filt$id[duplicated(padj.data.filt$id)]

  data.unq = padj.data.filt[!padj.data.filt$id %in% dup.id, ]

  data.dup = padj.data.filt[padj.data.filt$id %in% dup.id, ]
  data.dup = data.dup[order(data.dup$id), ]
}

{
  png('ANL001b.merged.enriched.motifs2tfs.all.by.cluster.padj.filt.top.per.cluster.dots.v1.png', 
      width = 1000, height = 1000)
  par(omi = c(0, 1.7, 0, 0))
  plot(c(1, ncol(padj.mtx.filt)), c(1, nrow(padj.mtx.filt)), type = 'n', 
       xlab = '', ylab = '', axes = F, frame.plot = T, 
       ylim = c(0, nrow(padj.mtx.filt)+1), yaxs = "i")
  axis(1, at = 1:ncol(padj.mtx.filt), labels = colnames(padj.mtx.filt), tick = F, las = 2, cex.axis = 1.5)
  axis(2, at = 1:nrow(padj.mtx.filt), labels = rev(rownames(padj.mtx.filt)), tick = F, las = 2, cex.axis = 1.5)
  abline(h = 1:nrow(padj.mtx.filt), v = 1:ncol(padj.mtx.filt), col = 'grey', lwd = .5)
  points(match(padj.data.filt$cluster, colnames(padj.mtx.filt)), match(padj.data.filt$gene, rev(rownames(padj.mtx.filt))), 
         pch = 19, 
         cex = log(abs(as.numeric(padj.data.filt$pval)) / min(abs(as.numeric(padj.data.filt$pval))), 5) + 1,
         col = sapply(as.numeric(padj.data$pval), function(x) if ( x < 0 ) 'orange' else 'blue'))
  dev.off()
}

{
  png('ANL001b.merged.enriched.motifs2tfs.all.by.cluster.padj.filt.top.per.cluster.dots.v2.png', 
      width = 1000, height = 1000)
  par(omi = c(0, 1.7, 0, 0))
  plot(c(1, ncol(padj.mtx.filt)), c(1, nrow(padj.mtx.filt)), type = 'n', 
       xlab = '', ylab = '', axes = F, frame.plot = T, 
       ylim = c(0, nrow(padj.mtx.filt)+1), yaxs = "i")
  axis(1, at = 1:ncol(padj.mtx.filt), labels = colnames(padj.mtx.filt), tick = F, las = 2, cex.axis = 1.5)
  axis(2, at = 1:nrow(padj.mtx.filt), labels = rev(rownames(padj.mtx.filt)), tick = F, las = 2, cex.axis = 1.5)
  abline(h = 1:nrow(padj.mtx.filt), v = 1:ncol(padj.mtx.filt), col = 'grey', lwd = .5)
  
  points(match(data.unq$cluster, colnames(padj.mtx.filt)), match(data.unq$gene, rev(rownames(padj.mtx.filt))), 
         pch = 19, 
         cex = log(abs(as.numeric(data.unq$pval)) / min(abs(as.numeric(data.unq$pval))), 5) + 1,
         col = sapply(as.numeric(data.unq$pval), function(x) if ( x < 0 ) 'orange' else 'blue'))
  
  data.dup.1 = data.dup[data.dup$pval < 0, ]
  head(data.dup.1)
  points(match(data.dup.1$cluster, colnames(padj.mtx.filt))-0.3, match(data.dup.1$gene, rev(rownames(padj.mtx.filt))), 
         pch = 19, 
         cex = log(abs(as.numeric(data.dup.1$pval)) / min(abs(as.numeric(data.dup.1$pval))), 5) + 1,
         col = 'orange')
  
  data.dup.2 = data.dup[data.dup$pval > 0, ]
  head(data.dup.2)
  points(match(data.dup.2$cluster, colnames(padj.mtx.filt))+0.3, match(data.dup.2$gene, rev(rownames(padj.mtx.filt))), 
         pch = 19, 
         cex = log(abs(as.numeric(data.dup.2$pval)) / min(abs(as.numeric(data.dup.2$pval))), 5) + 1,
         col = 'blue')
  dev.off()
}

{ ### Common and unique motifs/TFs ###
  dt.splits = split(top.tfs.all, top.tfs.all$cluster)
  dt.splits = lapply(dt.splits, function(x) split(x, x$Cond))

  dt.bygene = split(top.tfs.all, top.tfs.all$gene)

  (cluster.freq.bygene = sort(sapply(dt.bygene, function(x) length(unique(x$cluster))), decreasing = T))
}



###################################################################################################

##### Cluster-specific peaks : gene annotation within +/-5k for cluster 11 #####

###################################################################################################

#(sobj.filt = readRDS('ANL001b.merged.data.processed.motif.rds'))

cells.by.cluster = split(rownames(sobj.filt@meta.data), sobj.filt@meta.data$seurat_clusters)

CLUSTER = '11'

{
  atac.norm = data.matrix(sobj.filt[['ATAC']]@data)

  norm1 = atac.norm[, grep('WT', cells.by.cluster[[CLUSTER]], value = T)]
  norm2 = atac.norm[, grep('MUT', cells.by.cluster[[CLUSTER]], value = T)]

  norm.means = cbind(WT = rowMeans(norm1), MUT = rowMeans(norm2))
}

{
  annot = data.frame(Annotation(sobj.filt))
  annot$ID1 = apply(annot[, 1:3], 1, function(x) gsub(' +', '', paste(x, collapse = '-')))
  
  annot.unq = annot[!duplicated(annot$ID1), ]
}

{ ### Timepoint-specific peaks with no threshold ###
  peaks.unq1 = names(which(apply(norm.means, 1, function(x) x['WT'] > 0 & x['MUT'] == 0)))
  peaks.unq2 = names(which(apply(norm.means, 1, function(x) x['WT'] == 0 & x['MUT'] > 0)))
}

{ ### WT-unique peaks ###
  peaks.unq.tab1 = t(sapply(strsplit(peaks.unq1, '-'), function(x) as.numeric(x)))
  rownames(peaks.unq.tab1) = peaks.unq1
  colnames(peaks.unq.tab1) = c('Chr', 'Start', 'End')
  
  system.time({
    FLANKING = 5000
    tmp.cols = NULL
    for ( i in 1:nrow(peaks.unq.tab1) ) {
      cat(i, '\n')
      
      y1 = rep(0, nrow(annot.unq))
      y2 = which(annot.unq$seqnames == peaks.unq.tab1[i, 'Chr'])
      y1[y2] = as.numeric(peaks.unq.tab1[i, 'Start'] >= annot.unq$start[y2] - FLANKING & 
                            peaks.unq.tab1[i, 'End'] <= annot.unq$end[y2] + FLANKING)
      
      if ( sum(y1) == 0 ) next
      
      tmp.cols = c(tmp.cols, peaks.unq1[i])
      
      cat(y1, file = 'tmp.values', sep = '\n', append = T)
    }
    
    cat(tmp.cols, file = 'tmp.cols', sep = '\n')
  })

  system.time({
    tmp.values = scan('tmp.values')
  })

  system.time({
    peak.annot.mtx1 = matrix(tmp.values, nr = nrow(annot.unq), nc = length(tmp.cols), byrow = F)
    rownames(peak.annot.mtx1) = annot.unq$ID1
    colnames(peak.annot.mtx1) = tmp.cols
  })

  system.time({
    signalmap1.list = apply(peak.annot.mtx1, 2, function(x) annot.unq[x == 1, ])
  })

  signal.genes.freq1 = sort(table(unlist(sapply(signalmap1.list, 
                                                function(x) unique(x$gene_name)))), decreasing = T)

  signalmap1.genes = sapply(signalmap1.list, function(x) paste(unique(x$gene_name), collapse = '|'))
  
  system.time({
    save(peak.annot.mtx1, signalmap1.list, signalmap1.genes, 
         file = 'ANL001b.wt.mut.peaks.cluster11.unique.wt.annot.rdata')
  })

  file.remove('tmp.values', 'tmp.cols')

  peaks.unq1.int.genes = cbind(peaks.unq.tab1, 
                               norm.means[rownames(peaks.unq.tab1), ], 
                               Genes = signalmap1.genes[rownames(peaks.unq.tab1)])

  write.csv(peaks.unq1.int.genes, quote = F, 
            file = 'ANL001b.wt.mut.peaks.cluster11.unqiue.wt.allcells.mean.intensities.genes.csv')
}

{ ### Mutant-unique peaks ###
  peaks.unq.tab2 = t(sapply(strsplit(peaks.unq2, '-'), function(x) as.numeric(x)))
  rownames(peaks.unq.tab2) = peaks.unq2
  colnames(peaks.unq.tab2) = c('Chr', 'Start', 'End')
  
  system.time({
    FLANKING = 5000
    tmp.cols = NULL
    for ( i in 1:nrow(peaks.unq.tab2) ) {
      cat(i, '\n')
      
      y1 = rep(0, nrow(annot.unq))
      y2 = which(annot.unq$seqnames == peaks.unq.tab2[i, 'Chr'])
      y1[y2] = as.numeric(peaks.unq.tab2[i, 'Start'] >= annot.unq$start[y2] - FLANKING & 
                            peaks.unq.tab2[i, 'End'] <= annot.unq$end[y2] + FLANKING)
      
      if ( sum(y1) == 0 ) next
      
      tmp.cols = c(tmp.cols, peaks.unq2[i]) # peaks with gene annotations
      
      cat(y1, file = 'tmp.values', sep = '\n', append = T)
    }
    
    cat(tmp.cols, file = 'tmp.cols', sep = '\n') # For backup
  })

  system.time({
    tmp.values = scan('tmp.values')
  })

  system.time({
    peak.annot.mtx2 = matrix(tmp.values, nr = nrow(annot.unq), nc = length(tmp.cols), byrow = F)
    rownames(peak.annot.mtx2) = annot.unq$ID1
    colnames(peak.annot.mtx2) = tmp.cols
  })

  rm(tmp.values, tmp.cols)
  
  system.time({
    signalmap2.list = apply(peak.annot.mtx2, 2, function(x) annot.unq[x == 1, ])
  })

  signal.genes.freq2 = sort(table(unlist(sapply(signalmap2.list, 
                                                function(x) unique(x$gene_name)))), decreasing = T)

  signalmap2.genes = sapply(signalmap2.list, function(x) paste(unique(x$gene_name), collapse = '|'))

  system.time({
    save(peak.annot.mtx2, signalmap2.list, signalmap2.genes, 
         file = 'ANL001b.wt.mut.peaks.cluster11.unique.mut.annot.rdata')
  })

  file.remove('tmp.values', 'tmp.cols')

  peaks.unq2.int.genes = cbind(peaks.unq.tab2, 
                               norm.means[rownames(peaks.unq.tab2), ], 
                               Genes = signalmap2.genes[rownames(peaks.unq.tab2)])
  
  write.csv(peaks.unq2.int.genes, quote = F, 
            file = 'ANL001b.wt.mut.peaks.cluster11.unqiue.mut.allcells.mean.intensities.genes.csv')
}

