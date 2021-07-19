library(tximport)
library(DESeq2)
library(data.table)
library(dplyr)
library(tibble)
library(EnsDb.Mmusculus.v79)
library(GEOquery)
library(DT)
library(ggplot2)

# List count files
xu_cnt_files = list.files('/home/ssamberkar/work/projects/results/salmon/', pattern = 'quant.sf', full.names = T, recursive = T)

# Generate tx-gene map for mouse from Ensembl Bioconductor

edb = EnsDb.Mmusculus.v79
mm_tx2gene = transcriptsBy(edb, by = "gene", filter = SeqNameFilter(c("X", "Y"))) %>%  
  data.frame() %>% 
  dplyr::select('tx_name','gene_id')
colnames(mm_tx2gene) = c('TXNAME', 'GENEID')
xu_txi = tximport(xu_cnt_files, type = 'salmon',tx2gene = mm_tx2gene, ignoreTxVersion = T)

# To generate exp.design, utilise the GEO dataset

xu_geo_data = getGEO(GEO = 'GSE126242')
tmp_coldata = pData(xu_geo_data$GSE126242_series_matrix.txt.gz)[,c('title','genotype:ch1')]
colnames(xu_txi$abundance) = colnames(xu_txi$counts) = colnames(xu_txi$length) = tmp_coldata$title

# Now generate final exp.design table as required by DESeq2

xu_coldata = data.frame(genotype = tmp_coldata[,-1])
rownames(xu_coldata) = colnames(xu_txi$abundance)

# Diff. expression analysis
xu_dds = DESeqDataSetFromTximport(txi = xu_txi, colData = xu_coldata, design = ~genotype)
xu_dds = DESeq(xu_dds)

# Generate all pairwise contrasts from exp. design matrix

xu_contrasts_pairs = t(combn(unique(xu_coldata$genotype),2))
xu_DE_results.list = xu_DE_results_df.list = vector(mode = 'list', length = dim(xu_contrasts_pairs)[1])
names(xu_DE_results.list) = names(xu_DE_results_df.list) = paste(xu_contrasts_pairs[,1], xu_contrasts_pairs[,2], sep = ' ~ ')

for(i in 1:dim(xu_contrasts_pairs)[1]){
  
  xu_DE_results.list[[i]] = results(xu_dds, contrast = c('genotype', xu_contrasts_pairs[i,1], xu_contrasts_pairs[i,2])) 
  xu_DE_results_df.list[[i]] = xu_DE_results.list[[i]] %>% 
    data.frame() %>% 
    rownames_to_column('GeneID') %>%
    dplyr::filter(padj <= 0.05 & abs(log2FoldChange) >= 1) %>%
    mutate(`-log10_pval` = round(-log10(pvalue),digits = 3),
           `-log10_adj_pval` = round(-log10(padj),digits = 3),
           `log2FC` = round(log2FoldChange,digits = 3),
           `stat` = round(stat, digits = 3),
           `lfcSE` = round(lfcSE, digits = 3),
           `BaseMean` = round(baseMean, digits = 3)) %>%
    dplyr::select(c('GeneID', 'BaseMean', 'log2FC', 'lfcSE', 'stat', '-log10_pval', '-log10_adj_pval'))
  
  datatable(data.frame(number_of_DEGs = unlist(lapply(xu_DE_results_df.list, function(x)dim(x)[1]))) %>% rownames_to_column('Contrasts'), caption = 'Summary of DEGs for all pairwise contrasts')
}

# Applying rlog transform as number of samples < 30

xu_rlog_counts = rlogTransformation(xu_dds)

# Plot PCA for all contrasts

xu_pca_plots.list = vector(mode = 'list', length = length(xu_DE_results.list))
names(xu_pca_plots.list) = names(xu_DE_results.list)

for(i in 1:length(xu_pca_plots.list)) {
  
  xu_pca_plots.list[[i]] = plotPCA(xu_rlog_counts[,which(xu_rlog_counts$genotype %in% c(xu_contrasts_pairs[i,1], xu_contrasts_pairs[i,2]))], intgroup = c('genotype')) + 
  theme_bw() + 
  ggtitle(label = 'PCA - Xu et al.', subtitle = paste(xu_contrasts_pairs[i,1], '~', xu_contrasts_pairs[i,2], sep=' '))
  
}

# Plot MA for all contrasts

xu_ma_plots.list = vector(mode = 'list', length = length(xu_DE_results.list))
names(xu_ma_plots.list) = names(xu_DE_results.list)

for (i in 1:length(xu_ma_plots.list)){
  
  res = results(xu_dds, contrast = c('genotype', xu_contrasts_pairs[i,1], xu_contrasts_pairs[i,2])) %>% data.frame()
  mm_symbol = ensembldb::select(x = edb, keys = rownames(res), columns = 'SYMBOL', keytype = 'GENEID')
  res$name = mm_symbol$SYMBOL
  res$detection_call = 0
  res$detection_call[abs(res$log2FoldChange) > 1 & res$padj <= 0.05] = 1
  res = res %>% dplyr::select(c('name', 'baseMean', 'log2FoldChange', 'padj', 'detection_call'))
  
  xu_ma_plots.list[[i]] = ggmaplot(res, fdr = 0.05, fc = 2, size = 2,
           palette = c("#B31B21", "#1465AC", "darkgray"),
           genenames = as.vector(res$name),
           legend = "top", top = 20,
           font.label = c("bold", 11),
           font.legend = "bold",
           font.main = "bold",
           ggtheme = ggplot2::theme_minimal()) + 
    ggtitle(label = 'Xu et al. - MA plot', subtitle = paste(xu_contrasts_pairs[i,1] ,'~', xu_contrasts_pairs[i,2], sep = ' '))
  
}

# Compare with Dominissini et al.

dom_mouse_degs = fread('./rss_shared/mouseTxDEGs.tsv', sep = '\t', header = T, stringsAsFactors = F)
