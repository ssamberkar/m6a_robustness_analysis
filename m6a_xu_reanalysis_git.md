## m6a reanalysis - Xu METTL3 dataset

The following analysis is the summary of the differential expression
analysis performed for the [Xu et
al.](https://doi.org/10.1038/s41586-021-03210-1) dataset. In the first
instance, the transcriptomic dataset, which formed a subset of the many
datasets in this publication, was downloaded from NCBI SRA with project
ID [PRJNA521368](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA521368).
Total QC, preprocessing and quantification was done using the
[NF-core](http://nf-co.re) RNAseq pipeline. Salmon counts obtained from
the pipeline were used for subsequent differential expression (DE)
analysis using DESeq2 using custom R script.

### Preprocess prior to DE analysis

We first list the files quantified by salmon and using the `tximport`
package, import in the workspace. Meanwhile, to determine gene-level
counts, we also generate a transcript &lt;–&gt; gene map from the
Ensembl mouse annotation package.

``` r
library(tximport)
library(DESeq2)
library(data.table)
library(dplyr)
library(tibble)
library(EnsDb.Mmusculus.v79)
library(GEOquery)
library(DT)
library(ggpubr)
library(ggplot2)
library(EnsDb.Hsapiens.v86)
library(GEOquery)
library(DT)
library(ggplot2)
library(kableExtra)
library(ggvenn)

# List count files
xu_cnt_files = list.files('/home/ssamberkar/work/projects/results/salmon', pattern = 'quant.sf', full.names = T, recursive = T)

# Generate tx-gene map for mouse from Ensembl Bioconductor

edb = EnsDb.Mmusculus.v79
mm_tx2gene = transcriptsBy(edb, by = "gene", filter = SeqNameFilter(c("X", "Y"))) %>%  
  data.frame() %>% 
  dplyr::select('tx_name','gene_id')
colnames(mm_tx2gene) = c('TXNAME', 'GENEID')
xu_txi = tximport(xu_cnt_files, type = 'salmon',tx2gene = mm_tx2gene, ignoreTxVersion = T)
```

## DE analysis

To obtain the experimental design, we fetch the experiment `phenoData`
originally submitted to GEO and subset it to generate our design matrix.
Once generated, we then use the `genotype` variable to run DE analysis.

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in design formula are characters, converting to factors

### DE results

Since there are more than 1 `genotype` tested in this dataset, we
generate separate DE tables for each pairwise contrast:

``` r
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
```

#### DE genes

The individual DE genes table are as follows:

``` r
kable(data.frame(Number_of_DEGs = unlist(lapply(xu_DE_results_df.list, function(x) dim(x)[1]))))
```

<table>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
Number\_of\_DEGs
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
parental \~ METTL3 KO
</td>
<td style="text-align:right;">
28
</td>
</tr>
<tr>
<td style="text-align:left;">
parental \~ METTL3 KO+WT METTL3
</td>
<td style="text-align:right;">
115
</td>
</tr>
<tr>
<td style="text-align:left;">
parental \~ METTL3 KO+Mut METTL3
</td>
<td style="text-align:right;">
23
</td>
</tr>
<tr>
<td style="text-align:left;">
METTL3 KO \~ METTL3 KO+WT METTL3
</td>
<td style="text-align:right;">
159
</td>
</tr>
<tr>
<td style="text-align:left;">
METTL3 KO \~ METTL3 KO+Mut METTL3
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
METTL3 KO+WT METTL3 \~ METTL3 KO+Mut METTL3
</td>
<td style="text-align:right;">
168
</td>
</tr>
</tbody>
</table>

### MA plots

``` r
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

for (i in 1:length(xu_ma_plots.list)){
  plot(xu_ma_plots.list[[i]])
}
```

    ## Warning: Removed 1662 rows containing missing values (geom_point).

    ## Warning: ggrepel: 2 unlabeled data points (too many overlaps). Consider increasing max.overlaps

![](/home/ssa18/m6a_robustness_analysis/m6a_xu_reanalysis_git_files/figure-gfm/xu_ma-1.png)<!-- -->

    ## Warning: Removed 1662 rows containing missing values (geom_point).

    ## Warning: ggrepel: 7 unlabeled data points (too many overlaps). Consider increasing max.overlaps

![](/home/ssa18/m6a_robustness_analysis/m6a_xu_reanalysis_git_files/figure-gfm/xu_ma-2.png)<!-- -->

    ## Warning: Removed 1662 rows containing missing values (geom_point).

    ## Warning: ggrepel: 4 unlabeled data points (too many overlaps). Consider increasing max.overlaps

![](/home/ssa18/m6a_robustness_analysis/m6a_xu_reanalysis_git_files/figure-gfm/xu_ma-3.png)<!-- -->

    ## Warning: Removed 1662 rows containing missing values (geom_point).

    ## Warning: ggrepel: 8 unlabeled data points (too many overlaps). Consider increasing max.overlaps

![](/home/ssa18/m6a_robustness_analysis/m6a_xu_reanalysis_git_files/figure-gfm/xu_ma-4.png)<!-- -->

    ## Warning: Removed 1662 rows containing missing values (geom_point).

![](/home/ssa18/m6a_robustness_analysis/m6a_xu_reanalysis_git_files/figure-gfm/xu_ma-5.png)<!-- -->

    ## Warning: Removed 1662 rows containing missing values (geom_point).

    ## Warning: ggrepel: 4 unlabeled data points (too many overlaps). Consider increasing max.overlaps

![](/home/ssa18/m6a_robustness_analysis/m6a_xu_reanalysis_git_files/figure-gfm/xu_ma-6.png)<!-- -->

### PCA plots

``` r
xu_pca_plots.list = vector(mode = 'list', length = length(xu_DE_results.list))
names(xu_pca_plots.list) = names(xu_DE_results.list)

# Applying rlog transform as number of samples < 30

xu_rlog_counts = rlogTransformation(xu_dds)

for(i in 1:length(xu_pca_plots.list)) {
  
  xu_pca_plots.list[[i]] = plotPCA(xu_rlog_counts[,which(xu_rlog_counts$genotype %in% c(xu_contrasts_pairs[i,1], xu_contrasts_pairs[i,2]))], intgroup = c('genotype')) + 
  theme_bw() + 
  ggtitle(label = 'PCA - Xu et al.', subtitle = paste(xu_contrasts_pairs[i,1], '~', xu_contrasts_pairs[i,2], sep=' '))
  
}

for (i in 1:length(xu_ma_plots.list)){
  plot(xu_pca_plots.list[[i]])
}
```

![](/home/ssa18/m6a_robustness_analysis/m6a_xu_reanalysis_git_files/figure-gfm/xu_pca-1.png)<!-- -->![](/home/ssa18/m6a_robustness_analysis/m6a_xu_reanalysis_git_files/figure-gfm/xu_pca-2.png)<!-- -->![](/home/ssa18/m6a_robustness_analysis/m6a_xu_reanalysis_git_files/figure-gfm/xu_pca-3.png)<!-- -->![](/home/ssa18/m6a_robustness_analysis/m6a_xu_reanalysis_git_files/figure-gfm/xu_pca-4.png)<!-- -->![](/home/ssa18/m6a_robustness_analysis/m6a_xu_reanalysis_git_files/figure-gfm/xu_pca-5.png)<!-- -->![](/home/ssa18/m6a_robustness_analysis/m6a_xu_reanalysis_git_files/figure-gfm/xu_pca-6.png)<!-- -->

# Comparison with Dominissini et al.

To next compare the 2 similar METTL3 KO datasets, we consider the DE
genes derived from `control ~ HepG2` comparisons from the [Dominissini
et al.](https://www.nature.com/articles/nature11112) and analysed
[separately](https://github.com/BarryDigby/GSE37001). Since the Xu
dataset uses mouse, we need to first convert the DE genes obtained for
all our contrasts to their human homologs. This can be achieved through
the [OMA browser](https://omabrowser.org/oma/home/) and inputting Human
and Mouse in the [Genome Pair
View](https://omabrowser.org/oma/genomePW/) form.

``` r
# Compare with Dominissini et al.
# Read human-mouse homologs as extracted from OMA

edb_human = EnsDb.Hsapiens.v86
ens_human_univ = keys(edb_human, keytype = 'GENEID')
oma_human_mouse_homologs = fread(input = "/home/ssa18/m6a_robustness_analysis/Mouse_Human_orthologs_OMA_output_20072021.txt", sep = '\t', header = F, stringsAsFactors = F, fill = T)

dom_mouse_degs = fread('/home/ssa18/mouseTxDEGs.tsv', sep = '\t', header = T, stringsAsFactors = F)

# First comparison with `parental vs METTL3 KO` and Dominissini et al.

ggvenn(data = list(Xu = oma_human_mouse_homologs$V1[oma_human_mouse_homologs$V2 %in% xu_DE_results_df.list$`METTL3 KO+WT METTL3 ~ METTL3 KO+Mut METTL3`$GeneID], Dominissini = dom_mouse_degs$Gene.stable.ID), fill_alpha = 0.3)
```

![](/home/ssa18/m6a_robustness_analysis/m6a_xu_reanalysis_git_files/figure-gfm/xu_dom_comp-1.png)<!-- -->

``` r
# Second comparison with `METTL3 KO+WT METTL3 ~ METTL3 KO+Mut METTL3` and Dominissini et al.

ggvenn(data = list(Xu = oma_human_mouse_homologs$V1[oma_human_mouse_homologs$V2 %in% xu_DE_results_df.list$`METTL3 KO+WT METTL3 ~ METTL3 KO+Mut METTL3`$GeneID], Dominissini = dom_mouse_degs$Gene.stable.ID), fill_alpha = 0.3)
```

![](/home/ssa18/m6a_robustness_analysis/m6a_xu_reanalysis_git_files/figure-gfm/xu_dom_comp-2.png)<!-- -->

## Significance of overlaps

While the number of DE genes between these 2 studies is quite low, the
significance of these overlaps had to be determined. For this, we ran
100k permutations by randomly sampling the same amount of DE genes as Xu
and Dominissini et al., in both the comparisons above. For each
iteration, an overlap between the Xu and Dominissini random DE gene set
was determined. The histograms of these overlapping DE genes are as
follows:

``` r
# Permutation (100k) test between Xu (parental vs. METTL3 knockdown cells) and Dominissini et al. (mouse <--> human)
sim1 = replicate(n = 100000, expr = length(intersect(sample(x = ens_human_univ, size = length(xu_DE_results_df.list$`parental ~ METTL3 KO`$GeneID), replace = F), 
                                                    sample(x = ens_human_univ, size = length(dom_mouse_degs$Gene.stable.ID), replace = F))), simplify = T)

sim2 = replicate(n = 100000, expr = length(intersect(sample(x = ens_human_univ, size = length(xu_DE_results_df.list$`METTL3 KO+WT METTL3 ~ METTL3 KO+Mut METTL3`$GeneID), replace = F), 
                                                     sample(x = ens_human_univ, size = length(dom_mouse_degs$Gene.stable.ID), replace = F))), simplify = T)

hist(sim1, main = paste("Xu, parental vs METTL3_KO ~ Dominissini"), xlab = 'Overlaps of random DE genes', col = '#abebc6')
```

![](/home/ssa18/m6a_robustness_analysis/m6a_xu_reanalysis_git_files/figure-gfm/de_sim-1.png)<!-- -->

``` r
hist(sim2, main = paste("Xu, METTL3 KO+WT METTL3 ~ METTL3 KO+Mut METTL3 ~ Dominissini"), xlab = 'Overlaps of random DE genes', col = '#abebc6')
```

![](/home/ssa18/m6a_robustness_analysis/m6a_xu_reanalysis_git_files/figure-gfm/de_sim-2.png)<!-- -->
