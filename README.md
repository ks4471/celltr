
# Enrichment test for cell type markers identified using scRNA-seq in mouse hippocampus and cortex

To run enrichment, install addR - contains all relevant functions, including dependencies for celltr()
Data visualisation via Heat() - a wrapper for heatmap.2 function (included in addR), values='pval' converts P-values of enrichment to -log10(P-value) & adds R style significance values e.g. *** <0.0001

```
devtools::install_github("ks4471/addR")
library(adds)
```
Input is expected as a list of named gene sets, enrichment performed for each gene set. Example gene list available to download from dropbox https://www.dropbox.com/s/y2szo7ywloleh92/example.gene.list.Rdata?dl=0
```
input_gene_sets=list(gene_set1=c("ENSG00000121410","ENSG00000114779","ENSG00000168792"),gene_set_2=c("ENSG00000148584","ENSG00000198691","ENSG00000085563"))
cellt_enrich=celltr(input_gene_sets)
Heat(make.numeric(cellt_enrich),values='pval')
```




# cell type markers converted to human ensembl gene ids using one-2-one orthologs in biomart:
Zeisel, A. et al. Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq. Science (80-. ). 25, 279–284 (2015)

