
## Enrichment test for cell type markers identified using scRNA-seq in mouse hippocampus and cortex

To run enrichment, install addR - contains all relevant functions, including dependencies for celltr()
Data visualisation via Heat() - a wrapper for heatmap.2 function (included in addR), values='pval' converts P-values of enrichment to -log10(P-value) & adds R style significance values e.g. *** <0.0001

```
devtools::install_github("ks4471/addR")
library(adds)
```
Input is expected as a list of named gene sets, enrichment performed for each gene set. Example gene list available to download from dropbox https://www.dropbox.com/s/y2szo7ywloleh92/example.gene.list.Rdata?dl=0
```
input_gene_sets=list(gene_set1=c(,"ENSG00000167173","ENSG00000204161","ENSG00000120280","ENSG00000214212","ENSG00000121410","ENSG00000114779","ENSG00000168792"), gene_set_2=c("ENSG00000148584","ENSG00000198691","ENSG00000085563","ENSG00000136379","ENSG00000136754","ENSG00000177465","ENSG00000077080"))
cellt_enrich=celltr(input_gene_sets)
Heat(make.numeric(cellt_enrich),values='pval')
```
example output plot https://www.dropbox.com/s/e12xo5ugm4f5pyq/spatemp_dynamic.cellt.enrich.CDKL5mods.Pval.pdf?dl=0



### cell type markers converted to human ensembl gene ids using one-2-one orthologs in biomart:
Zeisel, A. et al. Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq. Science (80-. ). 25, 279–284 (2015)

### used in
1.  Srivastava, P. K. et al. A Systems-Level Framework For Drug Discovery Identifies Csf1R As A Novel Anti-Epileptic Drug Target. Nat. Comms. ISSN: 2041-1723 https://www.biorxiv.org/content/early/2017/05/22/140087
2.  Delahaye-Duriez, A. et al. Rare and common epilepsies converge on a shared gene regulatory network providing opportunities for novel antiepileptic drug discovery. Genome Biol. 17, 245 (2016).

