# RITAN
RITAN is an R package for Rapid Integration of Term Annotation and Network resources.

RITAN brings together multiple resources, primarily to help answer two questions:

1. "What functions are achieved by this group of genes?"
2. "How do these genes work together to achieve that function?"

Specifically, RITAN has been written to:

1. Identify terms/pathways/genesets represented by (enriched within) the input list of genes
2. Identify network interactions among the genes within the input list (protein-protein, metabolic, and regulatory)

We have also built RITAN to facilitate research questions such as:

1. What terms/pathways share similar genes?
2. How do two resource's definitions of each pathway compare? List their common and unique sets.
3. From my genes of interest, use network resources to add genes that connect them together or are neighbors.
4. Annotate a combined list of genes with multiple functions contributed to.
6. Generate a heatmap for pathway/geneset -by- condition; for each timepoint, treatment, etc., show how term enrichment or pathway activity changes.
5. Make linking these results to other tools, such as Cytoscape, simpler via generating the focused and annotated subnetwork.

<hr>

RITAN is to be applied to an unranked list of gene symbols and we perform false-discivery adjustment across all resources used. We encourage users to decide upon an analysis strategy prior to running their analysis, including consideration of which resources are most appropriate for their dataset/experiment, statistical significance thresholds, background geneset, etc. The ease with which RITAN facilitates multi-resource query may lead users to "try one more test," leading to an increased number of hypothesis tests made that may not be accounted for by multiple testing correction. To prevent this, we encourage users to "include one more resource" - to add an additional resource to a single query in RITAN so that false-dicovery correction is appropriately maintained.

These and other topics are discussed in our publication:
<ul>
<li>Zimmermann, M.T. <i>et al.</i> RITAN: Rapid Integration of Term Annotation and Network Resources. <i>in review</i>.</li>
</ul>

For additional examples of using RITAN, see our vignettes:

* `vignette('enrichment',package='RITAN')`
* `vignette('subnetworks',package='RITAN')`
* `vignette('choosing_resources',package='RITAN')`
* `vignette('resource_relationships',package='RITAN')`
* `vignette('multi_tissue_analysis',package='RITAN')`

<hr>

## For users interested in using RITAN, but who do not use R

We have made a standalone and self-contained executable available here:
https://mcw.box.com/s/bzmt12bgj2uygvcw6bk0hpg4p4w7nm2g
To use the app, download the .zip file, extract the zip file, and double-click the file, RITAN-electon.exe.

<p align="center">
  <object data="https://mcw.box.com/s/6y76mb22b06ovs3nxjocjz16al07o3up" width="350" type="image/pdf" title="Format for Group Comparison">
  <object data="https://mcw.box.com/s/zf3uqabksziu89uplabs2psnyonxfqwa" width="350" type="image/png" alt="Expected Result" title="Expected Result">
</p>

The standalone app uses Chromium-electron and stand-alone R packaged together using node.js.
<br>
For making the combination of Shiny, RPortable, and electron, thank you to: https://github.com/ksasso/Electron_ShinyApp_Deployment

<hr>

# Bringing More Data into RITAN

I make an R script that has all my setup within. That way, I can simply `source('~/setup_RITAN.R')` and I'm ready to go. If you use RITAN frequently, I suggest writing the `source()` command in your `.Rprofile`. Below, are excerpts from my setup script.

To honor data redistribution policies, we provide certain annotation and network data in the RITANdata package. RITAN is most useful when many resources are organized for use. In my lab, we have many genesets, gene panels, pathway resources, gene regulatory networks, etc. organized and loaded into RITAN.

As a middle-ground, we provide here recommendations for downloading other resources and excerpts from my resrouce loading script for adding them to RITAN. Eventually, we plan to add a setup() function to RITAN that would automate this process. For now, we hope the below examples are helpful to users. We welcome feedback on the package and how to work more efficiently with existing open-source solutions.

```{r eval=FALSE, echo=TRUE, tidy=TRUE}
require(RITANdata)
require(RITAN)

## -------------------------------- -
## Set a couple of utility functions

apath <- '.' # annotation data path - script assumes annotation data are within a common path

csv_split <- function(x){ sort(unique(strsplit( gsub( ' ', '', x ) ,",")[[1]] )) }

read_dlm_as_sif <- function( f = NA, sep = "\t",
                             label_column = 1, gene_column = 2, ...){

  warning(sprintf('Reading delimited file into SIF format.\nAssuming labels are in column %d and genes are in column %d.\n', label_column, gene_column))

  d <- read.table( file=f, sep=sep, ... )
  l <- unique( d[ , label_column ] )
  sif <- list()
  for (n in l){
    i <- ( d[ , label_column ] == n )
    sif[[n]] <- d[ i, gene_column ] %>% unique() %>% sort()
  }

  return(sif)

}

## -------------------------------- -
## Add additional networks to RITAN

network_list$RegNet  <- readSIF( paste(apath, 'RegNetwork/human/human.source.sif', sep="/" ))
network_list$TRRUST  <- readSIF( paste(apath, 'TRRUST/trrust_rawdata.txt', sep="/" ))
network_list$IID     <- readSIF( sprintf('%s/IID/iid.human.2016-03.sif.gz',apath))
network_list$BioPlex <- readSIF( sprintf('%s/BioPlex/BioPlex.sif.gz',apath), et='BioPlex', score=3)
network_list$HPRD    <- readSIF( sprintf('%s/HPRD/HPRD_Release9_062910/BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt', apath), p1=1, p2=4, et=7 )
network_list$BioGRID <- readSIF( sprintf('%s/BioGRID/BIOGRID-ORGANISM-Homo_sapiens-3.4.163.tab2.txt.gz', apath), p1=8, p2=9, et=12, score=19, quote='', skip=1, comment.char='' )
network_list$BioGRID$score <- as.numeric(network_list$BioGRID$score)

## -------------------------------- -
## Add additional genesets to RITAN

CPIC   <- read.csv( 'https://api.pharmgkb.org/v1/download/file/data/cpicPairs.csv', skip=1, header=TRUE, quote='', as.is=TRUE )
COSMIC <- read.table( sprintf("%s/COSMIC/cancer_gene_census_v83_GRCh38.csv",apath), sep=",", header=TRUE, as.is=TRUE)

geneset_list$MSigDB_C6     <- readGMT(sprintf('%s/genesets/c6.all.v6.1.symbols.gmt', apath))
geneset_list$CPIC_LevelABC <- sort(unique(CPIC$Gene[ ! grepl('D', CPIC$CPIC.Level) ])) %>% list(.)
geneset_list$HPO           <- read_dlm_as_sif( sprintf('%s/HPO/ALL_SOURCES_FREQUENT_FEATURES_phenotype_to_genes.2018-05-03.txt.gz', apath), label_column = 1, gene_column = 4, quote = '' )
geneset_list$monarch       <- read_dlm_as_sif( sprintf('%s/monarch/gene_disease.9606.tsv.gz', apath), label_column = 6, gene_column = 2, header = TRUE, quote = '' )
geneset_list$COSMIC        <- sort(unique(as.character( COSMIC$Gene.Symbol ))) %>% list(.)
geneset_list$ACMG56        <- sort(unique(strsplit("BRCA1 BRCA2 TP53 STK11 MLH1 MSH2 MSH6 PMS2 APC MUTYH VHL MEN1 RET PTEN RB1 SDHD SDHAF2 SDHC SDHB TSC1 TSC2 WT1 NF2 COL3A1 FBN1 TGFBR1 TGFBR2 SMAD3 ACTA2 MYLK MYH11 MYBPC3 MYH7 TNNT2 TNNI3 TPM1 MYL3 ACTC1 PRKAG2 GLA MYL2 LMNA RYR2 PKP2 DSP DSC2 TMEM43 DSG2 KCNQ1 KCNH2 SCN5A LDLR APOB PCSK9 RYR1 CACNA1S", " ")[[1]])) %>% list(.)

```
