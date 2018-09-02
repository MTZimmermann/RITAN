# RITAN
RITAN is an R package for Rapid Integration of Term Annotation and Network resources.

RITAN brings together multiple resources to help answer two questions:

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

For additional examples of RITAN use, see our vignettes:

* `vignette('enrichment',package='RITAN')`
* `vignette('subnetworks',package='RITAN')`
* `vignette('choosing_resources',package='RITAN')`
* `vignette('resource_relationships',package='RITAN')`
* `vignette('multi_tissue_analysis',package='RITAN')`

<hr>

<b>Bringing More Data into RITAN</b>

To honor data redistribution policies, we provide certain annotation and network data in the RITANdata package. RITAN is most useful when many resources are organized for use. In my lab, we have many genesets, gene panels, pathway resources, gene regulatory networks, etc. organized and loaded into RITAN. As a middle-ground, we provide here recommendations for downloading other resources and the script we use for loading them into RITAN.

<b>work in progress...</b>
