# RITAN
RITAN is an R package for Rapid Integration of Term Annotation and Network resources.

RITAN brings together multiple resources to help answer two questions:
  1) "What functions are achieved by this group of genes?"
  2) "How do these genes work together to achieve that function?"

Specifically, RITAN has been written to:
  1) Identify terms/pathways/genesets represented by (enriched within) the input list of genes
  2) Identify network interactions among the genes in the input list (protein-protein, metabolic, and regulatory)

RITAN is to be applied to an unranked list of gene symbols and we perform false-discivery adjustment across term resources. As always, we encourage users to decide upon an analysis strategy, including consideration of which resources are most appropriate for their dataset/experiment, statistical significance thresholds, background geneset, etc., prior to analysis. The ease with which RITAN facilitates multi-resource query makes it easy to "try one more" test, leading to an increased number of hypothesis tests made that may not be accounted for by multiple testing correction. To prevent this, we encourage users to "include one more" resource - to add an additional resource to a single query in RITAN so that false-dicovery correction is maintained.

For additional examples of RITAN use, see our vignettes:
  vignette('enrichment',package='RITAN')
  vignette('subnetworks',package='RITAN')

