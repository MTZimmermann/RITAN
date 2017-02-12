## ----first_load, echo=TRUE, warning=FALSE--------------------------------
library(RITANdata)
library(RITAN)

## ----quick_start, echo=TRUE, warning=FALSE-------------------------------
my_genes <- c('ABCA2','ACAT2','ACSS2','CD9','CPEB2','CTNNB1','FASN','LDLR','LPL','LSS')
e <- term_enrichment(my_genes)
summary(e)

## ----echo_geneset_list_names, echo=TRUE----------------------------------
names(geneset_list)

## ----echo_geneset_list_names2, echo=TRUE---------------------------------
head(geneset_list$MSigDB_C7$GSE9988_LPS_VS_LOW_LPS_MONOCYTE_UP)

## ----make_study_selection, echo = TRUE-----------------------------------
selection <- grepl( 'GSE9988_(LOW_)*LPS_VS_.+UP', names(geneset_list$MSigDB_C7), perl=TRUE )
study_set <- geneset_list$MSigDB_C7[selection]
str(study_set)

## ----apply_enrichment_to_selection, echo=TRUE----------------------------
e <- term_enrichment_by_subset( study_set, q_value_threshold = 1e-5, 
                                term_sources = c("KEGG_filtered_canonical_pathways",
                                                 "MSigDB_Hallmarks") )

## ----ep1, echo=TRUE, fig.width = 7, fig.height = 6-----------------------
plot( e, show_values = FALSE, label_size_y = 7, label_size_x = 7 )

## ----ep1_cap, echo=TRUE, fig.width = 7, fig.height = 6-------------------
plot( e, show_values = FALSE, label_size_y = 7, label_size_x = 7, cap=10 )

## ----ep1_data, echo=TRUE, fig.width = 7, fig.height = 6------------------
prmatrix( e[1:3,], collab = c('name','n.set',1:7) )

## ----make_ann_mat, echo=TRUE---------------------------------------------
mat <- matrix(c("LPS","LPS","LPS","LPS","LOW_LPS","LOW_LPS","LOW_LPS",
                "LOW_LPS","LPS_AND_ANTI_TREM1","CTRL_TREATED",
                "VEHICLE_TREATED","ANTI_TREM1_AND_LPS","CTRL_TREATED","VEHICLE_TREATED"),
                nrow = 2, byrow = TRUE )
rownames(mat) <- c('Condition1','Condition2')
colnames(mat) <- sprintf('Sample%s', 1:7)
print(mat)

## ----show_ann_mat, echo=TRUE, fig.width = 7, fig.height = 8--------------
plot( e, show_values = TRUE, label_size_y = 7, label_size_x = 7, cap=10, 
      annotation_matrix = mat, grid_line_color = 'black' )

## ----apply_enrichment_to_selection_n, echo=TRUE--------------------------
n <- term_enrichment_by_subset( study_set, q_value_threshold = 1e-5,
                                term_sources = c("KEGG_filtered_canonical_pathways",
                                                 "MSigDB_Hallmarks"),
                                display_type = 'n', phred = FALSE )

## ----apply_enrichment_to_selection_n_plot, echo=TRUE, fig.width = 7, fig.height = 6----
plot( n, show_values = TRUE, label_size_y = 7, label_size_x = 7, cap = 20 )

## ----term_enrichment1, eval=FALSE----------------------------------------
#  data("vac1.day0vs31.de.genes")
#  te <- term_enrichment( geneset = vac1.day0vs31.de.genes )

## ----term_enrichment2, eval=FALSE----------------------------------------
#  e <- term_enrichment( geneset = vac1.day0vs31.de.genes, verbose = TRUE,
#                        term_sources = c("Blood_Translation_Modules", "MSigDB_C7") )

## ----term_enrichment3, eval=FALSE----------------------------------------
#  e <- term_enrichment( geneset = vac1.day0vs31.de.genes, verbose = TRUE,
#                        term_sources = names(geneset_list) )

## ----term_enrichment_add_gmt, echo=TRUE----------------------------------
gs  <- geneset_list$MSigDB_C7[['GSE6269_HEALTHY_VS_FLU_INF_PBMC_UP']]
gmt <- system.file("extdata", "curated_gene_disease_associations.gmt.gz", package="RITAN")

# -->> Not running here for brevity
# geneset_list$DisGeNet <- readGMT(gmt)
# str(head(geneset_list$DisGeNet))

## ----term_enrichment_provided_gmt, echo=TRUE-----------------------------
e2  <- term_enrichment( gs, term_sources = gmt )
print( e2[1:3,-5] )

## ----term_enrichment_search_gmt, echo=TRUE-------------------------------
geneset_list$DisGeNet <- readGMT(gmt)
print( geneset_list$DisGeNet[['Influenza, Human']] )

## ----show_hist, echo=TRUE------------------------------------------------
show_active_genesets_hist()
length(active_genesets)

## ----LDL, echo=TRUE------------------------------------------------------
geneset_list$LDL = list( LDL_import = c('APOB','APOE','LDLR'), 
                         LDL_processing = c('HMGR','ACAT2','HMGCS1',
                                            'HMGCR','MVD','MVK',
                                            'PMVK','IDI1','IDI2') )

## ----LDLe, eval=FALSE----------------------------------------------------
#  e <- term_enrichment( gs, term_sources = c('GO','LDL') )

## ----load_geneset_symbols.example, echo=TRUE-----------------------------
load_geneset_symbols()

load_geneset_symbols( gmt="ReactomePathways" )

## ----enrichment_symbols.example, echo = TRUE-----------------------------
data("vac2.day0vs31.de.genes")
enrich <- enrichment_symbols( geneset = vac2.day0vs31.de.genes )
head(enrich$name)

## ----term_enrichment.example, eval = FALSE-------------------------------
#  e.together <- term_enrichment( geneset = vac2.day0vs31.de.genes )
#  e.separate <- term_enrichment( geneset = vac2.day0vs31.de.genes,
#                                 report_resources_separately = TRUE )

