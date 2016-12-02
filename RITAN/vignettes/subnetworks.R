## ----load_data from package, echo=TRUE-----------------------------------
library(RITAN)
library(RITANdata)
data(package="RITANdata", network_list)

## ----citation, echo=TRUE-------------------------------------------------
require(knitr)
kable( attr(network_list, 'network_data_sources') )

## ----example1, echo=TRUE-------------------------------------------------
my_genes <- geneset_list$MSigDB_C7[['GSE6269_HEALTHY_VS_FLU_INF_PBMC_DN']]
net <- subSIF( my_genes )

## ----example1.1, echo=TRUE-----------------------------------------------
net[ 1:5, ]

## ----example2, echo=TRUE-------------------------------------------------
net2 <- subSIF( my_genes, minStringScore = 700 )
str(net2)

## ----example3_1, echo=TRUE-----------------------------------------------
my_genes <- geneset_list$MSigDB_C2[['VERNOCHET_ADIPOGENESIS']]

net3.1 <- subSIF( my_genes, include_neighbors = FALSE, dedup = TRUE )
str(net3.1)
# 
# net3.2 <- subSIF( my_genes, include_neighbors = TRUE, dedup = TRUE )
# str(net3.2)
net3.2 <- subSIF( my_genes, include_neighbors = TRUE, dedup = TRUE, minStringScore = 990 )

## ----example3_3a, echo=TRUE----------------------------------------------
net3.3 <- subSIF( my_genes, include_neighbors = TRUE, dedup = TRUE, minStringScore = 990 )
str(net3.3)

## ----example3_3, echo=TRUE-----------------------------------------------
net3.3 <- subSIF( my_genes, include_neighbors = TRUE, dedup = TRUE, Net2Use = c('HPRD','CCSB') )
str(net3.3)

## ----example4, fig.height = 5, fig.width = 5, fig.align = 'center'-------
require(igraph)
edges <- as.matrix( net3.1[, c(1,3)] )
G <- graph( c(t(edges)), directed = FALSE )
par(mar=rep(0,4))
plot(G, vertex.size = 20, vertex.frame.color = 'white' )

## ----example5, echo=TRUE, eval=FALSE-------------------------------------
#  g <- unique(c( net3.1$p1, net3.1$p2 ))
#  
#  tab <- data.frame( gene      = c('FABP4',  'CEBPA','PPARG','ADRB3','RETN','AGT','HP',
#                                   'RARRES2','PANK3','FFAR2','LUM',  'MC2R','ADCYAP1R1'),
#                     TrogRatio = c( 1.8,      1.7,    0.6,    0.3 ,   0.3 ,  0.4  ,0.2,
#                                    0.3,      0.1,    0.5,    0.3 ,   0.5 ,  0.1),
#                     WAT_BAT   = c( 0.8,      1.0,    0.6,    10.0,   21.6,  215.4,2.4,
#                                    9.5,      3.9,    4.6,    4.0 ,   7.3 ,  2.6),
#                     initial   = g %in% my_genes
#                     )
#  
#  write_simple_table(net3.1, 'net_example.sif')
#  write_simple_table(tab,    'net_example.tab')

## ----example6, echo=TRUE, eval=FALSE-------------------------------------
#  
#  ## Add a new resource to "network_list"
#  ### For brevity, we
#  network_list[['BioGRID_Mouse']] <- readSIF( 'BIOGRID-ORGANISM-Mus_musculus-3.4.136.symbols.sif.gz', header = TRUE )
#  # > str(network_list[['BioGRID_Mouse']])
#  # 'data.frame':	38322 obs. of  3 variables:
#  #  $ p1       : chr  "SMAD2" "SMAD2" "SMAD2" "SMAD2" ...
#  #  $ edge_type: chr  "physical" "physical" "physical" "physical" ...
#  #  $ p2       : chr  "Rasd2" "Rab34" "Rhebl1" "Rab38" ...
#  
#  ## Short example from Tang's 2010 Nature paper
#  my_mouse <- c('Sost','Fxyd4','Tmprss6','Crtap','Thpo','Kcnn4','Osm','Slc29a3','ALB')
#  
#  ## First, check if these genes appear in the BioGRID network.
#  check_net_input( my_mouse, network_list[['BioGRID_Mouse']] )
#  #  Sost   Fxyd4 Tmprss6   Crtap    Thpo   Kcnn4     Osm Slc29a3     ALB
#  # "yes"    "no"    "no"    "no"    "no"    "no"    "no"    "no"    "no"
#  
#  ## After correcting a few gene names, get the induced subnetwork from mouse data.
#  my_mouse <- c('Sost','Fxyd4','Tmprss6','CRTAP','Thpo','KCNN4','Osm','Slc29a3','ALB')
#  net.m <- subSIF( my_mouse, include_neighbors = TRUE, Net2Use = c('BioGRID_Mouse') )
#  str(net.m)
#  # Generating undirected subnetwork...
#  # Total induced subnetwork from 9 genes has 17 nodes and 17 edges (17 unique).
#  # 'data.frame':	17 obs. of  3 variables:
#  #  $ p1       : chr  "Sf3a1" "Nphp1" "Iqcb1" "Invs" ...
#  #  $ edge_type: chr  "physical" "physical" "physical" "physical" ...
#  #  $ p2       : chr  "CRTAP" "Invs" "Nphp1" "ALB" ...
#  
#  ## Also, check within BioGRD's human network
#  check_net_input( my_mouse, network_list[['BioGRID_Human']] )
#  # Sost   Fxyd4 Tmprss6   CRTAP    Thpo   KCNN4     Osm Slc29a3     ALB
#  # "no"    "no"    "no"   "yes"    "no"   "yes"    "no"    "no"   "yes"
#  
#  ## Note that gene symbols are case sensitive
#  my_mouse <- c('SOST','Fxyd4','Tmprss6','CRTAP','THPO','KCNN4','OSM','Slc29a3','ALB')
#  check_net_input( my_mouse, network_list[['BioGRID_Human']] )
#  #  SOST   Fxyd4 Tmprss6   CRTAP    THPO   KCNN4     OSM Slc29a3     ALB
#  # "yes"    "no"    "no"   "yes"    "no"   "yes"   "yes"    "no"   "yes"
#  
#  ## Get the induced subnetowrk from human data
#  net.h <- subSIF( my_mouse, include_neighbors = TRUE, Net2Use = c('BioGRID_Human') )
#  str(net.h)
#  # Generating undirected subnetwork...
#  # Total induced subnetwork from 9 genes has 224 nodes and 755 edges (634 unique).
#  # 'data.frame':	755 obs. of  3 variables:
#  #  $ p1       : chr  "MBIP" "SH3GL1" "TNNT1" "GFAP" ...
#  #  $ edge_type: chr  "physical" "physical" "physical" "physical" ...
#  #  $ p2       : chr  "MBIP" "SH3GL1" "TNNT1" "GRAP2" ...
#  

## ----example7, echo=TRUE, eval=TRUE--------------------------------------
net <- subSIF( 'FOXP3', include_neighbors = TRUE, minStringScore = 990 )
genes <- unique(c( net$p1, net$p2 ))
e1 <- term_enrichment( genes, "Blood_Translaiton_Modules", verbose=FALSE )
summary(e1)
e2 <- term_enrichment( genes, "ReactomePathways", verbose=FALSE )
summary(e2)

