
#' write_simple_table
#'
#'This is a simple wrapper around "write.table" that writes a tab-delimited table with column names, no quoting, and no row names.
#'
#' @param d R data object
#' @param f file path
#' @param ... further options passed on to write.table
#'
#' @return invisible (nothing is returned)
#' @export
#'
#' @examples
#' \dontrun{
#' simple wrapper around write.table for writing a tab-delimieted, no row names, tab-seperated file
#' }
write_simple_table <- function(d=NULL,f=NULL,...){
  write.table(d,f, quote = FALSE, sep="\t", col.names = TRUE, row.names = FALSE,...)
  invisible(NULL)
}

### ---------------------------------------------------------------- -
### Get the network to use
### ---------------------------------------------------------------- -
#' readSIF
#'
#'This function reads a data table into R; the data table describes network interactions. It is named for the Simple Interaction Format (SIF), but can read any data table if the users identifies which columns contain the pertinent data (see below).
#'
#'The SIF file format is a 3-column format, with an optional 4th column:
#'<entity-1><tab><edge-type><tab><entity-2><tab><score>
#'
#'Entities may be genes, proteins, metabolites, etc. The edge type typically conveys the type of relationship that exists between the two entities, such as physical interaciton, phosphorylation, or activation.
#'
#' @param file location of file
#' @param header indicator of presense of header on file
#' @param sep file delimiter - used by read.table()
#' @param as.is logical (default TRUE)
#' @param p1 Column number for the 1st entity. Default = 1.
#' @param p2 Column number for the 2nd entity. Default = 2.
#' @param et Column number for the edge type. Default = 3.
#' @param score Column number for edge scores or weights. Default = NA (no score read).
#' @param ... Other options to read.table().
#'
#' @return Returns a data.frame with 3 (or 4) columns of data.
#' @export
#'
#' @examples
#' readSIF()
#' \dontrun{
#' sif <- readSIF("myFile.sif")
#' }
readSIF <- function( file = NA, header = FALSE, sep="\t", as.is=TRUE,
         p1=1, p2=2, et=3, score=NA, ... ){
  
  if(all(is.na(file))){
    print(sprintf('readSIF("myfile") where myfile is a tab-delimited file with 3 columns. Each row of the file describes an edge in a network. The first and third columns are genes/proteins/etc and the second column describes the interaction type. The optional fourth column contains a score or weight for the edge.'))
    return(invisible(NULL))
  }
  
  tmp <- read.table( file, header=header, sep=sep, as.is=as.is, ... )
  sif <- data.frame( p1 = tmp[,p1],
                     edge_type = tmp[,et],
                     p2 = tmp[,p2] )
  if (!all(is.na(score))){
    sif$score <- tmp[,score]
  }

  return(sif)

}

#' check_net_input
#'
#'A Quality Control function. This function will compare an input list of genes to a network reference and report if each member of the input is present in the resource.
#'
#' @param set An input list of genes to check against a reference.
#' @param ref A reference of network data. See readSIF().
#' @param check4similar Logical flag. If TRUE, a case-insensitive grep will be used for name matching. For genes in families with many related members (e.g. ABC*, FAM*, etc.), this will not be ideal. We intend this option as a QC screening method to identify if case, punctuaiton, etc is causing fewer than expected matches.
#' @param entity1name The column name in "ref" of the first entity. Default = "p1."
#' @param entity2name The column name in "ref" of the second entity. Default = "p2."
#'
#' @return Character vector of "yes/no" indicating "within-ref/not"
#' @export
#'
#' @examples
#' ## Return a "yes/no" vector indicating if each gene in myGeneSet is annotated with any term in GO
#' ## If no match, this function can attempt to suggest closest matches (check4similar = TRUE)
#' library(RITANdata)
#' myGeneSet <- c('BRCA1','RAD51C','VAV1','HRAS','ABCC1','CYP1B1','CYP3A5')
#' yorn <- check_net_input( myGeneSet, network_list[["CCSB"]] )
#' print(yorn)
#' 
#' yorn <- check_net_input( myGeneSet, network_list[["HPRD"]] )
#' print(yorn)
#' 
#' ## See check_any_net_input() for efficiently checking across all resources.
check_net_input <- function( set, ref, check4similar = FALSE, 
                             entity1name = "p1", entity2name = "p1" ){
  
  q <- sapply( set, function(x){
    h1 <- grepl( sprintf('^%s$',x), ref[[ entity1name ]], ignore.case = FALSE, perl=TRUE )
    h2 <- grepl( sprintf('^%s$',x), ref[[ entity2name ]], ignore.case = FALSE, perl=TRUE )
    if (any(h1|h2)){
      return('yes')
    } else {

      if (check4similar){
        h1 <- grepl( sprintf('%s',x), ref[[ entity1name ]], ignore.case = TRUE, perl=TRUE )
        h2 <- grepl( sprintf('%s',x), ref[[ entity2name ]], ignore.case = TRUE, perl=TRUE )
        r <- NULL

        if (any(h1)){ r <-      ref[[ entity1name ]][h1] ; }
        if (any(h2)){ r <- c(r, ref[[ entity2name ]][h2]); }

        if (!all(is.null(r))){
          return( sprintf('Potential Match: %s', paste( unique(r),sep=',',collapse = ',')) )
        }
      }

      return('no')
    }
  } )

  return(q)

}

#' check_any_net_input
#'
#'A Quality Control function. This function applies check_net_input() to all available resources (default).
#'
#' @param set An input list of genes to check against references.
#' @param Net2Use The collection of network resources to check within.
#' 
#' @return Logical vector indicating if the genes in "set" are within ANY of the resources.
#' @export
#'
#' @examples
#' #' ## Check if genes in myGeneSet are annotated by any resource in "network_list" (default).
#' library(RITANdata)
#' myGeneSet <- c('BRCA1','RAD51C','VAV1','HRAS','ABCC1','CYP1B1','CYP3A5')
#' yorn <- check_any_net_input( myGeneSet )
#' print(yorn)
check_any_net_input <- function(set, Net2Use = names(network_list) ){
  o <- lapply( network_list[Net2Use] , function(y){
      check_net_input( set, y, check4similar = FALSE )
    } )
  a <- do.call( cbind, o )
  y <- apply( a, 1, function(y){ any(y=='yes') } )
  return(y)
}


# ### ---------------------------------------------------------------- -
# ### Gene Network from shared GO term similarity
# ### ---------------------------------------------------------------- -
# # GOSimNET
# #
# # param g List of genes to generate symantic similarities among.
# # param symbols Logical flag for if "g" contains gene symbols or gene IDs.
# # param ont Which ontology within GO to use. Default is "BP" for Biologic Process.
# # param measure Distance measure to use. Defult is "Wang."
# # param mart BioMart mart object
# # param ... further arguments passed on to mgeneSim()
# # 
# # return The similarity matrix beteen all input pairs of "g"
# # 
# # 
# # examples
# # dontrun{
# # ## This is a wrapper around functions from the package "GOSemSim"
# # ## Example 1:
# # g <- c('BRCA1','BRCA2','GSTM1','PTEN','PIK3CA')
# # n <- GOSimNet( g )
# #
# # ## Example 2:
# # g <- c("835","5261","241","994")
# # n <- GOSimNet( g, symbols = FALSE )
# # }
# GOSimNet <- function( g = NA, symbols = TRUE, ont="BP", measure="Wang", mart = NULL, ... ){
# 
#     if ( (length(g)==1) & is.na(g) ){
#       stop('Please provide a list of genes, either as HGNC symbols or enterz IDs (set symbols=FALSE).')
#     }
# 
#     # source("https://bioconductor.org/biocLite.R")
#     # biocLite("GOSemSim")
#     # biocLite("org.Hs.eg.db")
#     require(GOSemSim)
#     ## If you use GOSemSim in published research, please cite: G Yu, F Li, Y Qin, X Bo, Y Wu, S Wang.
#     ##   GOSemSim: an R package for measuring semantic similarity among GO terms and gene products. Bioinformatcs 2010, 26(7):976-978.
# 
# 
#     if (symbols){
#       ## symbols given. Map to Entrez gene IDs
#       
#       require(biomaRt)
#       if (all(is.null(mart))){
#         ## get Mart, unless one was privided
#         mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="useast.ensembl.org")
#       }
#       
#       s <- g
#       to_get <- c('hgnc_symbol', 'entrezgene')
#       o <- getBM(attributes = to_get, filters = 'hgnc_symbol', values = s, mart = mart )
#       g <- as.character(unique(o$entrezgene))
#       
#     }
#   
#     m <- mgeneSim( genes = g, ont=ont, organism="human", measure=measure, verbose=FALSE, ...)
#     
#     if (symbols){
#       e <- as.character(o$entrezgene)
#       i <- match( e, g )
#       if(any( e[i] != g )){ stop('Name mismatch.') }
#       rownames(m) <- colnames(m) <- o$hgnc_symbol[i]
#     }
# 
#     return(m)
# 
#     # We  also  implemented clusterSim for  calculating  semantic  similarity  between  two  gene  clusters,  and mclusterSim for calculating semantic similarities among multiple gene clusters.
#     # GO  enrichment  analysis  can  be  supported  by  our  package clusterProfiler [10],  which  supports  hyper-geometric  test  and  Gene  Set  Enrichment  Analysis (GSEA).  Enrichment  results  across  different  gene clusters can be compared using compareCluster function.
# 
# }


### ---------------------------------------------------------------- -
### Get the network to use
### ---------------------------------------------------------------- -
#' SubSIF
#'
#' @param gene.list A list of genes to use. The function will identify edges across resources for or among these genes; identify the induced subnetwork around the gene.list.
#' @param Net2Use Name of network resource(s) to use.
#' @param minStringScore If STRING is among the Net2Use, only edges of at least the indicated score will be included. Default = 150.
#' @param minHumanNetScore If HumanNet is among the Net2Use, only edges of at least the indicated score will be included. Default = 0.4.
#' @param minScore Same as above, but used for any other networks where "score" is provided
#' @param verbose If TRUE (default), the function will update the user on what it is doing and how many edges are identified for each resource in Net2Use.
#' @param dedup If TRUE (Default = FALSE), 
#' @param directed.net Logical indicating if the network resources should be interpreted as directed. 
#' @param include_neighbors Logical to include 1st neighbors of "gene.list" (genes not in gene.list, but directly connected to them) in the induced subnetwork.
#'
#' @return Data table describing the induced subnetwork for "gene.list" across the requested resources.
#' @export
#'
#' @examples
#' ## Get interactions among a list of genes from the PID: Pathway Interaction Database
#' require(RITANdata)
#' myGeneSet <- c('BRCA1','RAD51C','VAV1','HRAS','ABCC1','CYP1B1','CYP3A5')
#' sif <- subSIF( myGeneSet, Net2Use = 'PID')
#' print(sif)
#' 
#' \dontrun{
#' ## Get the PPI network induced by genes within myGeneSet
#' ## Use 4 seperate resources, but trim STRING to only include more confident interactions
#' sif <- subSIF( myGeneSet, c('dPPI','HPRD','CCSB','STRING'), minStringScore = 500 )
#' }
subSIF <- function( gene.list      = NA,
                    Net2Use        = c('PID','TFe','dPPI','HPRD','CCSB','STRING'),
                    minStringScore = 150, # 7.8% have a score >= 0.7
                    minHumanNetScore = 0.4, # 0.5 ==> 84.8%, 1.0 ==> 44.7%, 1.5 ==> 25.8%, 2.0 ==> 12.9%, 2.5 ==> 6.6%, 3.0 ==> 3.0%
                    minScore       = 0,
                    verbose        = TRUE,
                    dedup          = FALSE,
                    directed.net   = FALSE,
                    include_neighbors = FALSE ) {

  return_full_network <- FALSE
  if (is.na(gene.list) || (class(gene.list) != 'character') || (length(gene.list) == 0) ){
    if (verbose){
      cat('** No input list provided.\n** The full edge list from requested networks will be returned.\n')
    }
    return_full_network <- TRUE
    net <- do.call( rbind, lapply( Net2Use, function(x){
      network_list[[x]][, c('p1','edge_type','p2') ] }) )
    return( net )
  }

  supportedNetworks <- names(network_list)
  i <- !(Net2Use %in% supportedNetworks)
  if (any(i)){
    w1 <- sprintf('  Unsupported Networks Requested: %s\n', paste(Net2Use[i], sep=',', collapse=',') )
    w2 <- sprintf('    Supported Networks Are: %s\n', paste(supportedNetworks, sep=',', collapse=',') )
    w  <- paste0( w1, w2 )
    warning( w )
  }

  if (verbose){
    cat(sprintf('\nGenerating %sdirected subnetwork...\n', ifelse( directed.net, '', 'un') ))
  }
  if (directed.net){ stop('directed networks not yet implemented.'); }
  
  n0 <- length(unique(gene.list))
  
  ## -------------------------- -
  select_edges <- function( x, inc_nei = include_neighbors ){
    
    y <- network_list[[x]][, c('p1','edge_type','p2') ]
    
    if (inc_nei){
      i <- (y$p1 %in% gene.list) | (y$p2 %in% gene.list)
    } else {
      i <- (y$p1 %in% gene.list) & (y$p2 %in% gene.list)
    }
    
    if ("score" %in% names(network_list[[x]]) ){
      if (x == 'HumanNet'){
        i <- i & (network_list[[x]]$score >= minHumanNetScore)
      } else if (x == 'STRING'){
        i <- i & (network_list[[x]]$score >= minStringScore)
      } else {
        i <- i & (network_list[[x]]$score >= minScore)
      }
    }
    
    return(y[i,])
    
  }
  
  ## -------------------------- -
  sif <- data.frame( p1=character(0), edge_type=character(0), p2=character(0))
  sif <- do.call( rbind, lapply( Net2Use , select_edges ) )
  
  if (include_neighbors){
    ## The first round is really to define which nodes/genes to consider.
    ## This second round will get all the edges.
    ## Without two rounds, we would just have "stars" around the input genes
    ##   and we'd miss the edges between neighbors.
    gene.list <- unique(c( sif$p1, sif$p2 ))
    sif <- do.call( rbind, lapply( Net2Use , function(x){
        select_edges(x, inc_nei = FALSE )
      } ) )
  }

  ## -------------------------- -
  if(dedup){
    if (verbose){
      cat('Removing duplicate edges and self-loops...\n')
    }

    ## easiest reduction first: full duplicates
    sif <- unique(sif)

    ## remove self loops
    ## NOTE: this can "remove" nodes from the network if the self-loop is the only edge they have
    i <- sif$p1 == sif$p2
    if (any(i)){
      sif <- sif[ !i, ]
    }

    ## collapse duplicate edges
    if (directed.net){
      # unique and directed. **** Loss of "source" information here...
      tmp <- sif[,c(1,3)]
    } else {
      # need to sort genes within each row - this preserves row order too
      tmp <- t(apply( sif, 1, function(x){ sort(x[c(1,3)]) } ))
    }
    i <- duplicated(tmp)
    if (any(i)){
      sif <- sif[ !i, ]
    }

  }

  if (verbose){
    n <- unlist(sif[,c(1,3)])
    if (!return_full_network){
      cat(sprintf('Total induced subnetwork from %d genes has %d nodes and %d edges (%d unique).\n',
                  n0, length(unique(n)), dim(sif)[1], dim(unique(sif[,c(1,3)]))[1] ))
    } else {
      cat(sprintf('Total network has %d nodes and %d edges.\n', length(unique(n)), dim(sif)[1] ))
    }
  }


  return(sif)

}



