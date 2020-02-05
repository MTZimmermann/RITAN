### TO IMPLEMENT
## - distributional test between random netowrks and a bootstrapped estimate from your data better than current interconnectivity test? (would just be a more confident estimate of observation)

if(getRversion() >= "2.15.1"){
  utils::globalVariables( names=c('active_genesets', 'all_symbols', 'geneset_list'), package='RITAN')
}
resources.default = c("GO", "ReactomePathways", "KEGG_filtered_canonical_pathways", "MSigDB_Hallmarks")

### --------------------------------------------------------------- -
## Declare datasets
### --------------------------------------------------------------- -

#' This dataset is included as an example in the package:
#' 
#' @examples
#' \dontrun{
#'  #data("vac1.day0vs31.de.genes")
#'  te <- term_enrichment(geneset = vac1.day0vs31.de.genes)
#' }
#' 
#' @return differentially expressed genes at 31 days post-vaccination with vaccine1
#' @references \url{https://www.ncbi.nlm.nih.gov/pubmed/26755593}
"vac1.day0vs31.de.genes"

#' This dataset is included as an example in the package:
#' 
#' @examples
#' \dontrun{
#'  #data("vac1.day0vs56.de.genes")
#'  te <- term_enrichment(geneset = vac1.day0vs56.de.genes)
#' }
#' 
#' @return differentially expressed genes at 56 days post-vaccination with vaccine1
#' @references \url{https://www.ncbi.nlm.nih.gov/pubmed/26755593}
"vac1.day0vs56.de.genes"

#' This dataset is included as an example in the package:
#' 
#' @examples
#' \dontrun{
#'  #data("vac2.day0vs31.de.genes")
#'  te <- term_enrichment(geneset = vac2.day0vs31.de.genes)
#' }
#' 
#' @return differentially expressed genes at 31 days post-vaccination with vaccine2
#' @references \url{https://www.ncbi.nlm.nih.gov/pubmed/26755593}
"vac2.day0vs31.de.genes"

#' This dataset is included as an example in the package:
#' 
#' @examples
#' \dontrun{
#'  #data("vac2.day0vs56.de.genes")
#'  te <- term_enrichment(geneset = vac2.day0vs56.de.genes)
#' }
#' 
#' @return differentially expressed genes at 56 days post-vaccination with vaccine2
#' @references \url{https://www.ncbi.nlm.nih.gov/pubmed/26755593}
"vac2.day0vs56.de.genes"


## Enrichment functions - primary function is term_enrichment.
### --------------------------------------------------------------- -
#' readGMT
#'
#' Created for simplification of reading .gmt files into RITAN.
#'
#' @param f GMT file name. Please provide a full path if the file is not in the current working directory.
#'
#' @return A list() where the name of each entry is the term (first column of GMT file) and the value is a chr array of genes associated with the term.
#' @export
#'
#' @examples
#' # Make an example list() to show the GMT format
#' set <- list( term1=c('gene_name1','gene_name2'),
#'              term2=c('gene_name3','gene_name4','gene_name5') )
#' \dontrun{
#' # Write a GMT file for "set"
#' writeGMT( set, 'my_file.gmt' )
#' 
#' # Reading GMT files
#' geneset <- readGMT( 'my_file.gmt' )
#' 
#' # GMT files are available from multiple sources including http://software.broadinstitute.org/gsea/msigdb/
#' }
#' 
readGMT <- function( f = NA ){
  
  if(all(is.na(f))){
    stop('Please provide an input file:
         readGMT("myfile") where myfile is a tab-delimited file. Each row of the file describes a geneset. The first column is the name of the geneset. The second is the source for the geneset. Each further column lists genes within the geneset.')
    return(invisible(NULL))
  }
  
  sets <- list()
  raw  <- read.table( f, as.is=TRUE, sep='|', fill=TRUE, quote='' )$V1
  for ( i in 1:length(raw) ){
    tmp <- strsplit( raw[[i]], '\t' )[[1]]
    sets[[i]] <- sort(unique(as.character( setdiff( tmp[3:length(tmp)], '' ) )))
    names(sets)[i] <- tmp[1]
  }
  return(sets)
}

#' writeGMT
#'
#' Created for future use and simplification of writing .gmt files from the package.
#'
#' @param s list of gene sets in current R session. Each entry will become a row in the GMT file.
#' @param file file name to write to
#' @param link default is "". This is the second column of a GMT file and is usually a hyperlink or note about the origin of the term
#'
#' @return Nothing is returned. A file is written.
#' @export
#' 
#' @examples
#' # Make an example list() to show the GMT format
#' set <- list( term1=c('gene_name1','gene_name2'),
#'              term2=c('gene_name3','gene_name4','gene_name5') )
#' \dontrun{
#' # Write a GMT file for "set"
#' writeGMT( set, 'my_file.gmt')
#' }
writeGMT <- function( s, file = NA, link=rep('', length(s)) ){
  
  if(all(is.na(file))){
    stop('Please provide data and a file name:
         writeGMT(geneset_def, "myfile") where "geneset_def" is an object of list class that describes which genes are contained in a list of genesets.')
    return(invisible(NULL))
  }
  
  if (is.na(file)){ stop('no file name given') }
  if (class(s) != 'list'){ stop('input should be a list of genesets') }
  cat('', file=file, append=FALSE)
  for (x in 1:length(s) ){
    cat( sprintf('%s\t%s\t%s\n', names(s)[x], link[x],
                 paste( unique(s[[x]]), sep='\t', collapse='\t' ) ),
         file=file, append=TRUE )
  }
  return(NULL)
}


#' load_geneset_symbols
#'
#' For most applications, this function is used internally by term_enrichment(). Users may call this function directly in some cases to force FDR adjustment to be across multiple resources. See Vignette for more details.
#'
#' load_geneset_symbols allows the user to specify an annotation resource (e.g. Gene Ontology terms)
#' to use in enrichment analysis. The expectation is that the annotation resource contains of at least
#' one set of genes in the form of a list. The RITAN package comes with 15 pre-loaded annotation
#' resources. The default active annotation resources are GO, ReactomePathways,
#' KEGG_filtered_canonical_pathways, and MSigDB_Hallmarks.
#'
#' The result of calling this function is to set the variable "active_genesets" which will be used by further functions.
#'
#' @param gmt Either 1) name of pre-loaded resource (i.e. names(geneset_list)) or 2) gmt file containing annotation resources for enrichment annotation
#' @param gmt_dir location of gmt file named in gmt parameter
#' @param verbose print results to screen
#'
#' @return R list object named active_genesets
#' @export
#'
#' @examples
#' ## Load generic GO-slim terms
#' require(RITANdata)
#' load_geneset_symbols("GO_slim_generic")
#' print(length(active_genesets))
#' print(head(active_genesets[[1]]))
#' 
#' \dontrun{
#' ## load the default set of resources into "active_genesets"
#' load_geneset_symbols()
#' 
#' ## Use only the Reactome Pathways annotation resource.
#' load_geneset_symbols(gmt="ReactomePathways")
#' 
#' ## Suppresses output message describing the annotation resource and size.
#' load_geneset_symbols(gmt="ReactomePathways", verbose=FALSE)
#' 
#' ## To list the available resources within RITAN:
#' print(names(geneset_list))
#' 
#' ## You can also load your own data
#' load_geneset_symbols(gmt="myFile.gmt")
#' }
load_geneset_symbols <- function( gmt = NA, gmt_dir = '', verbose = TRUE ){
  sets <- NULL

  ## Cleanup inputs
  gmt_dir <- sub( '[/\\]+$', '', gmt_dir )

  if ( all(is.na(gmt)) ){
    ## ---------------------------------- -
    ## Default genesets
    sets <- do.call( c, geneset_list[ resources.default ])
    if (verbose){
      warning( sprintf('\nNo input annotation resource was selected. Defaulting to:\n%s',
                       paste(resources.default, sep=', ', collapse=', ') ))
    }

  } else if ( all(is.character(gmt)) && all( gmt %in% names(geneset_list) ) ){
    ## ---------------------------------- -
    ## Check if it's a known gene list
    sets <- do.call( c, geneset_list[ gmt ] )
    if (verbose){
      if ( all(nchar(gmt) < 100) ){
        cat(sprintf('Loading the requested genesets of "%s"...\n', gmt))
      } else {
        cat('Loading the requested genesets...\n')
      }
    }

  } else {

    ## ---------------------------------- -
    ## Loop over the requests and handle each.
    fc <- file.exists(gmt) | file.exists(sprintf('%s/%s', gmt_dir, gmt)) # file check
    nc <- gmt %in% names(geneset_list) # name check

    if (any(!( fc | nc ))){
      stop(sprintf('some resources were not recognized:\n%s', gmt[ !(fc|nc) ]))
    }

    res <- list()
    for (i in 1:length(gmt)){
      if (verbose){ cat(sprintf('Loading %s\n', gmt[i] )); }

      if (fc[i]){

        ## Load the indicated file
        this.f <- gmt[i]
        if (!file.exists(this.f)){
          this.f <- paste0( gmt_dir, '/', this.f )
        }
        res[[i]] <- readGMT(this.f)

      } else if (nc[i]){

        ## Load the known resource
        res[[i]] <- geneset_list[[ gmt[i] ]]

      }

    }

    sets <- do.call( c, res )

  }
  
  if ( (length(sets) == 0) && all(is.null(sets)) ){
    stop('No genesets loaded. Please check that your requested genesets are either:
         1) "%in% names(geneset_list)"
         2) are valid .gmt format files.'
         )
  }
  
  ## ---------------------------------- -
  ## load the gene set into global environment
  assign("active_genesets", sets, envir = .GlobalEnv)
  if (verbose){
    cat(sprintf( '\n\tLoaded %d genesets.\n', length(active_genesets) ))
  }
  
  # return() # omit explicit call to avoid the NULL print

}

#' show_active_genesets_hist
#'
#' function to plot distribution of size of active_genesets object
#'
#' @param nbins Number of bins to include in histogram
#' @param ... further argumants are passed on to plot()
#'
#' @return NULL. The plot is shown.
#' @export
#'
#' @examples
#' require(RITANdata)
#' load_geneset_symbols('GO_slim_generic')
#' show_active_genesets_hist()
#' 
#' \dontrun{
#' ## Show the distribution of geneset sizes for the default set of geneset resources
#' load_geneset_symbols()
#' show_active_genesets_hist()
#' 
#' ## Show the distribution of geneset sizes for a specific resource
#' load_geneset_symbols(gmt="ReactomePathways")
#' show_active_genesets_hist()
#' 
#' }
show_active_genesets_hist <- function(nbins = 50, ...){

  if (! exists("active_genesets")) {
    stop('The variable "active_genesets" is not set. See load_geneset_symbols().')
  }

  h <- hist( sapply(active_genesets,length), nbins, plot = FALSE )
  z <- h$density <= 0
  for (n in c('counts','density','mids')){
    h[[n]][z] <- NA # make NA so that it works with log()
  }

  plot( h$mids, h$counts, log='y', type = 'n',
        ylab = '# genesets', xlab = '# genes contained',
        main = 'Size of "active_genesets"', ... )

  for (i in 1 : (length(h$breaks)-1) ){
    rect( xleft = h$breaks[i], xright = h$breaks[i+1], ybottom = 0.1,
         ytop = h$counts[i], border = 'black', col='gray80' )
  }
  box()

}


### --------------------------------------------------------------- -
###
#' load_all_protein_coding_symbols
#' 
#' The character array returned is, by default, all human protein coding gene symbols. This variable defines the "universe of possible genes" for use in enrichment. Users should load a different "universe" or filter this one down to the most appropriate setting for their current study. For example, if running RNA-Seq, genes are in the universie if they are detected in any sample.
#' 
#' @param file file name of a table containing gene symbols
#' @param col_name column name within "file" that contains symbols
#' 
#' @return A unique list of gene symbols from the current protein coding set at the EBI
#' 
load_all_protein_coding_symbols <- function(
      file = 'ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_groups/protein-coding_gene.txt',
      col_name = 'symbol'
      ){
  s <- read.table( file, sep="\t", header=TRUE, as.is=TRUE, quote="", comment.char = '' )
  u <- unique( s[, col_name ])
  return(u)
}

### --------------------------------------------------------------- -
### perform enrichment on gene symbols
#' enrichment_symbols
#'
#' This function is called by term_enrichment() and term_enrichment_by_subset(). The user may call it directly, but we suggest using term_enrichment(). The function uses the resources currently loaded into the active_genesets vector. See load_geneset_symbols().
#'
#' Outputs a data frame containing the gene set name, a hypergeometric-test p value, the number of genes from the
#' input gene list that occur in the gene set, the number of genes in the gene set, the gene symbols for
#' the genes in the input gene list, and the q value.
#'
#' @param geneset vector of gene symbols to be evaluated
#' @param term a list containing specific gene set term(s) and their corresponding gene symbols contained in one of the annotation resources, default is all gene set terms
#'  in the GO, ReactomePathways, KEGG_filtered_canonical_pathways, and MSigDB_Hallmarks libraries
#' @param all_symbols gene symbols to be evaluated, identified by gene symbol name. Default is all protein coding genes. This parameter should be manipulated to include only the gene symbols that pertain to the user's analysis.
#' @param ... additional arguments are not used
#'
#' @return results matrix of input gene list compared to active gene sets. Q value is calculated using entire group of active gene sets.
#' @export
#'
#' @examples
#' require(RITANdata)
#' myGeneSet <- c('BRCA1','RAD51C','VAV1','HRAS','ABCC1','CYP1B1','CYP3A5')
#' 
#' \dontrun{
#' ## We suggest using term_enrichment() instead. E.g.:
#' e <- enrichment_symbols(myGeneSet, 'GO')
#' }
#' 
#' ## But, you may use enrichment_symbols() directly for an individual term:
#' load_geneset_symbols('GO')
#' e <- enrichment_symbols(myGeneSet, 'DNA_repair')
#' print(e)
#' 
#' \dontrun{
#' ## Gene set enrichment using intersection of gene symbols 
#' ##   provided in geneset parameter and all protein coding genes.
#' enrichment_symbols(geneset = vac1.day0vs31.de.genes)
#' 
#' ## choose which terms to evaluate
#' t <- active_genesets[1:5]
#' 
#' ## Test enrichment of that set of terms
#' enrichment_symbols(geneset = vac1.day0vs31.de.genes, term = t) 
#' }

enrichment_symbols <- function( geneset, term = NULL, all_symbols = NA, ... ){
    
    if (! exists("active_genesets")) {
      load_geneset_symbols()
    }

    if (all(is.na(all_symbols))){
      all_symbols <- load_all_protein_coding_symbols()
    }

    geneset <- geneset[ geneset %in% all_symbols ]
    K <- length(all_symbols) # length( setdiff( kit, c('.') ) ) # assayable
    k <- length(geneset) # our observation
    
    if (is.null(term)){

      ## Test all genesets
      active_int_term <- names(active_genesets)
      p <- array( NA, dim=c( length(active_genesets), 4 ) )
      for (i in 1:length(active_genesets)){

          this.set.possible <- (active_genesets[[i]])[ active_genesets[[i]] %in% all_symbols ]

          q <- sum( geneset %in% this.set.possible )
          p[i,2] <- q
          p[i,3] <- paste( geneset[ geneset %in% this.set.possible ], collapse=';' )
          p[i,4] <- length(unique( this.set.possible ))
          if ( q >= 2 ){
              m <- as.numeric(p[i,4])
              n <- K - m
              p[i,1] <- phyper( q-1, m, n, k, lower.tail = FALSE )
          } else {
              p[i,1] <- 1
          }

      }

    } else {
      
      if (is.null(names(term))){
        active_int_term <- intersect(names(active_genesets), term)
        if ( length(active_int_term) == 0 ){
          active_int_term <- names(active_genesets)[ grepl( sprintf('[[:alnum:]]*[.]%s$',term), 
                                                            names(active_genesets) ) |
                                                     grepl( sprintf('[[:alnum:]]*[.]%s$', gsub(' ','_',term) ), 
                                                            names(active_genesets) )
                                                    ]
        }
      } else {
        active_int_term <- intersect(names(active_genesets), names(term))
      }
      
      if ( length(active_int_term) < max(c( length(term), length(names(term)) )) ) {
        cat(sprintf('Undefined Term: %s', setdiff(term, active_int_term) ), sep = "\n")
        stop('You have requested terms that have not been defined. Please check contents of "term" parameter.')
      } else {
        
        ## Test the requested term
        p <- array( NA, dim=c( length(active_int_term), 4 ) )
        for (i in 1:length(active_int_term)){
          
          j <- which(names(active_genesets) == active_int_term[i])
          this.set.possible <- (active_genesets[[j]])[ active_genesets[[j]] %in% all_symbols ]
          q <- sum( geneset %in% this.set.possible )
          p[i,2] <- q
          p[i,3] <- paste( geneset[ geneset %in% this.set.possible ], collapse=';' )
          p[i,4] <- length(unique( this.set.possible ))
          if ( q >= 2 ){
            m <- as.numeric(p[i,4])
            n <- K - m
            p[i,1] <- phyper( q-1, m, n, k, lower.tail = FALSE )
          } else {
            p[i,1] <- 1 # if only one gene, set the p-value to 1 - not quite right, but singletons were leading to numeric errors...
          }
          
        }
      }
        
    }
    
    # tmp <- names(active_genesets)
    # if (!is.null(term)){
    #     tmp <- names(term)
    # }

    res <- data.frame( name  = active_int_term,
                       p     = as.numeric(p[,1]),
                       n     = as.numeric(p[,2]),
                       n.set = as.numeric(p[,4]),
                       genes = p[,3],
                       q     = p.adjust( as.numeric(p[,1]), method = "BH" )  )
    
    #if (dim(res)[1] < 100){
    #  warning('Few terms were considered here.
    #The q-values returned may not be appropriate for this particular analysis')
    #}

    res <- res[order(res$p), ]

    return(res)
}


term.test <- function( set, sub_set, total ){
    K <- length(total) # space of possibilities
    k <- length(set) # our observation size
    q <- sum( set %in% sub_set ) # observed hits
    m <- length(sub_set) # possible hits
    n <- K - m
    return( phyper( q-1, m, n, k, lower.tail = FALSE ) )
}

### --------------------------------------------------------------- -
#' geneset_overlap
#' 
#' Return assymetric matrix of the fraction of genes shared between sets. E.G. The fraction of the first set that is "covered" by or "overlaps" the second set.
#'
#' @param s1 The first geneset
#' @param s2 the second geneset
#' @param s.size Denominator used in each comparison. The default is to determint the lengths of elements in "s1"
#'
#' @return results matrix of input gene list compared to active gene sets. Q value is calculated using entire group of active gene sets.
#' @export
#'
#' @examples
#' require(RITANdata)
#' r <- geneset_overlap( geneset_list$MSigDB_Hallmarks, geneset_list$NetPath_Gene_regulation )
#' heatmap(r, col = rev(gray(seq(0,1,length.out = 15))) )
#' summary(c(r))
geneset_overlap <- function( s1, s2 = s1, #threshold = 0.5,
                             s.size = unlist(lapply( s1, length )) ){
    
    n <- length(s1)
    m <- length(s2)
    r <- array( 0, dim=c(n,m) )
    for (i in 1:n){

        for (j in 1:m ){
            r[i,j] <- sum( s1[[i]] %in% s2[[j]] ) / s.size[i] #length( s[[i]] )
        }

    }
    
    if (!all(is.null(names(s1)))){
      rownames(r) <- names(s1)
    }
    
    if (!all(is.null(names(s2)))){
      colnames(r) <- names(s2)
    }
    
    return(r)

}


#' resource_reduce
#' Merge terms across resources to reduce the number of redundant and semi-redundant terms
#'
#' @param genesets the input genesets to consider. May be from one or multiple resources.
#' @param min_overlap terms that share at least this fraction of genes will be merged
#' @param verbose if TRUE, print status and summary output 
#' 
#' @return the list of terms, after merging to reduce redundnat and semi-redundant terms
#' @export
#' 
#' @examples
#' require(RITANdata)
#' r <- resource_reduce( geneset_list$DisGeNet )
#'
resource_reduce <- function( genesets = NULL, min_overlap = 0.8, verbose = TRUE ){
  
  if (all(is.null(genesets))){
    stop('no terms/genesets were provided')
  }
  
  if ( (length(genesets) > 2000) & verbose ){
    warning(sprintf( 'This function is considering all pairwise relationships among the %d input terms/genesets.
This may take some time...', length(genesets) ))
  }
  
  o <- geneset_overlap(genesets)
  diag(o) <- NA # ignore self-overlap
  i <- which( o > min_overlap, arr.ind = TRUE )
  
  all_terms <- unique(c(rownames(o), colnames(o)))
  terms_with_overlap <- unique(c( rownames(o)[ i[,1] ], colnames(o)[ i[,2] ] ))
  
  merge_list <- new_genesets <- list()
  seen <- rep(FALSE, length(all_terms))
  gi <- 0
  for (r in 1:dim(i)[1]){
    
    ## initilize new group
    if (seen[ i[r,1] ] == FALSE){
      gi <- gi + 1
      seen[i[r,1]] <- seen[i[r,2]] <- TRUE
      merge_list[[gi]] <- c( rownames(o)[ i[r,1] ], colnames(o)[ i[r,2] ] )
    } else {
      next
    }
    
    ## populate group with all terms/sets that are highly related
    j <- (!seen) & ( ( all_terms %in% colnames(o)[ i[ i[,1] == i[r,1] , 2] ] ) |
                     ( all_terms %in% rownames(o)[ i[ i[,2] == i[r,2] , 1] ] ) )
    while( any(j) ){
      
      merge_list[[gi]] <- sort(unique(c( merge_list[[gi]], all_terms[j] )))
      seen[j] <- TRUE
      j <- (!seen) & ( ( all_terms %in% colnames(o)[ i[ i[,1] == i[r,1] , 2] ] ) |
                       ( all_terms %in% rownames(o)[ i[ i[,2] == i[r,2] , 1] ] ) )
      
    }
    
    ## generate the merged geneset
    new_genesets[[gi]] <- sort(unique(c(unlist( genesets[ merge_list[[gi]] ] ))))
    names(new_genesets)[gi] <- paste( merge_list[[gi]], sep=';', collapse = ';' )
    
  }
  
  out <- c( new_genesets, genesets[ !(names(genesets) %in% terms_with_overlap) ] )
  
  if (verbose){
    cat(sprintf('
The input list had %d terms/genesets.
The %d terms with overlaps of %.2f were merged into %d composite terms.
The updated term list with %d terms will be returned.
',
                length(genesets), length(terms_with_overlap), min_overlap, length(new_genesets), length(out) ))
  }
  
  return(out)
  
}



### --------------------------------------------------------------- -
### Enrichment analysis
#' term_enrichment
#'
#' term_enrichment evaluates the input gene list for enrichment within each of the annotation resources. This differs from the enrichment_symbols function which evaluates the gene list for enrichment against all of the annotation resources grouped together.
#'
#' @param geneset vector of gene symbols to be evaluated
#' @param resources list containing the reference gene sets to test for enrichment
#' @param report_resources_separately logical (default FALSE) flag to report enrichments seperately for each requested resource, or to combine them and produce FDR adjustment across the combined set
#' @param verbose print the top results for each annotation resource
#' @param all_symbols the background/global set of gene symbols (study dependent; we provide all protien coding genes as a default)
#' @param filter_to_intersection [FALSE] should the background and foreground genesets be subsetted to one another?
#' @param ... further arguments are passed on to enrichment_symbols()
#' 
#' @return results matrix of input gene list compared to active gene sets. Q value is calculated within each of the active gene sets.
#' @import RITANdata
#' @export
#'
#' @examples
#' ## Check if there is enrichment for any "Hallmark" functions within a input set of genes
#' require(RITANdata)
#' myGeneSet <- c('BRCA1','RAD51C','VAV1','HRAS','ABCC1','CYP1B1','CYP3A5')
#' e <- term_enrichment(myGeneSet, "MSigDB_Hallmarks")
#' print( e[1:2, -6] )
#' 
#' \dontrun{
#' term_enrichment(geneset = vac1.day0vs31.de.genes)
#' term_enrichment(geneset = vac1.day0vs31.de.genes, resources = "MSigDB_Hallmarks")
#' vac1.day0v31.enrichment <- term_enrichment(geneset = vac1.day0vs31.de.genes, verbose = FALSE)
#' }
term_enrichment <- function( geneset, resources = resources.default,
                             report_resources_separately = FALSE,
                             verbose = TRUE, all_symbols = NA,
                             filter_to_intersection = FALSE, ... ){
  
  if ( (!exists('all_symbols')) || all(is.na(all_symbols)) ){
    all_symbols <- load_all_protein_coding_symbols()
  }
  
  enrich <- list()
  process_source <- function( s_file, v=TRUE, f_intersect=FALSE ){
    
    load_geneset_symbols( s_file, verbose = v ) # load into active_genesets
    
    if (f_intersect){
      all_symbols_active <- unique(c(unlist( active_genesets ))) # baseline background
      ## -->> Need to detect if user gave "all_symbols" and take intersection <<--
      ## -->> Also, limit "geneset" to intersection <<--
      if (!is.null(get('all_symbols'))){
        all_symbols_active <- all_symbols_active[ all_symbols_active %in% all_symbols ]
      }
      geneset <- geneset[ geneset %in% all_symbols_active ]
    } else {
      all_symbols_active <- all_symbols
    }

    # run enrichment for active_genesets
    p0 <- enrichment_symbols( geneset, term = names(active_genesets),
                              all_symbols = all_symbols_active, ... )
    p0 <- p0[ order(p0$p), ]

    if (v){
      #printBanner( s_file )
      print( head( p0[ ,-5], n=5 ) )
    }

    return(p0)

  }

  if( report_resources_separately ){

    enrich <- lapply( resources, function(x){
                      process_source(x, v = verbose,
                                     f_intersect = filter_to_intersection) })
    names(enrich) <- resources

  } else {

    enrich <- process_source( resources, v = verbose,
                              f_intersect = filter_to_intersection )

  }
  
  class(enrich) <- c( 'term_enrichment', class(enrich) )
  attr(enrich, 'resources') <- resources
  
  return(enrich)

}

### --------------------------------------------------------------- -
### Enrichment analysis
#' term_enrichment_by_subset
#'
#' Run enrichment simultaneously across a group of prioritized gene lists. For example, in a time course dataset, one may have a different list of genes that are differentially expressed at each time point. This function facilitates rapid evaluation of term enrichment across time point comparisons. Alternatively, one may have a different list of differentially expressed genes by drug treatment, environmental condition, ect.
#'
#' @param groups A list() of genes for enrichment. Each entry in the list() is an input set of genes. Enrichment is performed for each of these entries.
#' @param resources character vector for which resources to use in enrichment
#' @param q_value_threshold minimum q-value (FDR adjusted p-value) in any group for the term to be included in results
#' @param verbose print additional status updates on what the function is doing
#' @param display_type Flag for which data type will be returned. One of "q" (default) for q-values, "p" for unadjusted p-values, or "n" for the number of genes overlapping the term.
#' @param phred Logical flag (default TRUE) to return the -log10 of p/q values
#' @param ... Further arguments are passed on to enrichment_symbols()
#'
#' @return Returns a term-by-study matrix of enrichment values (value determined by "display_type")
#' @import RITANdata
#' @export
#'
#' @examples
#' ## Create list of gene sets to evaluate.
#' ##   This example is from a vaccine study where we pre-generated differentially expressed genes.
#' ##   This object will be passed to the groups parameter.
#' require(RITANdata)
#' vac1.de.genes <- list(vac1.day0vs31.de.genes, vac1.day0vs56.de.genes)
#' names(vac1.de.genes) <- c("Day0vs31", "Day0vs56")
#' print(str(vac1.de.genes))
#' 
#' \dontrun{
#' ## Run term_enrichment_by_subset on the two results.
#' ##   This function usually takes a few seconds to a minute to run.
#' m <- term_enrichment_by_subset(groups = vac1.de.genes, q_value_threshold = .9)
#' summary(m)
#' plot( m, label_size_y = 4, show_values = FALSE )
#' }
#' 
term_enrichment_by_subset <- function( groups = NA, resources = resources.default,
                                       q_value_threshold = 0.01, verbose = TRUE,
                                       display_type = 'q', phred=TRUE, ... ){

  ## Make a matrix for the term enrichment of each "group" (a list where each entry is a vector of genes)
  ## "Extra" arguments are passed to enrichment_symbols()

  if ( all(is.na(groups)) || (class(groups) != 'list') ){
      stop( '"groups" should be a list where each entry is a vector of gene symbols.' )
  }
  if ((display_type == 'n') && ( phred == TRUE ) ){
    warning('"display_type" is set to "n" and "phred" is TRUE. We suggest setting phred=FALSE.\n')
  }
  
  ## Load the requested resources
  ## This will load them as a group into "active_genesets"
  load_geneset_symbols(resources, verbose = verbose)
  
  ## get enrichments (each is sorted by p)
  ## This function acts on "active_genesets"
  if(verbose){ cat("\tRunning Enrichment...", sep=" ") }
  e <- lapply( groups, enrichment_symbols, ... ) # optional arguments are sent on

  ## sort them together
  term.order <- e[[1]]$name
  e <- lapply( e, function(x){ x[ match(term.order,x$name), ] })

  ## filter by q-value
  if(verbose){ cat("FDR adjustment...", sep=" ") }
  i <- apply( do.call( cbind, lapply( e, function(x){ getElement(x,'q') <= q_value_threshold } ) ), 1, any )
  m <- cbind( e[[1]][ , c('name', 'n.set' )],
              do.call( cbind, lapply( e, function(x){ getElement(x, display_type) } ) )
              )
  first.data.col <- 3
  data.cols <- first.data.col : dim(m)[2]
  m <- m[ i, ]
  m <- m[ order(apply(m[, data.cols ], 1, min)), ]

  if( phred ){
      m[ , data.cols ] <- -log10(m[ , data.cols ])
  #} else {
  #    m <- m[ order( m[,first.data.col], decreasing=FALSE ), ]
  }

  if ( phred || (display_type == 'n')){
    m <- m[ order( apply(m[,data.cols], 1, max), decreasing=TRUE ), ]
  }

  m[ , data.cols ] <- round( m[ , data.cols ], 1 )
  class(m) <- c( 'term_enrichment_by_subset', class(m) )
  attr(m, 'display_type') <- display_type
  attr(m, 'phred') <- phred
  
  if(verbose){ cat("done.", sep="\n") }
  return( m )

}


#' plot.term_enrichment
#' @param x data frame returned by term_enrichment
#' @param min_q Only q-values more significant than this threshold will be plotted. Default = 0.05.
#' @param max_terms Up to max_terms will be plotted. Default = 25.
#' @param extend_mar Term names can be long. We attempt to keep them readable by extending the left-hand-side margins automatically. Default = c(0,10,0,0) added to par()$mar.
#' @param ... Additional arguments are passed on to plot()
#' @return silent return from plot
#' @import RITANdata graphics utils gplots png
#' @rawNamespace import(plotrix, except=c('plotCI'))
#' @rawNamespace import(stats, except=c('decompose','lowess','spectrum'))
#' @export
#' 
#' @examples
#' require(RITANdata)
#' e <- term_enrichment(vac1.day0vs31.de.genes, resources = 'GO_slim_generic')
#' plot(e, min_q = .1)
#'
plot.term_enrichment <- function( x = NA, min_q = 0.05, max_terms = 25, extend_mar = c(0,10,0,0), ... ){
  
  if ( is.null(dim(x)) || all(is.na(x)) ){
    stop('input data for plotting not found')
  }
  
  if (max_terms > dim(x)[1]){
    max_terms <- dim(x)[1]
  }
  
  i <- (x$q < min_q)
  if (!any(i)){
    stop('No terms meet the requested min_q threshold')
  }
  
  if ( length(extend_mar==4) & (par()$mar[2] < extend_mar[2]) ){
    ## The term labels can be long strings.
    ## Extent the left-hand-side margins.
    par( mar = par()$mar + extend_mar )
  }
  
  if (sum(i) > max_terms){
    i <- rep(FALSE, dim(x)[1])
    i[ which(x$q < min_q)[ 1:max_terms] ] <- TRUE
  }
  
  plot( -log10(x$q[i]), which(i),
        xlab = expression(paste('-log'[10], '( q-value )')),
        ylab = '', yaxt='n', ... )
  
  axis( side = 2, at = which(i), labels = rep('', sum(i)) )
  
  text( par()$usr[1] - 0.05*( par()$usr[2] - par()$usr[1] ), # x placement is 5% to the left of the x-axis starting
        which(i), labels = sub( '.', '\n', x$name[i], fixed=TRUE ),
        cex = 0.6, srt = 0, xpd=TRUE, adj=c(1.0,0.5) )
  
}


#' plot.term_enrichment_by_subset
#'
#' @param x data frame returned by term_enrichment_by_subset
#' @param show_values True or False, plot values on the heatmap
#' @param annotation_matrix a matrix() of group-levle characteristics - same number of columns as "m"
#' @param annotation_palates Color palates (RColorBrewer) used for each row of the annotation matrix
#' @param annotation_legend_x offset for placing the legend
#' @param low color for low end of range
#' @param high color for high end of range
#' @param label_size_x size of text for x label. Default lable_size_x=16
#' @param label_angle_x angle for text for x label. Default is -30 degrees
#' @param label_size_y size of text for y label. Default label_size_y=9
#' @param grid_line_color color o grid lines between cells. Default is white.
#' @param mid sets lower threshold for color scale
#' @param wrap_y_labels Number of characters to wrap row labels
#' @param cap Clip numeric values to this maximum threshold
#' @param return_ggplot_object logical flag (default FALSE) that if TRUE, the ggplot object for the plot is returned
#' @param ... further areguments are not used at this time. If the user wants to modify the plot, use return_ggplot_object = TRUE.
#'
#' @return silent return, unless return_ggplot_object==TRUE. Then, the ggplot object for the plot is returned.
#' @export
#' @import ggplot2 reshape2 grid gridExtra RColorBrewer 
#' @examples
#' ## Create list of gene sets to evaluate.
#' ##   This example is from a vaccine study where we pre-generated differentially expressed genes.
#' ##   This object will be passed to the groups parameter.
#' require(RITANdata)
#' vac1.de.genes <- list(vac1.day0vs31.de.genes, vac1.day0vs56.de.genes)
#' names(vac1.de.genes) <- c("Day0vs31", "Day0vs56")
#' print(str(vac1.de.genes))
#' 
#' \dontrun{
#' ## Run term_enrichment_by_subset on the two results.
#' ##   This function usually takes a few seconds to a minute to run.
#' m <- term_enrichment_by_subset(groups = vac1.de.genes, q_value_threshold = .9)
#' summary(m)
#' plot( m, label_size_y = 4, show_values = FALSE )
#' }
#'
plot.term_enrichment_by_subset <- function( x, show_values = TRUE,
                                            annotation_matrix = NA,
                                            low = 'white', high="#2166AC",
                                            return_ggplot_object = FALSE,
                                            label_size_x = 16, label_angle_x=-30,
                                            label_size_y = 9, wrap_y_labels = 20,
                                            grid_line_color = 'white', mid = 0, cap = NA,
                                            annotation_palates = c('Reds','Greens','Purples','Greys','BuPu','RdPu','BrBG','PiYG','Spectral'),
                                            annotation_legend_x = -0.3, ... ){

  ## Check inputs
  if ( all(is.na(x)|is.null(x)) || is.null(dim(x)) || (dim(x)[1] < 1) ){
    stop(' insufficient data provided to plot.term_enrichment_by_subset()')
  }
  display_type <- attr(x, 'display_type')
  phred <- attr(x, 'phred')

  if ( (!all(is.na(annotation_matrix))) &&
       ((dim(x)[2] - dim(annotation_matrix)[2]) != 2 ) ){
    stop(sprintf('You have requested to add an annotation matrix, but the dimentions do not line up.\nExpected: %d',
                 dim(x)[2]-2, dim(annotation_matrix)[2] ))
  }

  ## Optionally, cap the values for plotting
  if ( (!all(is.na(cap))) & (class(cap) == 'numeric') ){
    if (length(cap) > 1){ stop('"cap" should be a single number') }
    x[, 3:dim(x)[2] ] <- apply( x[, 3:dim(x)[2] ], c(1,2), function(y){
      if (y > cap){ return(cap)
      } else {      return(y)
      }
    })
  }

  ## make matrix with row/column names
  mat <- as.matrix( x[ , 3:dim(x)[2] ] )
  colnames(mat) <- colnames(x)[ 3:dim(x)[2] ]

  ## trim off the "resource" from term names for plotting
  rownames(mat) <- as.character(sapply(
    getElement( x, 'name' ), function(x){
      sub( '.+[.]', '', x )
    }))

  ## basic 1-level name wrapping for plotting
  if ( !all(is.na(wrap_y_labels)) && is.numeric(wrap_y_labels) ){
    rownames(mat) <- as.character( sapply( rownames(mat), function(x){
        if( nchar(x) > wrap_y_labels ){
          y <- strsplit( x, '' )[[1]]
          i <- grep( '[-_ ]', y )
          j <- wrap_y_labels
          if(any(i)){
            j <- i[ which.min( abs(i-wrap_y_labels) ) ]
          }
          x <- paste( c(y[1:j], '\n', y[ (j+1) : length(y) ]), sep='', collapse = '' )
        }
        return(x)
      }))
  }


  ## make ggplot2 compatible
  #dat <- reshape2::melt( mat )
  #levels(dat$Var1) <- rownames(mat)
  #levels(dat$Var2) <- colnames(mat)
  
  ## make ggplot2 compatible - updated
  dat <- reshape2::melt( mat, varnames=c('Var1', 'Var2') )
  levels(dat$Var1) <- rownames(mat)
  levels(dat$Var2) <- colnames(mat)
  
  ## It would be nice to add a >= sign to the legend...
  #if ( show_values & (!all(is.na(cap))) & (class(cap) == 'numeric') ){
  #  dat$value[ dat$value == cap ] <- sprintf( '\u2265%f', dat$value )
  #}

  ## make ggplot heatmap
  g <- ggplot( dat, aes(Var2,Var1)) +
    geom_tile( aes_string(fill = 'value' ), colour=grid_line_color )
  if (show_values){
     g <- g + geom_text( label = round(dat$value, 1) ) # aes_string(fill = 'value' )
  }
  g <- g + scale_fill_gradient2(midpoint=mid, low=low, high=high)
  g <- g + theme( axis.ticks = element_blank(),
                  axis.text.x = element_text(size = label_size_x, angle = label_angle_x,
                                             colour = "black", hjust = 0),
                  axis.text.y = element_text(size = label_size_y,
                                             colour = "black"),
                  panel.background = element_blank()
  ) + xlab('') + ylab('')


  ## Check if the user has supplied a sample/study/column annotation matrix
  if (!all(is.na(annotation_matrix))){

    mc <- max(apply(annotation_matrix, 1, function(x){ length(unique(x)) }))
    if (mc > 8){
      stop('Generation of the annotation_matrix visualization currently uses RColorBrewer.
These color palates are only defined for up to 9 colors (including white).
Please modify your annotation_matrix to have up to 8 values per row.')
    }

    ## set color palates to use. Repeat Spectral as necessary
    nr <- dim(annotation_matrix)[1]
    if ( nr > length(annotation_palates)){
      d <- nr - length(annotation_palates)
      annotation_palates <- c(annotation_palates, rep('Spectral', d) )
    }

    ## Translate data into color values
    dat <- reshape2::melt(annotation_matrix); names(dat) <- c('Var1','Var2','value')
    levels(dat$Var1) <- rownames(annotation_matrix)
    levels(dat$Var2) <- colnames(annotation_matrix)
    dat$col <- dat$colID <- NA
    uc <- ul <- NULL
    ci <- 0
    for (i in 1:nr){

      nc <- length(unique(annotation_matrix[i,]))
      row_col <- brewer.pal( n = nc+1, name = annotation_palates[i] )[-1]
      uc <- c(uc, row_col)
      ul <- c(ul, unique(annotation_matrix[i,]) )
      for (j in 1:nc){
        ci <- ci+1

        k <- (dat$Var1 == rownames(annotation_matrix)[i]) &
             (dat$value == unique(annotation_matrix[i,])[j] )
        dat$col[k] <- row_col[j]
        dat$colID[k] <- ci
      }

    }
    dat$colID <- factor(dat$colID)

    ## Make ggplot
    ggmat <- ggplot( dat, aes(Var2,Var1)) +
      geom_tile( aes_string(fill = 'colID' ), colour=grid_line_color ) +
      scale_fill_manual( values = uc, labels = ul ) +
      theme( axis.ticks = element_blank(),
             axis.text.y = element_blank(),
             axis.text.x = element_text(size = 8, angle = -45, hjust = 0),
             panel.background = element_blank(),
             legend.spacing = unit(0, "cm"), legend.text = element_text(size=rel(0.5)),
             title = element_blank(), legend.key.size = unit(.3,'cm'),
             legend.position = c(annotation_legend_x, 0.5)
      ) + xlab('') + ylab('')

    ## Combine with heatmap
    p1 <- ggplotGrob(g)
    p2 <- ggplotGrob(ggmat)
    p2$widths <- p1$widths
    g <- gridExtra::arrangeGrob( p2, p1, heights = c(1,4), ncol=1 )
    
    ## Return or Plot
    if (return_ggplot_object){
      #g <- arrangeGrob( p2, p1, layout_matrix = matrix(c(1, rep(2,4) ), ncol=1) )
      return(g)
    } else {
      plot(g)
    }

  } else {

    ## Return or Plot
    if (return_ggplot_object){
      return(g)
    } else {
      plot(g)
    }

  }

}


#' summary.term_enrichment
#'
#' @param object data frame returned by term_enrichment()
#' @param ... Further arguments are passed on to head()
#' 
#' @return the data.frame of top enrichment results
#' @export
#' 
#' @examples
#' require(RITANdata)
#' e <- term_enrichment( vac1.day0vs31.de.genes, "MSigDB_Hallmarks" )
#' summary(e, n=3)
#' 
summary.term_enrichment <- function(object, ...){
  
  if ( 'list' %in% class(object) ){
    lapply( object, function(y){
      print(head(y[,-5], ...))
    })
  } else {
    print(head(object[,-5], ...))
  }
  
}


#' summary.term_enrichment_by_subset
#'
#' @param object data frame returned by term_enrichment_by_subset()
#' @param verbose if TRUE (default), print a header describing the data type
#' @param ... Further arguments are passed on to head()
#' 
#' @return the data.frame of top enrichment results
#' @export
#' 
#' @examples
#' require(RITANdata)
#' vac1.de.genes <- list(vac1.day0vs31.de.genes, vac1.day0vs56.de.genes)
#' names(vac1.de.genes) <- c("Day0vs31", "Day0vs56")
#' e <- term_enrichment_by_subset(vac1.de.genes, "MSigDB_Hallmarks", q_value_threshold = 0.1 )
#' summary(e)
#'
summary.term_enrichment_by_subset <- function(object, verbose = TRUE, ...){
  
  if (verbose){
    # attributes(object)
    
    dt <- attr(object, 'display_type')
    ph <- attr(object, 'phred')
    
    if (dt == 'p'){
      
      ph_tag <- ifelse( ph, '-log10(p)', 'p' )
      print(sprintf('The data matrix contains %s-values which have not been adjusted for multiple testing', ph_tag))
      
    } else if (dt == 'q'){
      
      ph_tag <- ifelse( ph, '-log10(q)', 'q' )
      print(sprintf('The data matrix contains %s-values which have been adjusted for multiple testing (FDR)', ph_tag))
      
    } else if (dt == 'n'){
      
      print('The data matrix contains the number of genes in the query that are contained within each subset')
      
    } else {
      print( sprintf('Data type "%s" not recognized.', as.character(dt)) )
    }
    
    print('Each term has its own row and each subset/substudy is a column')
    
  }
  
  head(object, ...)
  
}

