#' icon_test
#' 
#' "icon" is an abbreviation for the "interconnectivity" of a network or graph.
#' 
#' This function handles different inputs and directs them to the appropriate "icon" testing method. Depending on the values given to "nodes1" and "nodes2," a different specific test is performed.
#' 
#' Note that the specific functions called make use of the "param" attribute of each input. These parameters are populated by network_overlap() so that the permutation reflects the exact procedure that was done to generate "nodes1" and/or "nodes2."
#' 
#' @param nodes1 [NULL] the first network. See network_overlap().
#' @param nodes2 [NULL] the second network. See network_overlap().
#' @param s [100] the number of random permutations to make.
#' @param verbose [TRUE] If optional text describing what the algorithm is doing should be shown in the console.
#' @param ... Additional argumetns are passed on to the specific test performed
#' 
#' @export
#' @return metrics and significance of the network overlap
#' 
#' @examples
#' \dontrun{
#' icon_test( nodes1=n, s=10)
#' }
#' 
icon_test <- function( nodes1 = NULL, nodes2 = NULL, s=100, verbose=TRUE, ... ){
  
  stop('TO DO: update all tests to use the cached parameters from each network_overlap() result')
  param1 <- attr(nodes1,'param')
  
  test_mode <- 'dual_between'
  if (all(is.null(c(nodes1,nodes2)))){
    stop('input required')
  }
  #print(class(nodes1))

  o <- NULL
  if ( all(is.null(nodes2)) || all(is.na(nodes2)) ){

    if (class(nodes1) == 'character'){
      test_mode <- 'single_within'
      o <- icon_single_within( nodes1, s )
    } else if (class(nodes1) == 'list'){
      test_mode <- 'list_between'
      o <- icon_list_between( nodes1, s, verbose=verbose, ... )
    }

  } else {
    o <- icon_dual_between( nodes1, nodes2, s, verbose=verbose, ... )
  }

  o$test_mode <- test_mode
  return(o)

}

#' cov_undirected
#' function to show the un-directed coverabe between two nodes lists, from two networks
#' @param this_nodes1 list of nodes for first network
#' @param this_nodes2 list of nodes for second network
#' @param this_net1 the first network
#' @param this_net2 the second network
cov_undirected <- function( this_nodes1, this_nodes2, this_net1, this_net2 ){
  cov12 <- ( this_nodes1 %in% this_net2$p1 ) | ( this_nodes1 %in% this_net2$p2 )
  cov21 <- ( this_nodes2 %in% this_net1$p1 ) | ( this_nodes2 %in% this_net1$p2 )
  gin <- unique(c( this_nodes1[ cov12 ], this_nodes2[ cov21 ] ))
  ag  <- unique(c( this_nodes1, this_nodes2 ))
  obs <- length( gin ) / length(ag)
  return(obs)
}

#' icon_single_within
#' interconnectivity score within a network
#' @param nodes the node labels to use
#' @param net the network to use
#' @param s [10] the number of repeated random draws to make
#' @param verbose [TRUE] if more verbose output should be shown
icon_single_within <- function( nodes = NULL, net = NULL, s=10, verbose=TRUE ){

  if (all(is.null(net))){
    stop('you must input a network to search within')
  }
  if (all(is.null(nodes))){
    stop('you must input a list of node labels to search for, within the network')
  }

  int_within <- (nodes %in% net$p1) | (nodes %in% net$p2)
  obs <- sum( int_within ) / length(nodes)

  ran <- rep(NA,s)
  if (verbose){ pb <- txtProgressBar( style = 3 ); }
  for ( k in 1:s ){
    if (verbose){ setTxtProgressBar(pb = pb, value = k/s) }

    r <- sample( all_symbols, length(nodes), replace=FALSE )
    n_ran <- all_net[ (all_net$p1 %in% r) & (all_net$p2 %in% r), ]
    i_ran <- (r %in% n_ran$p1) | (r %in% n_ran$p2)
    ran[k] <- sum( i_ran ) / length(r)
  }

  out <- list()
  out$obs <- obs
  out$ran <- ran
  out$p <- sum(ran) >= obs
  out$z <- (obs - median(ran)) / mad(ran, constant = 1)

  if (verbose){
    d <- density(ran)
    d$y <- d$y / sum(d$y)
    plot( 0, type='n', xlim=range(c(obs,ran)), ylim=range(d$y),
          main = sprintf('I = %.2f. Z = %.1f', obs, out$z ),
          xlab = 'Interconnectivity', ylab = 'Probability Density' )
    polygon(d, col = 'gray')
    abline( v=obs, col='red')
    box()
  }

  return(out)

}

icon_dual_between <- function( nodes1=NULL, nodes2=NULL, s,
                               verbose=TRUE, net1=NULL, net2=NULL ){

  if ( all(is.null(net1)) || all(is.null(net2)) ){
    stop('you must input 2 networks to search between')
  }
  
  if ( all(is.null(nodes1)) || all(is.null(nodes2)) ){
    stop('you must input two lists of node labels to search for, between the two networks')
  }
  
  obs <- cov_undirected( nodes1, nodes2, net1, net2 )
  
  ran <- rep(NA,s)
  if (verbose){ pb <- txtProgressBar( style = 3 ); }
  for ( k in 1:s ){
    if (verbose){ setTxtProgressBar(pb = pb, value = k/s) }
    r1 <- sample( all_symbols, length(nodes1), replace=FALSE )
    r2 <- sample( all_symbols, length(nodes2), replace=FALSE )
    nodes1_ran <- all_net[ (all_net$p1 %in% r1) | (all_net$p2 %in% r1), ]
    nodes2_ran <- all_net[ (all_net$p1 %in% r2) | (all_net$p2 %in% r2), ]
    ran[k] <- cov_undirected( r1, r2, nodes1_ran, nodes2_ran )
  }

  out <- list()
  out$obs <- obs
  out$ran <- ran
  out$p <- sum(ran) >= obs
  out$z <- (obs - median(ran)) / mad(ran, constant = 1)

  if (verbose){
    d <- density(ran)
    d$y <- d$y / sum(d$y)
    plot( 0, type='n', xlim=range(c(obs,ran)), ylim=range(d$y),
          main = sprintf('I = %.2f. Z = %.1f', obs, out$z ),
          xlab = 'Interconnectivity', ylab = 'Probability Density' )
    polygon(d, col = 'gray')
    abline( v=obs, col='red')
    box()
  }

  return(out)

}


icon_list_between <- function( n_list, s, verbose=TRUE, net_list = NULL ){

  preCompLists <- FALSE
  if ( all(!is.null(net_list)) & (length(n_list) == length(net_list)) ){
    preCompLists <- TRUE
  }
  #cat(sprintf('Using pre-compiled lists? %s', preCompLists))

  m <- length(n_list)
  z_mat <- array( 0, dim=c(m,m) )
  for (i in 1:m){
    if(verbose){
      cat(sprintf('\nCalculating interactions with input set %d.',i))
    }

    ## NOTE: the "within" lists shoudl not have neighbors, but the "between" should...
    #if(preCompLists){
    #  z_mat[i,i] <- icon_single_within( n_list[[i]], s=s, verbose = verbose,
    #                                    net_within = net_list[[i]] )$z
    #} else {
      z_mat[i,i] <- icon_single_within( n_list[[i]], s=s, verbose = verbose )$z
    #}

    if ( (i+1) <= m ){
    for (j in (i+1):m){

      if(preCompLists){
        z_mat[i,j] <- z_mat[j,i] <-
          icon_dual_between( n_list[[i]], n_list[[j]], s=s, verbose=verbose,
                             net_list[[i]], net_list[[j]] )$z
      } else {
        z_mat[i,j] <- z_mat[j,i] <-
          icon_dual_between( n_list[[i]], n_list[[j]], s=s, verbose=verbose )$z
      }

    }}
  }

  out <- list()
  out$z_mat <- z_mat
  return(out)

}



