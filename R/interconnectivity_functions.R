
# ss2 <- function(x, ...){
#   s <- network_overlap( x, Net2Use = c('ISHQ','HPRD','CCSB','HumanNet','STRING'),
#                minStringScore = 700, # 7.8% have a score >= 0.7
#                minHumanNetScore = 0.0, # 0.5 ==> 84.8%, 1.0 ==> 44.7%
#                dedup = TRUE, directed.net = FALSE, ... ) 
#   return(s)
# }


#' icon_test
#' 
#' "icon" is an abbreviation for the "interconnectivity" of a network or graph.
#' 
#' This function handles different inputs and directs them to the appropriate "icon" testing method. Depending on the values given to "n1" and "n2," a different specific test is performed.
#' 
#' Note that the specific functions called make use of the "param" attribute of each input. These parameters are populated by network_overlap() so that the permutation reflects the exact procedure that was done to generate "n1" and/or "n2."
#' 
#' @param n1 [NULL] the first network. See network_overlap().
#' @param n2 [NULL] the second network. See network_overlap().
#' @param s [100] teh number of random permutations to make.
#' @param verbose [TRUE] If optional text describing what the algorithm is doing should be shown in the console.
#' @param ... Additional argumetns are passed on to the specific test performed
#' 
#' @export
#' @return metrics and significance of the network overlap
#' 
#' @examples
#' \dontrun{
#' icon_test( n1=n, s=10)
#' }
#' 
icon_test <- function( n1 = NULL, n2 = NULL, s=100, verbose=TRUE, ... ){
  
  stop('TO DO: update all tests to use the cached parameters from each network_overlap() result')
  param1 <- attr(n1,'param')
  
  
  
  
  test_mode <- 'dual_between'
  if (all(is.null(c(n1,n2)))){
    stop('input required')
  }
  #print(class(n1))

  o <- NULL
  if ( all(is.null(n2)) || all(is.na(n2)) ){

    if (class(n1) == 'character'){
      test_mode <- 'single_within'
      o <- icon_single_within( n1, s )
    } else if (class(n1) == 'list'){
      test_mode <- 'list_between'
      o <- icon_list_between( n1, s, verbose=verbose, ... )
    }

  } else {
    o <- icon_dual_between( n1, n2, s, verbose=verbose, ... )
  }

  o$test_mode <- test_mode
  return(o)

}

cov_undirected <- function( this_n1, this_n2, this_net1, this_net2 ){
  cov12 <- ( this_n1 %in% this_net2$p1 ) | ( this_n1 %in% this_net2$p2 )
  cov21 <- ( this_n2 %in% this_net1$p1 ) | ( this_n2 %in% this_net1$p2 )
  gin <- unique(c( this_n1[ cov12 ], this_n2[ cov21 ] ))
  ag  <- unique(c( this_n1, this_n2 ))
  obs <- length( gin ) / length(ag)
  return(obs)
}

icon_single_within <- function( n1, s, verbose=TRUE, net_within = NULL ){

  if (all(is.null(net_within))){
    net_within <- ss2( n1, verbose = verbose )
  }

  int_within <- (n1 %in% net_within$p1) | (n1 %in% net_within$p2)
  obs <- sum( int_within ) / length(n1)

  ran <- rep(NA,s)
  if (verbose){ pb <- txtProgressBar( style = 3 ); }
  for ( k in 1:s ){
    if (verbose){ setTxtProgressBar(pb = pb, value = k/s) }

    r <- sample( all_symbols, length(n1), replace=FALSE )
    #n_ran <- ss2(r, verbose=FALSE)
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

icon_dual_between <- function( n1, n2, s, verbose=TRUE, net1=NULL, net2=NULL ){

  if(all(is.null(net1))){
    net1 <- ss2( n1, verbose = verbose, include_neighbors = TRUE )
  }
  if(all(is.null(net2))){
    net2 <- ss2( n2, verbose = verbose, include_neighbors = TRUE )
  }
  #print(str(net1))
  #print(str(net2))

  obs <- cov_undirected( n1, n2, net1, net2 )
  #print(obs)

  ran <- rep(NA,s)
  if (verbose){ pb <- txtProgressBar( style = 3 ); }
  for ( k in 1:s ){
    if (verbose){ setTxtProgressBar(pb = pb, value = k/s) }
    r1 <- sample( all_symbols, length(n1), replace=FALSE )
    r2 <- sample( all_symbols, length(n2), replace=FALSE )
    n1_ran <- all_net[ (all_net$p1 %in% r1) | (all_net$p2 %in% r1), ]
    n2_ran <- all_net[ (all_net$p1 %in% r2) | (all_net$p2 %in% r2), ]
    ran[k] <- cov_undirected( r1, r2, n1_ran, n2_ran )
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



