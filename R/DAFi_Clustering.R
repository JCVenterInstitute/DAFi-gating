#' Clustering algorithm wrapper function
#'
#' Will generate the clustering results necessary for downstream filtering/gating and plotting
#'
#' @param fcs_in    compensated and transformed FCS in the flowFrame format as input
#' @param colsToUse select the channels to use for the clustering: c(7,8,9,10)
#' @param clusterSize   select the number of desired clusters for the Kmeans-type clustering algorithms (default 100)
#' @param xdims         number of SOM nodes in x direction for FlowSOM algorithm (default 10)
#' @param ydims         number of SOM nodes in y direction for FlowSOM algorithm (default 10)
#' @param method        set clustering algorithm from the following choices: 'Kmeans' (default), 'fSOM'
#' @param initializer   set initializer for Kmeans methods from the following choices: 'optimal_init', 'quantile_init', 'kmeans++' (default), and 'random'; see \href{https://cran.r-project.org/web/packages/ClusterR/ClusterR.pdf}{ClusterR}
#'
#' @return A \code{list} with four items: the first is the clustering result object directly from ClusterR or FlowSOM 
#'         containing all information as originally constructed in the respective clustering packages, the second is the cluster membership mapping 
#'         according the the original sequence from the input flowFrame, the third is the centroids in terms of the medianvalues of the columns 
#'         from the flowFrame based on the membership mapping, and the fourth is the new flowFrame containing just the k (or xdims X ydims) number of centroids.
#'         This is a wrapper function for clustering methods included in
#'         \href{https://cran.r-project.org/web/packages/ClusterR/ClusterR.pdf}{ClusterR} and
#'         \href{https://bioconductor.org/packages/release/bioc/html/FlowSOM.html}{FlowSOM}
#'         Please see their documentations for more options.
#'
#' @seealso \code{\link{DAFi}}
#' 
#' @examples
#' # Kmeans with random initializer
#' clusterResult <- cluster_run(sampleflowFrame, colsToUse=c(9,12,14:18), clusterSize=500,
#'                       method="Kmeans", initializer="random")
#'                       
#' # Kmeans with kmeans++ initializer
#' clusterResult <- cluster_run(sampleflowFrame, colsToUse=c(9,12,14:18), clusterSize=500)
#' 
#' # FlowSOM
#' clusterResult <- cluster_run(sampleflowFrame, colsToUse=c(9,12,14:18), xdims=20, ydims=25, method='fSOM')
#' @export
cluster_run <-
  function(fcs_in,
           colsToUse,
           clusterSize = 100,
           xdims = 10,
           ydims = 10,
           method = "Kmeans",
           initializer = "kmeans++",
           debugmode = FALSE) {
    
    cat("Clustering via", method, "\n")
    
    clusterResult <-
      list(
        Object = NULL,
        mapping = NULL,
        centroids = NULL,
        flowFrame = NULL
      )
    if (method == "fSOM") {
      clusterResult$Object <-
        FlowSOM(
          fcs_in,
          compensate = FALSE,
          transform = FALSE,
          scale = FALSE,
          colsToUse = colsToUse,
          xdim = xdims,
          ydim = ydims,
          nClus = 10,
          seed = 42
        )
      if(isTRUE(debugmode)){print("create mapping")}
      clusterResult$mapping <-
        clusterResult$Object$FlowSOM$map$mapping[, 1]
      if(isTRUE(debugmode)){print("create centroids")}
      clusterResult$centroids <-
        clusterResult$Object$FlowSOM$map$medianValues
      if(isTRUE(debugmode)){print("create cluster flowFrame")}
      clusterResult$flowFrame <- flowFrame(clusterResult$centroids)
    } else if (method == "Kmeans") {
      data = exprs(fcs_in)
      clusterResult$Object <-
        KMeans_rcpp(
          data[, colsToUse],
          clusters = clusterSize,
          num_init = 1,
          max_iters = 100,
          initializer = initializer,
          fuzzy = FALSE,
          verbose = FALSE,
          CENTROIDS = NULL
        )
      clusterResult$mapping <- clusterResult$Object$clusters
      clusterResult$centroids <- t(sapply(seq(clusterSize), function(i) {
        apply(subset(data, clusterResult$mapping == i),
              2,
              stats::median)
      }))
      clusterResult$flowFrame <- flowFrame(clusterResult$centroids)
    }
    clusterResult
  }