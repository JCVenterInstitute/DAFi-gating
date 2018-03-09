#' Run the DAFi algorithm
#'
#' Method to run general DAFi workflow.
#' Will transform the data via FCSTrans and scale each channel to the range 0-4096 by default.
#'
#' @param rawfcs_path   path to raw fcs file
#' @param config_path   path to DAFi cell gating hierarchy configuration file (inclusion gates)
#' @param reverse_path  path to DAFi cell gating hierarchy configuration file (exclusion gates)
#' @param k             specify number of clusters for Kmeans type algorithms (default 100)
#' @param xdims         number of SOM nodes in x direction for FlowSOM algorithm (default 10)
#' @param ydims         number of SOM nodes in y direction for FlowSOM algorithm (default 10)
#' @param excludeGated  *experimental* set TRUE to exclude previously gated channels in reclustering (default FALSE)
#' @param method        set clustering algorithm from the following choices: 'Kmeans' (default), 'fSOM'
#' @param initializer   set initializer for Kmeans methods from the following choices: 'optimal_init', 'quantile_init', 'kmeans++' (default), and 'random'; see \href{https://cran.r-project.org/web/packages/ClusterR/ClusterR.pdf}{ClusterR}
#' @param plotCentroids whether or not to plot diagnostic plots to see centroid locations
#' @param writeOutput   whether or not to generate output events and percentage tables
#'
#' @seealso \code{\link{cluster_run}},\code{\link{scatterTransform}},\code{\link{timeTransform}}
#'
#' @examples
#' # Minimal parameter run (all defaults)
#' DAFi(SampleData,inclusionConfig)
#'
#' # Or specify method and options for FlowSOM clustering
#' DAFi(SampleData,inclusionConfig,exclusionConfig, xdims=20, ydims=25, method='fSOM')
#'
#' #' # Or specify method and options for Kmeans++ clustering
#' DAFi(SampleData,inclusionConfig",exclusionConfig", k=500, method='Kmeans', initializer='kmeans++')
#'
#'
#' @import methods
#' @import utils
#' @import flowCore
#' @import flowViz
#' @import FlowSOM
#' @import ClusterR
#' @importFrom flowStats curv1Filter
#'
#' @export
DAFi <-
  function(rawfcs_path,
           config_path,
           reverse_path = NULL,
           k = 100,
           xdims = 10,
           ydims = 10,
           excludeGated = FALSE,
           method = "Kmeans",
           initializer = "kmeans++",
           plotCentroids = TRUE,
           writeOutput = TRUE,
           smoothPlot = FALSE,
           debugmode = FALSE) {
    #Parse dafi config file
    config_table <- read.table(config_path,
                               sep = "\t", header = FALSE)
    colnames(config_table) <-
      c(
        "PopIndex",
        "ChannelX",
        "ChannelY",
        "minX",
        "maxX",
        "minY",
        "maxY",
        "ParentIndex",
        "FilterType",
        "CreateSubpop",
        "Recluster",
        "Phenotype"
      )
    cat("Parsed configuration table\n")
    
    if (!is.null(reverse_path)) {
      reverse_table <- read.table(reverse_path,
                                  sep = "\t", header = FALSE)
      colnames(reverse_table) <-
        c("PopIndex",
          "ChannelX",
          "ChannelY",
          "minX",
          "maxX",
          "minY",
          "maxY",
          "ParentIndex")
      cat("Parsed reverse configuration table\n")
    }
    
    #read in, compensate, transform, and format fcs file
    fcsTransInput <- FCSTransInput(rawfcs_path, debugmode = debugmode)
    
    #fcs_raw = fcsTransInput$rawflowFrame
    colsToUse = fcsTransInput$colsToUse
    fcs = fcsTransInput$processedFlowFrame
    sampleName = tools::file_path_sans_ext(basename(fcs@description$FILENAME))
    
    #Locate columns (marker channels) used for clustering
    colsToUseList <- list()
    colsToUseList[[1]] <- colsToUse
    
    #Setup file connections for output
    percentfileConn <- file("pop_percentage.txt")
    writeLines(c(
      paste(
        "Predefined_Population_ID",
        "PercentageBasedonParent",
        sep = "\t"
      )
    ), percentfileConn)
    eventsfileConn <- file("pop_events.txt")
    writeLines(c(paste(
      "Predefined_Population_ID", "NumofEvents", sep = "\t"
    )), eventsfileConn)
    close(percentfileConn)
    close(eventsfileConn)
    
    #setup data vectors to store intermediate results
    popEvent <- list()
    popPercent <- list()
    pop_filt_list <- list()
    pop_cluster_list <- list()
    popGatingMatrixList <- list()
    reverse_popGatingMatrixList <- list()
    popGatingFilterList <- list()
    popFlowFrameList <- list()
    
    #events clustering info
    clusterResults <- list()
    
    finalPopMatrix <- matrix(0, nrow(fcs), 1)
    colnames(finalPopMatrix) <- c("Population")
    
    ##############DAFi main iterative filtering-gating loop (loop through each line of the configuration file to find cells of each population)
    for (i in 1:nrow(config_table)) {
      #Get phenotype name for the population gate
      phenotype = ""
      if (!config_table$Phenotype[[i]] == "") {
        phenotype = as.character(config_table$Phenotype[[i]])
      } else {
        phenotype = paste("Gate", i)
      }
      
      #resolve label names of current and parent populations
      currentPopID = paste("pop", i, sep = "")
      parentPopID = paste("pop", config_table$ParentIndex[i], sep = "")
      #display current gating level to standard output
      print(
        paste(
          "Current cell population:",
          paste(currentPopID, " (", phenotype, ")", sep = ""),
          "Parent:",
          paste(parentPopID, " (", config_table$Phenotype[config_table$ParentIndex[i]], ")", sep =
                  ""),
          "Gating on:",
          colnames(fcs)[config_table$ChannelX[i]],
          "/",
          colnames(fcs)[config_table$ChannelY[i]]
        )
      )
      
      
      #Check if the parent is whole sample. Retrieve filtered fcs object if it is a subset.
      if (config_table$ParentIndex[i] == 0) {
        currentFlowFrame = fcs
        popMatrix <- matrix(1, nrow(fcs), 1)
        colnames(popMatrix) <- c(currentPopID)
      } else{
        currentFlowFrame = popFlowFrameList[[config_table$ParentIndex[i]]]
        popMatrix<-`colnames<-`(cbind(popMatrix, matrix(1,nrow(popMatrix),1)), c(colnames(popMatrix), currentPopID))
      }
      
      
      #Check filtering type (0=filter by cluster centroids, 1=simple rectangular bisecting filter, 2=slope bisecting filter)
      if (config_table$FilterType[i] == 0 ||
          config_table$FilterType[i] == 1) {
        #Use rectangular gating boundary for both bisecting and cluster filtering methods
        popGatingMatrixList[[i]] <-
          matrix(c((config_table$minX[i] * 4096 / 200),
                   (config_table$maxX[i] * 4096 / 200),
                   (config_table$minY[i] * 4096 / 200),
                   (config_table$maxY[i] * 4096 / 200)
          ),
          ncol = 2,
          dimnames = list(c("min", "max"),
                          c(colnames(fcs)[config_table$ChannelX[i]], colnames(fcs)[config_table$ChannelY[i]])))
        
        #check and apply reverse filtering if specified
        if (!is.null(reverse_path)) {
          if (i <= nrow(reverse_table)) {
            #Get rectangular reverse gating boundary
            reverse_popGatingMatrixList[[i]] <-
              matrix(c((reverse_table$minX[i] * 4096 / 200),
                       (reverse_table$maxX[i] * 4096 / 200),
                       (reverse_table$minY[i] * 4096 / 200),
                       (reverse_table$maxY[i] * 4096 / 200)
              ),
              ncol = 2,
              dimnames = list(c("min", "max"),
                              c(
                                colnames(fcs)[config_table$ChannelX[i]], colnames(fcs)[config_table$ChannelY[i]]
                              )))
            popGatingFilterList[[i]] <-
              rectangleGate(filterId = phenotype, .gate = popGatingMatrixList[[i]]) &
              !rectangleGate(filterId = paste(phenotype, "_rev", sep = ""),
                             .gate = reverse_popGatingMatrixList[[i]])
            
          } else {
            popGatingFilterList[[i]] <-
              rectangleGate(filterId = phenotype, .gate = popGatingMatrixList[[i]])
          }
        } else {
          popGatingFilterList[[i]] <-
            rectangleGate(filterId = phenotype, .gate = popGatingMatrixList[[i]])
        }
        
        #check if parent is whole sample and initiate clustering algorithm
        if (config_table$ParentIndex[i] == 0)
        {
          print("Creating initial clustering on whole sample")
          clusterResults[[parentPopID]] <-
            cluster_run(
              currentFlowFrame,
              colsToUseList[[1]],
              k = k,
              xdims = xdims,
              ydims = ydims,
              method = method,
              initializer = initializer,
              debugmode = debugmode
            )
        }
        
        #For cluster filtering
        if (config_table$FilterType[i] == 0) {
          if (isTRUE(debugmode)) {
            cat("Filtering by cluster centroids...")
          }
          
          pop_cluster_list[[i]] <-
            filter(clusterResults[[parentPopID]]$flowFrame, popGatingFilterList[[i]])
          
          temp_include_table = pop_cluster_list[[i]]@subSet
          
          if (isTRUE(debugmode)) {
            cat("define logical vector via gating on centroids...")
          }
          pop_filt_list[[i]] <-
            apply(as.array(clusterResults[[parentPopID]]$mapping), 1, function(x)
              temp_include_table[x])
          
          if (isTRUE(debugmode)) {
            cat("events subsetting...\n")
          }
          popFlowFrameList[[i]] <-
            Subset(currentFlowFrame, pop_filt_list[[i]])
          
          #plotting 2D scatter plots for filtered population and cluster centroids
          if (plotCentroids == TRUE) {
            if (nrow(popFlowFrameList[[i]]) > 0) {
              if (isTRUE(debugmode)) {
                cat("Plotting 2D scatter plot...")
              }
              png(
                paste(sampleName, "_", currentPopID, ".png", sep = ""),
                width = 4,
                height = 4,
                units = "in",
                res = 600
              )
              par(bg = 'grey')
              plot(
                fcs,
                main = sampleName,
                sub = config_table$Phenotype[i],
                c(colnames(fcs)[config_table$ChannelX[i]], colnames(fcs)[config_table$ChannelY[i]]),
                smooth =
                  smoothPlot,
                col = "white",
                pch = '.',
                cex = 0.1,
                xlim = c(0, 4096),
                ylim = c(0, 4096)
              )
              if (config_table$ParentIndex[i] > 0) {
                par(new = TRUE)
                plot(
                  popFlowFrameList[[config_table$ParentIndex[i]]],
                  c(colnames(fcs)[config_table$ChannelX[i]], colnames(fcs)[config_table$ChannelY[i]]),
                  smooth =
                    smoothPlot,
                  col = "black",
                  pch = '.',
                  cex = 0.1,
                  xlim = c(0, 4096),
                  ylim = c(0, 4096)
                )
              }
              par(new = TRUE)
              plot(
                popFlowFrameList[[i]],
                c(colnames(fcs)[config_table$ChannelX[i]], colnames(fcs)[config_table$ChannelY[i]]),
                smooth =
                  smoothPlot,
                col = "red",
                pch = '.',
                cex = 0.1,
                xlim = c(0, 4096),
                ylim = c(0, 4096)
              )
              #xaxis=as.character(colnames(fcs)[config_table$ChannelX[i]])
              #yaxis=as.character(colnames(fcs)[config_table$ChannelY[i]])
              #xyplot(data=popFlowFrameList[[i]], channel.x=xaxis, channel.y=yaxis, smooth = FALSE)
              dev.off()
              if (!is.null(pop_cluster_list[[i]])) {
                if (isTRUE(debugmode)) {
                  cat("Plotting cluster centroids overlay...\n")
                }
                png(
                  paste(sampleName, "_", currentPopID, "_c.png", sep = ""),
                  width = 4,
                  height = 4,
                  units = "in",
                  res = 300
                )
                par(bg = NA)
                plot(
                  popFlowFrameList[[i]],
                  main = sampleName,
                  c(colnames(fcs)[config_table$ChannelX[i]], colnames(fcs)[config_table$ChannelY[i]]),
                  smooth =
                    smoothPlot,
                  xlim = c(0, 4096),
                  ylim = c(0, 4096)
                )
                par(new = TRUE)
                plot(
                  Subset(clusterResults[[parentPopID]]$flowFrame, pop_cluster_list[[i]]),
                  c(colnames(fcs)[config_table$ChannelX[i]], colnames(fcs)[config_table$ChannelY[i]]),
                  pch = "+",
                  col = "red",
                  smooth = FALSE,
                  xlim = c(0, 4096),
                  ylim = c(0, 4096)
                )
                dev.off()
              }
            } else if (!is.null(clusterResults[[parentPopID]]$flowFrame)) {
              if (isTRUE(debugmode)) {
                cat(
                  "No cell found for ",
                  currentPopID,
                  ";plotting cluster centroids overlay of parent population...\n"
                )
              }
              png(
                paste(currentPopID, "_noevents_from_parent.png", sep = ""),
                width = 4,
                height = 4,
                units = "in",
                res = 300
              )
              plot(
                popFlowFrameList[[config_table$ParentIndex[i]]],
                c(colnames(fcs)[config_table$ChannelX[i]], colnames(fcs)[config_table$ChannelY[i]]),
                smooth =
                  smoothPlot,
                xlim = c(0, 4096),
                ylim = c(0, 4096)
              )
              par(new = TRUE)
              plot(
                clusterResults[[parentPopID]]$flowFrame,
                c(colnames(fcs)[config_table$ChannelX[i]], colnames(fcs)[config_table$ChannelY[i]]),
                pch = "+",
                col = "red",
                smooth = smoothPlot,
                xlim = c(0, 4096),
                ylim = c(0, 4096)
              )
              dev.off()
            }
          }
          
          print(
            paste(
              "#Parent events: ",
              nrow(currentFlowFrame),
              "     #Filtered events: ",
              nrow(popFlowFrameList[[i]]),
              " (",
              format(
                100 * nrow(popFlowFrameList[[i]]) / nrow(currentFlowFrame),
                digits = 2,
                nsmall = 2
              ),
              "%)",
              sep = ""
            )
          )
          
        } else {
          #using bisecting approach to filter on the cell events directly via the specified gates
          if (isTRUE(debugmode)) {
            cat("Filtering by events (bisecting)...\n")
          }
          pop_filt_list[[i]] <-
            filter(currentFlowFrame, popGatingFilterList[[i]])
          pop_cluster_list[[i]] <-
            filter(clusterResults[[parentPopID]]$flowFrame, popGatingFilterList[[i]])
          temp_include_table = pop_cluster_list[[i]]@subSet
          popFlowFrameList[[i]] <-
            Subset(currentFlowFrame, pop_filt_list[[i]])
          
          print(summary(pop_filt_list[[i]]))
          print(paste(
            "Parent events:",
            nrow(currentFlowFrame),
            "Filtered list events:",
            nrow(popFlowFrameList[[i]])
          ))
          popFlowFrameList[[i]] <-
            Subset(currentFlowFrame, pop_filt_list[[i]])
          
          print(
            paste(
              "filter cluster mapping",
              parentPopID,
              " using logicalVector from pop",
              i
            )
          )
          
          print(paste(
            "cluster map",
            currentPopID,
            "has length",
            nrow(clusterResults[[currentPopID]]$mapping)
          ))
          
        }
        
      } else if (config_table$FilterType[i] == 2) {
        #simulating the slope method from the original DAFi implementation (to be updated with polygon method in the future)
        if (isTRUE(debugmode)) {
          cat("Filtering by events slope definition...\n")
        }
        popGatingMatrixList[[i]] <-
          matrix(c(
            0,
            4096 * (config_table$maxY[i] / config_table$maxX[i]),
            4096,
            4096,
            0,
            4096,
            4096,
            4096 * (config_table$minY[i] / config_table$minX[i])
          ),
          ncol = 2,
          nrow = 4)
        
        colnames(popGatingMatrixList[[i]]) <-
          c(colnames(fcs)[config_table$ChannelX[i]], colnames(fcs)[config_table$ChannelY[i]])
        
        popGatingFilterList[[i]] <-
          polygonGate(filterId = phenotype, .gate = popGatingMatrixList[[i]])
        
        pop_filt_list[[i]] <-
          filter(currentFlowFrame, popGatingFilterList[[i]])
        print(summary(pop_filt_list[[i]]))
        pop_cluster_list[[i]] <-
          filter(clusterResults[[parentPopID]]$flowFrame, popGatingFilterList[[i]])
        
        popFlowFrameList[[i]] <-
          Subset(currentFlowFrame, pop_filt_list[[i]])
        
        
        print(paste(
          "cluster map",
          currentPopID,
          "has length",
          nrow(clusterResults[[currentPopID]]$mapping)
        ))
        
        print(paste(
          "Parent events:",
          nrow(currentFlowFrame),
          "Filtered list events:",
          nrow(popFlowFrameList[[i]])
        ))
      }
      
      cat("processing popMatrix...\n")
      currentFilterExpr <- exprs(popFlowFrameList[[i]])
      currentFilterVector <- currentFilterExpr[,"Event"]
      print(paste("event# length",length(currentFilterVector)))
      for (k in 1:length(currentFilterVector)) {
        popMatrix[currentFilterVector[k],currentPopID]<-0
        finalPopMatrix[currentFilterVector[k]]<-i
        #print(paste("event#",currentFilterVector[k]))
      }
      
      
      #clone clustering results from parent population
      clusterResults[[currentPopID]] <-
        clusterResults[[parentPopID]]
      
      #reset clustering output object
      clusterResults[[currentPopID]]$Object <- NULL
      
      
      if (length(pop_filt_list[[i]]) == 1) {
        logicalVector = pop_filt_list[[i]]@subSet
      } else {
        logicalVector = pop_filt_list[[i]]
      }
      oldMapping <- clusterResults[[parentPopID]]$mapping
      
      clusterResults[[currentPopID]]$mapping <-
        oldMapping[logicalVector]
      
      #exclude the current x and y markers for downstream clustering/filtering
      if ((config_table$ParentIndex[i] == 0) |
          excludeGated == FALSE) {
        colsToUseList[[i]] <- colsToUse
      } else{
        colsToUseList[[i]] <-
          setdiff(colsToUseList[[config_table$ParentIndex[i]]], as.integer(
            which(
              fcs@parameters[[2]] == colnames(fcs)[config_table$ChannelX[i]] |
                fcs@parameters[[2]] == colnames(fcs)[config_table$ChannelY[i]]
            )
          ))
      }
      #Check to see if reclustering is specified for this population
      if ((config_table$Recluster[i] == 1) &
          (nrow(popFlowFrameList[[i]]) > 100)) {
        if (is.null(clusterResults[[currentPopID]]$Object))
        {
          print(paste("Reclustering on", currentPopID))
          clusterResults[[currentPopID]] <-
            cluster_run(
              popFlowFrameList[[i]],
              colsToUseList[[i]],
              k = k,
              xdims = xdims,
              ydims = ydims,
              method = method,
              initializer = initializer,
              debugmode = debugmode
            )
        }
      }
      
      #Write population events and percentages to file
      if (writeOutput == TRUE) {
        if (!(config_table$ParentIndex[i] == 0)) {
          parentEvents = nrow(popFlowFrameList[[config_table$ParentIndex[i]]])
        } else {
          parentEvents = nrow(fcs)
        }
        write(
          paste(i, 100 * nrow(popFlowFrameList[[i]]) / parentEvents, sep = "\t"),
          file = "pop_percentage.txt",
          append = TRUE
        )
        write(paste(i, nrow(popFlowFrameList[[i]]), sep = "\t"),
              file = "pop_events.txt",
              append = TRUE)
      }
    }

    if (writeOutput == TRUE) {
      
      
      cdata <- exprs(fcs)
      cdata <- cbind(cdata[,1:target,drop=F], popMatrix, cdata[,(target+1):ncol(cdata),drop=F], finalPopMatrix)
      
      dafi_filename=paste(strsplit(tail(unlist(strsplit(rawfcs_path, .Platform$file.sep)),n=1),"[.]")[[1]][1],"dafi","txt",sep=".")
      
      write.table(cdata, file=dafi_filename, quote=F, row.names=F,
                  col.names=T, sep='\t', append=F)
     
      
    }
  }
