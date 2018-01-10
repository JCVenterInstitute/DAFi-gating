#Accessory codes for FCSTrans Transformation of FCS

#' Set output to 0 when input is less than cutoff value
#' 
#' @param x
#' @param cutoff
#' @param target
#' 
ipfloor <- function (x, cutoff = 0, target = 0) {
  y = x
  if (x <= cutoff)
    y = target
  y
}

#' Set output to 0 when input is less than cutoff value
#'
#' @param x
#' @param cutoff
#' @param target
#' 
ipceil <- function (x, cutoff = 0, target = 0) {
  y = x
  if (x >= cutoff)
    y = target
  y
}

#' immport linear function - convert scatter values to channel output
#' linear transformation
#' @param x
#' @param range
#' @param cutoff
#' @param channelrange
#' 
ipscatter <-
  function (x,
            range = 4096.0,
            cutoff = 4096,
            channelrange = 262144) {
    y = range * x / channelrange
    y = sapply(y, ipfloor)
    y = sapply(y, ipceil, cutoff = cutoff, target = range)
    y
  }

# immport time function - convert time values to channel output
# linear transformation
iptime <- function (x, channelrange) {
  # use simple cutoff for now
  y = sapply(x, ipfloor)
  y
}

#linear transformation for the scatter channels
scatterTransform <-
  function(transformationId = "defaultScatterTransform",
           channelrange = 262144,
           range = 4096,
           cutoff = 4096,
           rescale = TRUE) {
    t <- new(
      "transform",
      .Data = function(x) {
        x <- ipscatter(x, channelrange = channelrange)
      }
    )
    t
  }

#linear transformation for the time channel
timeTransform <- function(transformationId = "defaultTimeTransform",
                          channelrange = 262144,
                          range = 4096,
                          rescale = TRUE) {
  t <- new(
    "transform",
    .Data = function(x) {
      x <- iptime(x, channelrange = channelrange)
    }
  )
  t
}

FCSTransInput <- function(rawfcs_path) {
  
  #Construct output list
  output <- list(
    rawflowFrame = NULL,
    processedFlowFrame = NULL,
    colsToUse = NULL
  )
  
  #Read in input raw fcs file
  fcs_raw <-
    suppressWarnings(flowCore::read.FCS(rawfcs_path))
  
  output$rawflowFrame <- fcs_raw
  
  #Check and apply compensation
  spill = keyword(fcs_raw)$SPILL
  
  if (is.null(spill) == FALSE) {
    tryCatch({
      fcs_raw = compensate(fcs_raw, spill)
    }, error = function(ex) {
      str(ex)
      
    })
  }
  
  #Locate columns (marker channels) used for clustering
  colsToUseList <- list()
  colsToUse = as.integer(which(!is.na(fcs_raw@parameters[[2]]), TRUE))
  colsToUseList[[1]] <- colsToUse
  #rename channels to marker type
  colnames(fcs_raw)[colsToUse] <-
    as.character(fcs_raw@parameters[[2]])[colsToUse]
  
  output$colsToUse <- colsToUse
  
  # transformation: linear on time and scatter channels, FCSTrans on marker channels
  listC = colnames(fcs_raw)[colsToUse]
  listS = colnames(fcs_raw)[1:colsToUse[1] - 1]
  #timeC = colnames(fcs_raw)[17]
  
  scattT <- transformList(listS, scatterTransform())
  #timeT <- transformList(timeC, timeTransform())
  lgcl <- transformList(listC, FCSTransTransform())
  
  fcs <- transform(fcs_raw, scattT)
  #fcs <- transform(fcs, timeT)
  fcs <- transform(fcs, lgcl)
  output$processedFlowFrame <- fcs
  
  output
}