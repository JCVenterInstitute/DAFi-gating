#Accessory codes for FCSTrans Transformation of FCS

#' Set output to 0 when input is less than cutoff value
#' 
#' @param x input array
#' @param cutoff cutoff number
#' @param target target number
#' 
ipfloor <- function (x, cutoff = 0, target = 0) {
  y = x
  if (x <= cutoff)
    y = target
  y
}

#' Set output to 0 when input is less than cutoff value
#'
#' @param x input array
#' @param cutoff cutoff number
#' @param target target number
#' 
ipceil <- function (x, cutoff = 0, target = 0) {
  y = x
  if (x >= cutoff)
    y = target
  y
}

#' immport linear function - convert scatter values to channel output
#' linear transformation
#' @param x input array
#' @param range range number
#' @param cutoff cutoff number
#' @param channelrange original channel range
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

#' immport time function - convert time values to channel output
#' linear transformation
#' 
#' @param x input array
#' @param channelrange range of channel
iptime <- function (x, channelrange) {
  # use simple cutoff for now
  y = sapply(x, ipfloor)
  y
}

#' linear transformation for the scatter channels
#' 
#' @param transformationId name of transformation
#' @param channelrange original channel range
#' @param range total range
#' @param cutoff cutoff value
#' @param rescale whether to rescale or not
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

#' linear transformation for the time channel
#' 
#' @param transformationId name of transformation
#' @param channelrange original channel range
#' @param range total range
#' @param rescale whether to rescale or not
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

#' Read in raw FCS file and transform it using FCSTrans
#'
#' also automatically relabel the columns using marker labels and find colsToUse
#' 
#' @param rawfcs_path     path to the raw fcs file
#' @export
FCSTransInput <- function(rawfcs_path, debugmode=FALSE) {
  
  
  #Construct output list
  output <- list(
    rawflowFrame = NULL,
    processedFlowFrame = NULL,
    colsToUse = NULL
  )
  
  if(isTRUE(debugmode)){print("Reading FCS file")}
  #Read in input raw fcs file
  fcs_raw <-
    suppressWarnings(read.FCS(rawfcs_path))
  
  if(isTRUE(debugmode)){print("Running compensation...")}
  #Check and apply compensation
  spill = keyword(fcs_raw)$SPILL
  
  if (is.null(spill) == FALSE) {
    tryCatch({
      fcs_raw = compensate(fcs_raw, spill)
    }, error = function(ex) {
      str(ex)
      
    })
  }
  
  if(isTRUE(debugmode)){print("Locating columns used for clustering...")}
  #Locate columns (marker channels) used for clustering
  colsToUse = as.integer(which(!is.na(fcs_raw@parameters[[2]]), TRUE))
  #rename channels to marker type
  if(isTRUE(debugmode)){
  print("renaming channels to marker type...")
  print(colnames(fcs_raw)[colsToUse])
  print(as.character(fcs_raw@parameters[[2]])[colsToUse])}
  colnames(fcs_raw)[colsToUse] <-
    as.character(fcs_raw@parameters[[2]])[colsToUse]
  
  if(isTRUE(debugmode)){print("Transformation via FCSTrans...")}
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
  
  output$rawflowFrame <- fcs_raw
  output$colsToUse <- colsToUse
  output$processedFlowFrame <- fcs
  
  output
}