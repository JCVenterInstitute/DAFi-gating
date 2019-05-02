#
# FCS conversion program
#
#
suppressMessages(library(parallel))
suppressMessages(library(flowCore))
suppressMessages(library(openxlsx))
suppressMessages(library(tools))

anonymous <- function(name) {
  if (is.null(nametable[[name]])) {
    assign("gbcount", gbcount + 1, envir = .GlobalEnv)
    assign(name, paste("Sub", formatC(
      gbcount, width = 3, flag = "0"
    ), sep = ""), envir = nametable)
  }
  nametable[[name]]
}

# set output to 0 when input is less than cutoff value
ipfloor <- function (x, cutoff = 0, target = 0) {
  y = x
  if (x <= cutoff)
    y = target
  y
}

# set output to 0 when input is less than cutoff value
ipceil <- function (x, cutoff = 0, target = 0) {
  y = x
  if (x >= cutoff)
    y = target
  y
}

# Convert fluorescent values to channel output using log transformation
# log10 transformation, not needed because of the use of FCSTransTransform
#iplog <- function(x) {
#  x = sapply(x, ipfloor, cutoff = 1, target = 1)
#  y = 1024 * log10(x) - 488.6
#  y
#}

# immport linear function - convert scatter values to channel output
# linear transformation
ipscatter <- function (x, channelrange = 262144) {
  y = 4095.0 * x / channelrange
  y = sapply(y, ipfloor)
  y = sapply(y, ipceil, cutoff = 4095, target = 4095)
  y
}

# immport time function - convert time values to channel output
# linear transformation for time parameter, not needed because of the use of FCSTransTransform
#iptime <- function (x, channelrange) {
  # use simple cutoff for now
#  y = sapply(x, ipfloor)
#  y
#}


# get marker type
getMarkerType <- function(name) {
  type = ""
  prefix2 = toupper(substr(name, 1, 2))
  prefix3 = toupper(substr(name, 1, 3))
  prefix4 = toupper(substr(name, 1, 4))
  if (prefix2 == "FS" || prefix2 == "SS") {
    type = "SCATTER"
  } else if (prefix3 == "FSC" || prefix3 == "SSC") {
    type = "SCATTER"
  } else if (prefix4 == "TIME") {
    type = "TIME"
  } else {
    pieces = unlist(strsplit(name, "-"))
    if (toupper(pieces[length(pieces)]) == "A") {
      type = "FLUO_AREA"
    } else {
      type = "FLUO_NON_AREA"
    }
  }
  #print(paste("marker:", name, ", type:", type))
  type
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

# immport convert function - convert flow cytometry values to channel output
# iterate columns name and treat data in three categories:
#   scatter      linear
#   time         linear
#   fluorescent  logicle
# example:
#   > ad008 = read.FCS('ad008.fcs', transformation=F)
#   > ad008_ipc = convertfcs(ad008r)
# @param fcs    flowFrame object
convertfcs <- function(fcs_raw, compen = "internal") {
  debug = TRUE
  #check file type and FCS version
  if (class(fcs_raw)[1] != "flowFrame") {
    print("convertfcs requires flowFrame object as input")
    return(FALSE)
  }
  keywords = keyword(fcs_raw)
  
  if (debug)
    print(paste("FCS version:", keywords$FCSversion))
  if (debug)
    print(paste("File data type:", keywords['$DATATYPE']))
  if (keywords$FCSversion == "2" ||
      keywords$FCSversion == "3" || keywords$FCSversion == "3.1") {
    datatype = unlist(keywords['$DATATYPE']) #check data type
    if (datatype == 'F') {
      #apply compensation if available
      if (compen == "internal") {
        print(paste("Applying compensation from fcs ", compen))
        spill = keyword(fcs_raw)$SPILL
      } else if (file_ext(compen) == "xlsx") {
        print(paste("Applying compensation matrix file: ", compen))
        spill = as.matrix(read.xlsx(compen,
                                    rowNames = TRUE, colNames = TRUE))
      } else if (file_ext(compen) == "csv") {
        print(paste("Applying compensation matrix file: ", compen))
        spill = as.matrix(read.csv(
          compen,
          header = TRUE,
          row.names = 1,
          check.names = FALSE
        ))
      } else {
        spill = NULL
      }
      
      if (debug)
        print("check spill")
      if (is.null(spill) == FALSE) {
        colnames(spill) <- gsub(
          x = colnames(spill),
          pattern = "\\.",
          replacement = " "
        )
        tryCatch({
          fcs_raw = compensate(fcs_raw, spill)
        }, error = function(ex) {
          str(ex)
          
        })
        print ("applied compensation!")
      }
      
      
      markers = colnames(fcs_raw)
      markerCols = list()
      if (debug)
        print("loop through markers")
      scatterlist = list()
      markerlist = list()
      
      if (debug) {
        print("Transformation via FCSTrans...")
      }
      
      fcs = fcs_raw
      
      for (i in 1:length(markers)) {
        rangekeyword = paste("$P", i, "R", sep = "")
        #if (debug) print(paste("  range keyword:", rangekeyword))
        print(markers[i])
        channelrange = as.numeric(keywords[rangekeyword])
        if (debug)
          print(paste(" range value:", as.character(channelrange)))
        markertype = getMarkerType(markers[i])
        print(markertype)
        if ((markertype != "SCATTER") & (markertype != "TIME")) {
          print("transforming marker channel")
          markerCols <- c(markerCols, as.integer(i))
          markerlist <-
            transformList(markers[i],
                          FCSTransTransform(channelrange = channelrange))
          fcs <- transform(fcs, markerlist)
          print("done transforming marker channel")
        } else if (markertype == "SCATTER") {
          print("transforming scatter channel")
          scatterlist <-
            transformList(markers[i],
                          scatterTransform(channelrange = channelrange))
          fcs <- transform(fcs, scatterlist)
          print("done transforming scatter channel)")
        } else if (markertype == "TIME") {
          fcs <-
            transform(fcs, transformList(
              markers[i],
              timeTransform(channelrange = channelrange)
            ))
        }
      }
      
      colsToUse = as.integer(which(!is.na(fcs_raw@parameters[[2]]), TRUE))
      
      if (length(colsToUse) == length(markerCols)) {
        #rename channels to marker type
        if (debug) {
          print("renaming channels to marker type...")
          print(colnames(fcs_raw)[colsToUse])
          print(as.character(fcs_raw@parameters[[2]])[colsToUse])
        }
        colnames(fcs_raw)[colsToUse] <-
          as.character(fcs_raw@parameters[[2]])[colsToUse]
      } else {
        colsToUse = as.integer(markerCols)
      }
      
      
      # transformation: linear on time and scatter channels, FCSTrans on marker channels
      #listC = colnames(fcs_raw)[colsToUse]
      #listS = colnames(fcs_raw)[1:colsToUse[1] - 1]
      #timeC = colnames(fcs_raw)[17]
      
      #scattT <- transformList(listS, scatterlist)
      #timeT <- transformList(timeC, timeTransform())
      #lgcl <- transformList(listC, markerlist)
      
      #fcs <- transform(fcs_raw, scattT)
      #print("Done with transformation of scatter channels")
      #fcs <- transform(fcs, timeT)
      #fcs <- transform(fcs, lgcl)
      #print("Done with transformation of non-scatter channels")
    } else if (datatype == 'I') {
      fcs = fcs_raw
    } else {
      print(paste("Data type", datatype, "in FCS 3 is not supported"))
      fcs = NULL
    }
  } else {
    print(paste("FCS version", keyword(fcs)$FCSversion, "is not supported"))
    fcs = NULL
  }
  fcs
}

# convert fcs into channel output and save result
# @param fcsfile
# @param verbose       Prints file and folder names when TRUE
# @param overwrite     Convert and overwrite existing channel files
# @param renameSource  Rename fcs file to have .fcs suffix
# @param convert       Flag for convert FCS file
# @param savekeys      Flag for save keywords to file
#
# example:
#   > convertfcs('~/data/ad008')
#   This will convert ad008.fcs and create ad008_ip_channel.txt
processfcsfile <-
  function(fcsfile,
           compen = "internal",
           verbose = F,
           overwrite = F,
           renameSource = F,
           removeName = T,
           copyRaw = T,
           convert = T,
           savekeys = T) {
    print("in processfcsfile")
    print(
      paste(
        "parameters:",
        "verbose",
        verbose,
        "overwrite",
        overwrite,
        "renameSource",
        renameSource,
        "removeName",
        removeName,
        "copyRaw",
        copyRaw,
        "convert",
        convert,
        "savekeys",
        savekeys,
        sep = " "
      )
    )
    
    isvalid = F
    tryCatch({
      isvalid = isFCSfile(fcsfile)
    }, error = function(ex) {
      print (paste("    ! Error in isFCSfile", ex))
    })
    
    if (isvalid) {
      pieces = unlist(strsplit(fcsfile, .Platform$file.sep))
      filename = pieces[length(pieces)]
      
      tryCatch({
        fcs_raw <- read.FCS(fcsfile,
                            transformation = F,
                            emptyValue = FALSE)
        
        
        tempFIL <- fcs_raw@description$`$FIL`
        tempSRC <- fcs_raw@description$`$SRC`
        tubeName <- fcs_raw@description$`TUBE NAME`
        srcformat = F
        realname = ""
        
        if (!is.null(tempSRC)) {
          srcpieces <- unlist(strsplit(tempSRC, " "))
          if (length(srcpieces) > 1) {
            if (grepl("^[A-Z][a-z]+$", srcpieces[1])) {
              realname = srcpieces[1]
              srcformat = T
              if (grepl("^[A-Z][a-z]+$", srcpieces[2])) {
                realname = paste(realname, srcpieces[2], sep = " ")
              }
            } else if (grepl("^[A-Z][a-z]*[A-Z][-][0-9][0-9]", srcpieces[1])) {
              realname = sub("[A-Z][-][0-9][0-9]", "", srcpieces[1])
              secondTerm = sub("^[A-Z][a-z]*", "", srcpieces[1])
              srcpieces[1] = secondTerm
              srcpieces <- c(c(realname), srcpieces)
              srcformat = T
            }
          }
        } else {
          srcpieces = NULL
        }
        sampleNumber = ""
        
        tryCatch({
          if (grepl("_", tempFIL)) {
            filpieces <- unlist(strsplit(tempFIL, "_"))
            sampleNumber <-
              (unlist(strsplit(filpieces[length(filpieces)], ".")))[1]
          } else {
            filpieces <- unlist(strsplit(tempFIL, " "))
            sampleNumber <-
              (unlist(strsplit(filpieces[length(filpieces)], ".")))[1]
          }
          
          sampleNumber <- filpieces[length(filpieces)]
          
        }, error = function(ex) {
          if (verbose)
            print (ex)
          print ("ignore and continue")
        })
        
        
        if (removeName) {
          if (verbose)
            print(paste("    Removing Names", fcsfile))
          
          if (srcformat) {
            outputfilename <-
              paste(tubeName,
                    paste(srcpieces[2:length(srcpieces)], collapse = '_'),
                    sampleNumber,
                    sep = "_")
            fcs_raw@description$`$SRC` <- outputfilename
            fcs_raw@description$`$FIL` <- outputfilename
            fcs_raw@description$FILENAME <- outputfilename
          } else {
            outputfilename <- fcs_raw@description$ORIGINALGUID
            fcs_raw@description$`$SRC` <- NULL
            fcs_raw@description$`$FIL` <- NULL
            fcs_raw@description$FILENAME <- "redacted"
          }
        } else {
          #filepieces = unlist(strsplit(filename, '\\.'))
          
          #outputfilename <- filepieces[length(filepieces)-1]
          
          outputfilename <- filename
        }
        
        pdata <- keyword(fcs_raw)  #extract keywords from fcs file
        if (copyRaw) {
          write.FCS(fcs_raw, paste(outputfilename, "raw.fcs", sep = "_"))
        }
        if (convert) {
          if (verbose)
            print(paste("    Converting", fcsfile))
          
          fcs <- convertfcs(fcs_raw, compen = compen)
          
          if (!is.null(fcs)) {
            cdata <- exprs(fcs)
            
            #write.FCS(fcs, paste(outputfilename,"fcs",sep="."))
            #if (verbose) print(paste("      Converted to", outputfilename))
            channelfile <- paste(outputfilename, "txt", sep = ".")
            write.table(
              cdata,
              file = channelfile,
              quote = F,
              row.names = F,
              col.names = T,
              sep = '\t',
              append = F
            )
            if (verbose)
              print(paste("      Converted to", channelfile))
          }
        }
        if (savekeys) {
          keylistfile <- paste(outputfilename, "lst", sep = ".")
          write.table(
            as.matrix(pdata),
            file = keylistfile,
            quote = F,
            row.names = T,
            col.names = F,
            sep = '=',
            append = F
          )
          if (verbose)
            print(paste("      Keywords saved to", keylistfile))
        }
        
        tempName = NULL
        if (srcformat) {
          tempName <- realname
        }
        output = c(filename, outputfilename, tempName)
      }, error = function(ex) {
        if (verbose)
          print (paste("    ! Error in reading in FCS ", ex))
        output = NULL
      })
    } else {
      output = NULL
    }
    
    output
  }

#' Main recursive function to borrow through directory trees to process fcs files
#'
#'
#' @param entry         Can be fcs object, fcs file or folder name
#' @param verbose       Prints file and folder names when TRUE
#' @param overwrite     Convert and overwrite existing channel files
#' @param renameSource  Rename fcs file to have .fcs suffix
#' @param removeName    Remove patient name from $Fil and $Src fields
#'
# example:
#   > ipconvert('~/data/test')
#   > ipconvert(c('/tmp/a.fcs', 't101.fcs', '/data')
ipconvert <-
  function(entry,
           compen = "internal",
           verbose = F,
           overwrite = F,
           renameSource = F,
           removeName = F,
           copyRaw = F,
           convert = T,
           savekeys = T,
           level = 1) {
    success = F
    print(paste("Running ipconvert in level", level, " and entry ", entry, sep =
                  " "))
    if (level == 1) {
      assign("filetable", new.env(has = TRUE), envir = .GlobalEnv)
      assign("nametable", new.env(has = TRUE), envir = .GlobalEnv)
      assign("gbcount", 0, envir = .GlobalEnv)
      originalwd <- getwd()
      # dir.create(file.path(".", "export_fcs"), showWarnings = FALSE)
      # setwd(file.path(".", "export_fcs"))
      assign("no_cores", detectCores() - 1, envir = .GlobalEnv)
      assign("cl", makeCluster(no_cores, outfile = "cltest.out"), envir = .GlobalEnv)
      clusterExport(cl, ls(.GlobalEnv), envir = .GlobalEnv)
      clusterEvalQ(cl, library(flowCore))
      clusterEvalQ(cl, library(tools))
      clusterEvalQ(cl, library(openxlsx))
      
      #write(paste("Original","New", sep = "\t"), file="filetable.out",append=TRUE)
      #write(paste("Original","New", sep = "\t"), file="nametable.out",append=TRUE)
      clusterExport(cl, ls(environment()), envir = environment())
      
    }
    
    if (compen == "internal") {
      if (file.exists(compen)) {
      print(paste("Applying compensation from fcs ", compen))
      }else {
        print("comp file not found!")
      }
      
    } else if (file_ext(compen) == "xlsx") {
      if (file.exists(compen)) {
        print(paste("Applying compensation matrix file: ", compen))
      }else {
        print("comp file not found!")
      }
      
    } else if (file_ext(compen) == "csv") {
      if (file.exists(compen)) {
        print(paste("Applying compensation matrix file: ", compen))
      }else {
        print("comp file not found!")
      }
    }
    
    entryclass = class(entry)
    
    
    
    if (is.vector(entry) && length(entry) > 1) {
      if (verbose)
        print (paste(level, ":", "vector - size", length(entry)))
      sapply(
        entry,
        ipconvert,
        verbose = verbose,
        overwrite = overwrite,
        renameSource = renameSource,
        removeName = removeName,
        copyRaw = copyRaw,
        convert = convert,
        savekeys = savekeys,
        level = level + 1
      )
    } else if (entryclass == 'flowFrame') {
      if (verbose)
        print (paste(level, ":", "flowFrame"))
      #convertfcs(entry)
    } else if (entryclass == 'character') {
      if (file.exists(entry)) {
        if (file.info(entry)$isdir) {
          if (verbose)
            print (paste(level, ":", "Folder :", entry))
          #get folder name without ending separator
          pattern = paste(.Platform$file.sep, '$', sep = '')
          entry = sub(pattern, '', entry)
          files = list.files(entry)
          if (length(files) > 0) {
            files = paste(entry, files, sep = .Platform$file.sep)
            if (file.info(files[1])$isdir) {
              sapply(
                files,
                ipconvert,
                level = level + 1,
                verbose = verbose,
                overwrite = overwrite,
                renameSource = renameSource,
                removeName = removeName,
                copyRaw = copyRaw,
                convert = convert,
                savekeys = savekeys
              )
            } else {
              print("Processing fcs files")
              print(file.info(files[1]))
              if (length(files) > 1) {
                fileinfolist <-
                  parLapply(cl, files, function(x)
                    processfcsfile(
                      x,
                      compen = compen,
                      verbose = verbose,
                      overwrite = overwrite,
                      renameSource =
                        renameSource,
                      removeName = removeName,
                      copyRaw = copyRaw,
                      convert = convert,
                      savekeys = savekeys
                    ))
              } else {
                processfcsfile(
                  files[1],
                  compen = compen,
                  verbose = verbose,
                  overwrite = overwrite,
                  renameSource = renameSource,
                  removeName = removeName,
                  copyRaw = copyRaw,
                  convert = convert,
                  savekeys = savekeys
                )
              }
              
              # lapply(fileinfolist, function(x){
              #   if(!is.null(x)){
              #     if(!is.null(x[3])){
              #       #newoutputname=paste(anonymous(x[3]),x[2],sep="_")
              #       newoutputname=x[2]
              #       #if(savekeys) file.rename(paste(x[2],"lst",sep="."), paste(newoutputname,"lst",sep="."))
              #       #if(convert) file.rename(paste(x[2],"fcs",sep="."), paste(newoutputname,"fcs",sep="."))
              #       #file.rename(paste(x[2],"_raw.fcs",sep=""), paste(newoutputname,"_raw.fcs",sep=""))
              #
              #       write(paste(x[1],newoutputname, sep = "\t"), file="filetable.txt",append=TRUE)
              #     }else{
              #       write(paste(x[1],x[2], sep = "\t"), file="filetable.txt",append=TRUE)
              #     }
              #   }
              # })
              print("done processing files!")
            }
          }
        } else {
          print("processing fcs")
          
          processfcsfile(
            entry,
            compen = compen,
            verbose = verbose,
            overwrite = overwrite,
            renameSource = renameSource,
            removeName = removeName,
            copyRaw = copyRaw,
            convert = convert,
            savekeys = savekeys
          )
          
        }
      } else {
        if (verbose)
          print(paste(level, ":", entry, "does not exist"))
      }
    } else {
      if (verbose)
        print(paste(level, ":", entry, "is not supported"))
    }
    
    
    if (level == 1) {
      setwd(originalwd)
      # setwd(file.path(".", "export_fcs"))
      
      # namekeylist<-ls(nametable)
      # if (length(namekeylist)>0){
      #
      #   lapply(namekeylist, function(x){write(paste(x,nametable[[x]], sep = "\t"), file="nametable.txt",append=TRUE)})
      # }
      stopCluster(cl)
      closeAllConnections()
      # setwd(originalwd)
      success = T
    }
    success
  }
