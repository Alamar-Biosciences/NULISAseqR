#' NULISAseq Inter-Plate Normalization
#'
#' Does inter-plate normalization using one or more of the following methods:
#' platewise total count (TC; sets total plate counts equal to the number of 
#' samples on the plate * 10^6), inter-plate control (IPC; divides each target 
#' by the median IPC count for that target), intensity normalization (IN; 
#' sets median of each target on each plate to equal the global median). Default
#' is to do IPC only. If multiple options are chosen, they are implemented 
#' in the order above (TC, IPC, and finally IN).
#'   
#' Input is a list of data matrices, one for each plate. 
#' These would typically be normData matrix from the intraPlateNorm function 
#' or the Data matrix of readNULISAseq.R (if no intra-plate normalization is 
#' needed).
#'
#' @param data_list A list of data matrices, one for each plate, 
#' to be normalized together. 
#' @param TC Logical TRUE/FALSE indicating whether a platewise total count 
#' normalization should be done. Default is FALSE.
#' @param IPC Logical TRUE/FALSE indicating whether a IPC 
#' normalization should be done. Default is TRUE.
#' @param IN Logical TRUE/FALSE indicating whether intensity normalization
#' should be done. Default is FALSE.
#' @param IPC_wells List of vectors of either the column indices (numeric) 
#' or the column names (character string) of the inter-plate control wells 
#' for each plate. IPC wells are omitted from the 
#' intensity normalization median calculation. 
#' Required if argument IPC=TRUE and/or IN=TRUE.
#' @param NC_wells List of vectors of either the row indices (numeric) or the 
#' row names (character string) 
#' of the negative control wells for each plate. NC wells are omitted from 
#' normalization methods. Required for all methods. 
#' @param IN_samples optional argument. A list of row
#' indices specifying which subset of samples to use for 
#' intensity normalization step. Default uses all samples
#' except for NCs and IPCs.
#' @param TC_omit optional...
#' @param scaleFactor. Optional numeric value used to rescale all data 
#' after normalizing. Default is 1. This may be desirable if there are 
#' many targets with read counts less than the IC counts, to avoid
#' normalized quantities between 0 and 1 (which will be negative
#' in the log scale).
#'
#' @return 
#' @param interNormData A matrix of normalized count data (not log-transformed).
#' @param normFactor A vector of the normalization factors used for 
#' each sample (column) in the data matrix.
#'
#' @examples
#' 
#'
#' @export
#' 
interPlateNorm - function(countMatrices,
                          IPCs,
                          IPCmethod='targetwise',
                          intensityNorm=FALSE,
                          intensityNormMethod='targetwise',
                          intensityNormSamples=NULL,
                          allICs=NULL){
  counts <- lapply(countMatrices, function(x){
    as.matrix(x[,3:ncol(x)])
  })
  if (TC==TRUE){
    if (is.null(totalCountCols)){
      counts <- lapply(countMatrices, function(x){
        as.matrix(x[,3:ncol(x)])
      })
    } else if (!is.null(totalCountCols)){
      counts <- lapply(countMatrices, function(x){
        as.matrix(x[,totalCountCols])
      })
    }
    counts_combined <- as.matrix(do.call(rbind, counts))
    plate_counts <- unlist(lapply(counts, sum))
    total_counts <- sum(counts_combined)
    total_plate_count <- total_counts / length(plate_counts)
    normFactors <- total_plate_count / plate_counts
    normCountsList <- list()
    for (i in 1:length(counts)){
      normCountsList[[i]] <- countMatrices[[i]][,3:ncol(countMatrices[[i]])]*normFactors[i]
    }
    normCounts <- as.matrix(do.call(rbind, normCountsList))
    # rescale data to equal original total counts
    # but exclude all ICs
    allCounts_combined <- lapply(countMatrices, function(x){
      as.matrix(x[,3:ncol(x)])
    })
    allCounts_combined <- as.matrix(do.call(rbind, allCounts_combined))
    normCounts <- normCounts*(sum(allCounts_combined[,-(allICs-2)], na.rm=TRUE)/sum(normCounts[,-(allICs-2)], na.rm=TRUE))
    normData <- do.call(rbind, countMatrices)
    normData[,3:ncol(normData)] <- normCounts
    normData[,3:ncol(normData)] <- apply(normData[,3:ncol(normData)], 2, as.numeric)
    plateNs <- unlist(lapply(countMatrices, nrow))
    plate <- NULL
    for (i in 1:length(plateNs)){
      plate <- c(plate, rep(i, plateNs[i]))
    }
    normData$plate <- plate
    normData <- normData[,c(1,2,ncol(normData),3:(ncol(normData)-1))]
  }
  
  
  IPC_m <- list()
  counts_IPC <- list()
  if (IPCmethod=='targetwise'){
    for (i in 1:length(counts)){
      IPC_rows <- counts[[i]][IPCs[[i]],]
      IPC_m_i <- apply(IPC_rows, 2, median, na.rm=TRUE)
      # replace 0s with 1s
      IPC_m_i[IPC_m_i==0] <- 1
      IPC_m[[i]] <- IPC_m_i
      counts_IPC[[i]] <- counts[[i]] / rep(IPC_m_i, each = nrow(counts[[i]])) 
      counts_IPC[[i]][!is.finite(counts_IPC[[i]])] <- NA
      # line below makes count total equal to original total
      # counts_IPC_rescaled[[i]] <- counts_IPC[[i]]*(sum(counts[[i]], na.rm=TRUE)/sum(counts_IPC[[i]], na.rm=TRUE))
    }
  } else if (IPCmethod=='targetwise.GM'){
    IPC_medians <- list()
    for (i in 1:length(counts)){
      IPC_rows <- counts[[i]][IPCs[[i]],]
      IPC_m_i <- apply(IPC_rows, 2, median, na.rm=TRUE)
      IPC_medians[[i]] <- IPC_m_i
    }
    IPC_medians_matrix <- as.data.frame(do.call(rbind, IPC_medians))
    IPC_global_medians <- apply(IPC_medians_matrix, 2, median, na.rm=TRUE)
    # replace 0s with 1s
    IPC_global_medians[IPC_global_medians==0] <- 1
    for (i in 1:length(counts)){
      IPC_rows <- counts[[i]][IPCs[[i]],]
      IPC_m_i <- apply(IPC_rows, 2, median, na.rm=TRUE)
      # replace 0s with 1s
      IPC_m_i[IPC_m_i==0] <- 1
      IPC_m[[i]] <- IPC_m_i
      IPC_GM_m_i <- IPC_global_medians / IPC_m_i
      counts_IPC[[i]] <- counts[[i]] * rep(IPC_GM_m_i, each = nrow(counts[[i]])) 
      counts_IPC[[i]][!is.finite(counts_IPC[[i]])] <- NA
      # line below makes count total equal to original total
      # counts_IPC_rescaled[[i]] <- counts_IPC[[i]]*(sum(counts[[i]], na.rm=TRUE)/sum(counts_IPC[[i]], na.rm=TRUE))
    }
  } else if (IPCmethod=='acrossTarget'){
    for (i in 1:length(counts)){
      IPC_rows <- counts[[i]][IPCs[[i]],]
      IPC_m_i <- apply(IPC_rows, 2, median, na.rm=TRUE)
      # replace 0s with 1s
      IPC_m_i[IPC_m_i==0] <- 1
      # take across-target median 
      IPC_m[[i]] <- median(IPC_m_i[-(allICs-2)], na.rm=TRUE)
      counts_IPC[[i]] <- counts[[i]] / IPC_m[[i]]
      counts_IPC[[i]][!is.finite(counts_IPC[[i]])] <- NA
      # line below makes count total equal to original total
      # counts_IPC_rescaled[[i]] <- counts_IPC[[i]]*(sum(counts[[i]], na.rm=TRUE)/sum(counts_IPC[[i]], na.rm=TRUE))
    }
  } else if (IPCmethod=='acrossTarget.GM'){
    all_counts <- do.call(rbind, counts)
    global_medians <- apply(all_counts, 2, median, na.rm=TRUE)
    # replace zeros with 1s
    global_medians[global_medians==0] <- 1
    acrossTarget_global_median <- median(global_medians[-(allICs - 2)], na.rm=TRUE)
    for (i in 1:length(counts)){
      IPC_rows <- counts[[i]][IPCs[[i]],]
      IPC_m_i <- apply(IPC_rows, 2, median, na.rm=TRUE)
      # replace 0s with 1s
      IPC_m_i[IPC_m_i==0] <- 1
      # take across-target median 
      IPC_m[[i]] <- median(IPC_m_i[-(allICs - 2)], na.rm=TRUE)
      IPC_GM_m_i <- acrossTarget_global_median / IPC_m[[i]]
      counts_IPC[[i]] <- counts[[i]] * IPC_GM_m_i
      counts_IPC[[i]][!is.finite(counts_IPC[[i]])] <- NA
      # line below makes count total equal to original total
      # counts_IPC_rescaled[[i]] <- counts_IPC[[i]]*(sum(counts[[i]], na.rm=TRUE)/sum(counts_IPC[[i]], na.rm=TRUE))
    }
  }
  normCountsList <- counts_IPC
  if (intensityNorm==TRUE){
    if (is.null(intensityNormSamples)){
      intensityNormSamples <- list()
      for (i in 1:length(counts)){
        # use all sample rows except IPCs
        intensityNormSamples[[i]] <- c(1:nrow(counts[[i]]))[!1:nrow(counts[[i]]) %in% IPCs[[i]]]
      }
    }
    counts_intensityNorm <- list()
    for (i in 1:length(counts)){
      counts_intensityNorm[[i]] <- counts_IPC[[i]][intensityNormSamples[[i]],]
    }
    if (intensityNormMethod=='targetwise'){
      plate_medians <- lapply(counts_intensityNorm, function(x){
        apply(x, 2, median, na.rm=TRUE)
      })
      # change any 0s to 1s
      plate_medians <- lapply(plate_medians, function(x){
        x[x==0] <- 1
        return(x)
      })
      all_counts <- as.matrix(do.call(rbind, counts_intensityNorm))
      global_medians <- apply(all_counts, 2, median, na.rm=TRUE)
      # change any 0s to 1s
      global_medians[global_medians==0] <- 1
      counts_IPC_global_rescale <- list()
      for (i in 1:length(counts)){
        global_scaling <- global_medians / plate_medians[[i]]
        counts_IPC_global_rescale[[i]] <- counts_IPC[[i]] * rep(global_scaling, each = nrow(counts_IPC[[i]]))
      }
      normCountsList <- counts_IPC_global_rescale
    } else if (intensityNormMethod=='global'){
      plate_medians <- lapply(counts_intensityNorm, function(x){
        median(x[,-(allICs-2)], na.rm=TRUE)
      })
      all_counts <- as.matrix(do.call(rbind, counts_intensityNorm))
      global_median <- median(all_counts[,-(allICs-2)], na.rm=TRUE)
      counts_IPC_global_rescale <- list()
      for (i in 1:length(counts)){
        global_scaling <- global_median / plate_medians[[i]]
        counts_IPC_global_rescale[[i]] <- counts_IPC[[i]] * global_scaling
      }
      normCountsList <- counts_IPC_global_rescale
    }
  } # end intensityNorm step
  normCounts <- as.matrix(do.call(rbind, normCountsList))
  normData <- do.call(rbind, countMatrices)
  normData[,3:ncol(normData)] <- normCounts
  normData[,3:ncol(normData)] <- apply(normData[,3:ncol(normData)], 2, as.numeric)
  plateNs <- unlist(lapply(countMatrices, nrow))
  plate <- NULL
  for (i in 1:length(plateNs)){
    plate <- c(plate, rep(i, plateNs[i]))
  }
  normData$plate <- plate
  normData <- normData[,c(1,2,ncol(normData),3:(ncol(normData)-1))]
  return(normData)
}
