#' Calculate Inter-plate Coefficients of Variation
#'
#' Calculates CV% for the specified replicates across multiple plates. 
#'
#' @param data_list A list of data matrices that contain data for 
#' the sample replicates for calculating inter-plate CV.
#' @param samples A list of vectors of length equal to 
#' the number of columns in the corresponding data matrix which encodes
#' the sets of sample replicates that will be used to calculate inter-plate CVs. 
#' Can be a list of either numeric or character vectors, where matching entries 
#' represent replicates for a given sample. 
#' Any NAs in the sample vectors will be excluded (it may be desirable to 
#' exclude negative controls or inter-plate controls, for example). 
#' Unique entries of these vectors will be used to identify 
#' unique sample replicate sets in the results output. 
#' @param aboveLOD A list of logical TRUE/FALSE matrices output by the 
#' NULISAseqR::lod() function. Indicates which values to include / exclude 
#' in the CV calculation. 
#' Default is NULL, which includes all values in CV calculation.
#' @param exclude_targets A list of vectors of row names or row indices 
#' representing the 
#' targets that should be excluded from intra-plate CV calculation. For example,
#' one might want to exclude internal controls. Default is NULL, which includes
#' all targets in the data_matrix.
#' @param useMean A logical TRUE / FALSE that indicates whether or not 
#' to calculate inter-CV using the intra-plate means of replicates (TRUE), or 
#' to use the pooled replicate data (FALSE). 
#'
#' @return matrix with inter-plate %CV values for each target for each sample.
#'
#' @examples
#' 
#'
#' @export
#'
interCV <- function(norm_data,
                         techReps,
                         techRepOrder_alpha,
                         techRepOrder_kf){
  CV <- matrix(nrow=21, ncol=49)
  colnames(CV) <- c('sample', 'alpha_rep', 'kf_rep',
                    colnames(norm_data[[1]])[3:48])
  CV[,1] <- rep(1:7, each=3)
  CV[,2] <- c(t(techRepOrder_alpha))
  CV[,3] <- c(t(techRepOrder_kf))
  for(i in 1:21){
    techRep_i <- CV[i,1]
    alpha_row <- which(techReps==techRep_i)[CV[i,2]]
    kf_row <- which(techReps==techRep_i)[CV[i,3]]
    data_i <- rbind(
      norm_data[[1]][alpha_row,3:48],
      norm_data[[2]][kf_row,3:48]
    )
    means_i <- colMeans(data_i)
    sds_i <- apply(data_i, 2, sd)
    cvs_i <- sds_i/means_i*100
    CV[i,4:ncol(CV)] <- cvs_i
  }
  # take mean CV over combos for each techRep
  avgCV <- matrix(nrow=7, ncol=47)
  colnames(avgCV) <- c('techRep', colnames(norm_data[[1]])[3:48])
  avgCV[,1] <- 1:7
  for(i in 1:7){
    CV_data_i <- CV[CV[,1]==i,4:ncol(CV)]
    avgCV[i,2:ncol(avgCV)] <- colMeans(CV_data_i, na.rm=TRUE)
  }
  avgCV[!is.finite(avgCV)] <- NA
  return(list(CV=CV, avgCV=avgCV))
}

CV <- function(interPlateNormData,
               techReps){
  CV <- matrix(nrow=30, ncol=40)
  colnames(CV) <- colnames(interPlateNormData)[4:43]
  for(i in 1:30){
    data_i <- interPlateNormData[techReps==i,4:43]
    means_i <- colMeans(data_i)
    sds_i <- apply(data_i, 2, sd)
    cvs_i <- sds_i/means_i*100
    CV[i,] <- cvs_i
  }
  meanCV <- colMeans(CV, na.rm=TRUE)
  medianCV <- apply(CV, 2, median, na.rm=TRUE)
  return(list(CV=CV,
              meanCV=meanCV,
              medianCV=medianCV))
}