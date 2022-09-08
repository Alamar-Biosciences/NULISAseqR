#' Calculate Inter-plate Coefficients of Variation
#'
#' Calculates CV% for the specified replicates across multiple plates. 
#'
#' @param plates A list of count data matrices. 
#
#'
#' @return matrix with %CV values for each target for each sample.
#'
#' @examples
#' 
#'
#' @export
#' 
# NEEDS EDITING
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