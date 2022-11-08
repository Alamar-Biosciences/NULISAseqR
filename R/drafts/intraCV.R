#' Calculate Intra-plate Coefficients of Variation
#'
#' Calculates CV% for the specified replicates on a given plate. 
#'
#' @param data_matrix The Data matrix output from readNULISAseq.R
#' or normalized data from normalization functions.
#
#'
#' @return matrix with %CV values for each target for each 
#' multiple-replicate sample.
#'
#' @examples
#' 
#'
#' @export
#' 
intraCV <- function(countMatrix, samples){
  countMatrix <- countMatrix[-c(NIC, MCherry),]
  unique_samples <- unique(samples)[order(unique(samples))][c(16:21,1,7:14,2:6,15,22,28:35,23:27,36:44)]
  cv_results <- data.frame(matrix(nrow=nrow(countMatrix), ncol=length(unique_samples)))
  rownames(cv_results) <- countMatrix$BarcodeName
  colnames(cv_results) <- unique_samples
  for (i in 1:length(unique_samples)){
    sample_data <- countMatrix[,3:98][,samples==unique_samples[i]]
    sample_means <- rowMeans(sample_data)
    sample_sds <- apply(sample_data,1,sd)
    sample_cv <- sample_sds/sample_means*100
    cv_results[,i] <- sample_cv
  }
  return(cv_results)
}

intraPlateCV <- function(interPlateNormData,
                         techReps){
  CV_C1 <- matrix(nrow=30, ncol=41)
  CV_C2 <- matrix(nrow=30, ncol=41)
  CV_D1 <- matrix(nrow=30, ncol=41)
  CV_D2 <- matrix(nrow=30, ncol=41)
  CV_C1[,1] <- 1:30
  CV_C2[,1] <- 1:30
  CV_D1[,1] <- 1:30
  CV_D2[,1] <- 1:30
  plate <- interPlateNormData[,'plate']
  plateC1 <- interPlateNormData[plate==1,4:43][techRepsC1[techRepsC1!=100] %in% 1:30,]
  plateC2 <- interPlateNormData[plate==2,4:43][techRepsC2[techRepsC2!=100] %in% 1:30,]
  plateD1 <- interPlateNormData[plate==3,4:43][techRepsD1[techRepsD1!=100] %in% 1:30,]
  plateD2 <- interPlateNormData[plate==4,4:43][techRepsD2[techRepsD2!=100] %in% 1:30,]
  colnames(CV_C1) <- c('techRep', 
                       colnames(interPlateNormData)[4:43])
  colnames(CV_C2) <- colnames(CV_C1)
  colnames(CV_D1) <- c('techRep', 
                       colnames(interPlateNormData)[4:43])
  colnames(CV_D2) <- colnames(CV_D1)
  techRepsC1_intraCV <- techRepsC1[techRepsC1 %in% 1:30]
  techRepsC2_intraCV <- techRepsC2[techRepsC2 %in% 1:30]
  techRepsD1_intraCV <- techRepsD1[techRepsD1 %in% 1:30]
  techRepsD2_intraCV <- techRepsD2[techRepsD2 %in% 1:30]
  for(i in 1:30){
    techRep_i <- i
    # plate C1
    data_C1_i <- plateC1[techRepsC1_intraCV==techRep_i,]
    means_C1_i <- colMeans(data_C1_i)
    sds_C1_i <- apply(data_C1_i, 2, sd)
    cvs_C1_i <- sds_C1_i/means_C1_i*100
    CV_C1[i,2:41] <- cvs_C1_i
    # plate C2
    data_C2_i <- plateC2[techRepsC2_intraCV==techRep_i,]
    means_C2_i <- colMeans(data_C2_i)
    sds_C2_i <- apply(data_C2_i, 2, sd)
    cvs_C2_i <- sds_C2_i/means_C2_i*100
    CV_C2[i,2:41] <- cvs_C2_i
    # plate D1
    techRep_D1_i <- CV_D1[i,1]
    data_D1_i <- plateD1[techRepsD1_intraCV==techRep_D1_i,]
    means_D1_i <- colMeans(data_D1_i)
    sds_D1_i <- apply(data_D1_i, 2, sd)
    cvs_D1_i <- sds_D1_i/means_D1_i*100
    CV_D1[i,2:41] <- cvs_D1_i
    # plate D2
    techRep_D2_i <- CV_D2[i,1]
    data_D2_i <- plateD2[techRepsD2_intraCV==techRep_D2_i,]
    means_D2_i <- colMeans(data_D2_i)
    sds_D2_i <- apply(data_D2_i, 2, sd)
    cvs_D2_i <- sds_D2_i/means_D2_i*100
    CV_D2[i,2:41] <- cvs_D2_i
  }
  # calculate plate median CVs by target
  CV_medians <- matrix(nrow=4, ncol=41)
  colnames(CV_medians) <- c('plate', colnames(CV_C1)[2:41])
  CV_medians[,1] <- 1:4
  CV_medians[1,2:41] <- apply(CV_C1[,2:41], 2, median, na.rm=TRUE)
  CV_medians[2,2:41] <- apply(CV_C2[,2:41], 2, median, na.rm=TRUE)
  CV_medians[3,2:41] <- apply(CV_D1[,2:41], 2, median, na.rm=TRUE)
  CV_medians[4,2:41] <- apply(CV_D2[,2:41], 2, median, na.rm=TRUE)
  median_CV_medians <- apply(CV_medians[,2:41], 2, median, na.rm=TRUE)
  return(list(CV_C1=CV_C1,
              CV_C2=CV_C2,
              CV_D1=CV_D1,
              CV_D2=CV_D2,
              CV_medians=CV_medians,
              median_CV_medians=median_CV_medians))
}