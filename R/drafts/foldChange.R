#' Observed vs Expected Fold Change 
#'
#' Calculates fold change across targets for pairs of samples. 
#' Option to calculate error from expected fold change values 
#' if those values are provided.
#'
#' @param 
#
#'
#' @return Draws a plot.
#'
#' @examples
#' 
#'
#' @export
#' 
# NEEDS EDITING
FC <- function(interPlateNormData){
  FC2 <- matrix(nrow=60, ncol=12)
  FC2[,1] <- rep(3:4, each=30)
  plate <- interPlateNormData[,'plate']
  techRepsD1_spikein <- 1:30
  techRepsD2_spikein <- 1:30
  FC2[,2] <- c(techRepsD1_spikein, techRepsD2_spikein)
  techRepsD1_spikein_ordered <- techRepsD1[techRepsD1 %in% techRepsD1_spikein]
  techRepsD2_spikein_ordered <- techRepsD2[techRepsD2 %in% techRepsD2_spikein]
  groupD_spikein <- groupD[techRepsD2 %in% techRepsD2_spikein]
  plateD1 <- interPlateNormData[plate==3,44:53][techRepsD2 %in% techRepsD2_spikein,]
  plateD2 <- interPlateNormData[plate==4,44:53][techRepsD2 %in% techRepsD2_spikein,]
  # replace any zeros with 1
  # plateD1[plateD1==0] <- 1
  # plateD2[plateD2==0] <- 1
  plateD1 <- cbind(techRepsD1_spikein_ordered, groupD_spikein, plateD1)
  plateD2 <- cbind(techRepsD2_spikein_ordered, groupD_spikein, plateD2)
  colnames(FC2) <- c('plate', 'techRep', colnames(plateD1)[3:12])
  for(i in 1:30){
    plateD1_techRep_i <- techRepsD1_spikein[i]
    plateD2_techRep_i <- techRepsD2_spikein[i]
    plateD1_data_i <- plateD1[plateD1[,1]==plateD1_techRep_i,]
    plateD2_data_i <- plateD2[plateD2[,1]==plateD2_techRep_i,]
    # plate D1
    plateD1_FC2 <- (plateD1_data_i[plateD1_data_i$groupD_spikein==2,3:12] - plateD1_data_i[plateD1_data_i$groupD_spikein==0,3:12])/
      (plateD1_data_i[plateD1_data_i$groupD_spikein==1,3:12] - plateD1_data_i[plateD1_data_i$groupD_spikein==0,3:12])
    # plate D2
    plateD2_FC2 <- (plateD2_data_i[plateD2_data_i$groupD_spikein==2,3:12] - plateD2_data_i[plateD2_data_i$groupD_spikein==0,3:12])/
      (plateD2_data_i[plateD2_data_i$groupD_spikein==1,3:12] - plateD2_data_i[plateD2_data_i$groupD_spikein==0,3:12])
    FC2[i,3:12] <- unlist(plateD1_FC2)
    FC2[i+30,3:12] <- unlist(plateD2_FC2)
  }
  # additional function code
  trimmed_MAE <- function(x){
    x05 <- quantile(x, 0.05, na.rm=TRUE)
    x95 <- quantile(x, 0.95, na.rm=TRUE)
    x90 <- x[(x >= x05) & (x <= x95)]
    mean(abs(x90-2), na.rm=TRUE)}
  
  MAE <- apply(FC2, 2, function(x) {
    x05 <- quantile(x, 0.05, na.rm=TRUE)
    x95 <- quantile(x, 0.95, na.rm=TRUE)
    x90 <- x[(x >= x05) & (x <= x95)]
    mean(abs(x90-2), na.rm=TRUE)})
  
  interPlateFC <- function(interPlateNormData){
    FC2 <- matrix(nrow=60, ncol=12)
    FC2[,1] <- rep(3:4, each=30)
    plate <- interPlateNormData[,'plate']
    techRepsD1_spikein <- 1:30
    techRepsD2_spikein <- 1:30
    FC2[,2] <- c(techRepsD1_spikein, techRepsD2_spikein)
    techRepsD1_spikein_ordered <- techRepsD1[techRepsD1 %in% techRepsD1_spikein]
    techRepsD2_spikein_ordered <- techRepsD2[techRepsD2 %in% techRepsD2_spikein]
    groupD_spikein <- groupD[techRepsD2 %in% techRepsD2_spikein]
    plateD1 <- interPlateNormData[plate==3,44:53][techRepsD2 %in% techRepsD2_spikein,]
    plateD2 <- interPlateNormData[plate==4,44:53][techRepsD2 %in% techRepsD2_spikein,]
    # replace any zeros with 1
    # plateD1[plateD1==0] <- 1
    # plateD2[plateD2==0] <- 1
    plateD1 <- cbind(techRepsD1_spikein_ordered, groupD_spikein, plateD1)
    plateD2 <- cbind(techRepsD2_spikein_ordered, groupD_spikein, plateD2)
    colnames(FC2) <- c('plate', 'techRep', colnames(plateD1)[3:12])
    for(i in 1:30){
      plateD1_techRep_i <- techRepsD1_spikein[i]
      plateD2_techRep_i <- techRepsD2_spikein[i]
      plateD1_data_i <- plateD1[plateD1[,1]==plateD1_techRep_i,]
      plateD2_data_i <- plateD2[plateD2[,1]==plateD2_techRep_i,]
      # plate D1
      plateD1_delta2 <- (plateD1_data_i[plateD1_data_i$groupD_spikein==2,3:12] - plateD1_data_i[plateD1_data_i$groupD_spikein==0,3:12])
      plateD1_delta1 <- (plateD1_data_i[plateD1_data_i$groupD_spikein==1,3:12] - plateD1_data_i[plateD1_data_i$groupD_spikein==0,3:12])
      # plate D2
      plateD2_delta2 <- (plateD2_data_i[plateD2_data_i$groupD_spikein==2,3:12] - plateD2_data_i[plateD2_data_i$groupD_spikein==0,3:12])
      plateD2_delta1 <- (plateD2_data_i[plateD2_data_i$groupD_spikein==1,3:12] - plateD2_data_i[plateD2_data_i$groupD_spikein==0,3:12])
      # calculate interPlate FCs
      plateD1_FC2 <- plateD1_delta2/plateD2_delta1
      plateD2_FC2 <- plateD2_delta2/plateD1_delta1
      FC2[i,3:12] <- unlist(plateD1_FC2)
      FC2[i+30,3:12] <- unlist(plateD2_FC2)
    }
    return(list(FC2=FC2))
  }
  return(list(FC2=FC2))
}
