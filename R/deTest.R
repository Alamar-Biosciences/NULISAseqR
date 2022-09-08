#' Differential Expression Test for all Targets.
#'
#' Within or between-sample t-test.
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
withinSampleTest <- function(interPlateNormData){
  techReps <- techRepsD1[techRepsD1!=100]
  group <- groupD[groupD!=100]
  plateD1 <- log((interPlateNormData[interPlateNormData$plate==3,4:53] + 1), base=2)
  plateD2 <- log((interPlateNormData[interPlateNormData$plate==4,4:53] + 1), base=2)
  plateD1_0x_vs_1x <- matrix(nrow=30, ncol=50)
  plateD1_0x_vs_2x <- matrix(nrow=30, ncol=50)
  plateD1_1x_vs_2x <- matrix(nrow=30, ncol=50)
  plateD2_0x_vs_1x <- matrix(nrow=30, ncol=50)
  plateD2_0x_vs_2x <- matrix(nrow=30, ncol=50)
  plateD2_1x_vs_2x <- matrix(nrow=30, ncol=50)
  for (i in 1:30){
    plateD1_0x_vs_1x[i,] <- unlist(plateD1[techReps==i & group==1,] - plateD1[techReps==i & group==0,])
    plateD1_0x_vs_2x[i,] <- unlist(plateD1[techReps==i & group==2,] - plateD1[techReps==i & group==0,])
    plateD1_1x_vs_2x[i,] <- unlist(plateD1[techReps==i & group==2,] - plateD1[techReps==i & group==1,])
    plateD2_0x_vs_1x[i,] <- unlist(plateD2[techReps==i & group==1,] - plateD2[techReps==i & group==0,])
    plateD2_0x_vs_2x[i,] <- unlist(plateD2[techReps==i & group==2,] - plateD2[techReps==i & group==0,])
    plateD2_1x_vs_2x[i,] <- unlist(plateD2[techReps==i & group==2,] - plateD2[techReps==i & group==1,])
  }
  tResults_D1_0x_vs_1x <- data.frame(matrix(nrow=50, ncol=2))
  tResults_D1_0x_vs_2x <- data.frame(matrix(nrow=50, ncol=2))
  tResults_D1_1x_vs_2x <- data.frame(matrix(nrow=50, ncol=2))
  tResults_D2_0x_vs_1x <- data.frame(matrix(nrow=50, ncol=2))
  tResults_D2_0x_vs_2x <- data.frame(matrix(nrow=50, ncol=2))
  tResults_D2_1x_vs_2x <- data.frame(matrix(nrow=50, ncol=2))
  colnames(tResults_D1_0x_vs_1x) <- c('effectSize', 'pval')
  colnames(tResults_D1_0x_vs_2x) <- c('effectSize', 'pval')
  colnames(tResults_D1_1x_vs_2x) <- c('effectSize', 'pval')
  colnames(tResults_D2_0x_vs_1x) <- c('effectSize', 'pval')
  colnames(tResults_D2_0x_vs_2x) <- c('effectSize', 'pval')
  colnames(tResults_D2_1x_vs_2x) <- c('effectSize', 'pval')
  for (i in 1:50){
    # do t tests
    t_D1_0x_vs_1x <- t.test(plateD1_0x_vs_1x[,i])
    t_D1_0x_vs_2x <- t.test(plateD1_0x_vs_2x[,i])
    t_D1_1x_vs_2x <- t.test(plateD1_1x_vs_2x[,i])
    t_D2_0x_vs_1x <- t.test(plateD2_0x_vs_1x[,i])
    t_D2_0x_vs_2x <- t.test(plateD2_0x_vs_2x[,i])
    t_D2_1x_vs_2x <- t.test(plateD2_1x_vs_2x[,i])
    # save results
    tResults_D1_0x_vs_1x[i,] <- c(t_D1_0x_vs_1x$estimate, t_D1_0x_vs_1x$p.value)
    tResults_D1_0x_vs_2x[i,] <- c(t_D1_0x_vs_2x$estimate, t_D1_0x_vs_2x$p.value)
    tResults_D1_1x_vs_2x[i,] <- c(t_D1_1x_vs_2x$estimate, t_D1_1x_vs_2x$p.value)
    tResults_D2_0x_vs_1x[i,] <- c(t_D2_0x_vs_1x$estimate, t_D2_0x_vs_1x$p.value)
    tResults_D2_0x_vs_2x[i,] <- c(t_D2_0x_vs_2x$estimate, t_D2_0x_vs_2x$p.value)
    tResults_D2_1x_vs_2x[i,] <- c(t_D2_1x_vs_2x$estimate, t_D2_1x_vs_2x$p.value)
  }
  # do FDR adjustment
  tResults_D1_0x_vs_1x$pvalFDR <- p.adjust(tResults_D1_0x_vs_1x$pval, method='fdr')
  tResults_D1_0x_vs_2x$pvalFDR <- p.adjust(tResults_D1_0x_vs_2x$pval, method='fdr')
  tResults_D1_1x_vs_2x$pvalFDR <- p.adjust(tResults_D1_1x_vs_2x$pval, method='fdr')
  tResults_D2_0x_vs_1x$pvalFDR <- p.adjust(tResults_D2_0x_vs_1x$pval, method='fdr')
  tResults_D2_0x_vs_2x$pvalFDR <- p.adjust(tResults_D2_0x_vs_2x$pval, method='fdr')
  tResults_D2_1x_vs_2x$pvalFDR <- p.adjust(tResults_D2_1x_vs_2x$pval, method='fdr')
  # define row names
  rownames(tResults_D1_0x_vs_1x) <- colnames(plateD1)
  rownames(tResults_D1_0x_vs_2x) <- colnames(plateD1)
  rownames(tResults_D1_1x_vs_2x) <- colnames(plateD1)
  rownames(tResults_D2_0x_vs_1x) <- colnames(plateD1)
  rownames(tResults_D2_0x_vs_2x) <- colnames(plateD1)
  rownames(tResults_D2_1x_vs_2x) <- colnames(plateD1)
  return(list(tResults_D1_0x_vs_1x=tResults_D1_0x_vs_1x,
              tResults_D1_0x_vs_2x=tResults_D1_0x_vs_2x,
              tResults_D1_1x_vs_2x=tResults_D1_1x_vs_2x,
              tResults_D2_0x_vs_1x=tResults_D2_0x_vs_1x,
              tResults_D2_0x_vs_2x=tResults_D2_0x_vs_2x,
              tResults_D2_1x_vs_2x=tResults_D2_1x_vs_2x))
}

withinSampleBetweenPlateTest <- function(interPlateNormData){
  techReps <- techRepsD1[techRepsD1!=100]
  group <- groupD[groupD!=100]
  plateD1 <- log((interPlateNormData[interPlateNormData$plate==3,4:53] + 1), base=2)
  plateD2 <- log((interPlateNormData[interPlateNormData$plate==4,4:53] + 1), base=2)
  plateD1_D2_0x_vs_1x <- matrix(nrow=30, ncol=50)
  plateD1_D2_0x_vs_2x <- matrix(nrow=30, ncol=50)
  plateD1_D2_1x_vs_2x <- matrix(nrow=30, ncol=50)
  plateD2_D1_0x_vs_1x <- matrix(nrow=30, ncol=50)
  plateD2_D1_0x_vs_2x <- matrix(nrow=30, ncol=50)
  plateD2_D1_1x_vs_2x <- matrix(nrow=30, ncol=50)
  for (i in 1:30){
    plateD1_D2_0x_vs_1x[i,] <- unlist(plateD1[techReps==i & group==1,] - plateD2[techReps==i & group==0,])
    plateD1_D2_0x_vs_2x[i,] <- unlist(plateD1[techReps==i & group==2,] - plateD2[techReps==i & group==0,])
    plateD1_D2_1x_vs_2x[i,] <- unlist(plateD1[techReps==i & group==2,] - plateD2[techReps==i & group==1,])
    plateD2_D1_0x_vs_1x[i,] <- unlist(plateD2[techReps==i & group==1,] - plateD1[techReps==i & group==0,])
    plateD2_D1_0x_vs_2x[i,] <- unlist(plateD2[techReps==i & group==2,] - plateD1[techReps==i & group==0,])
    plateD2_D1_1x_vs_2x[i,] <- unlist(plateD2[techReps==i & group==2,] - plateD1[techReps==i & group==1,])
  }
  tResults_D1_D2_0x_vs_1x <- data.frame(matrix(nrow=50, ncol=2))
  tResults_D1_D2_0x_vs_2x <- data.frame(matrix(nrow=50, ncol=2))
  tResults_D1_D2_1x_vs_2x <- data.frame(matrix(nrow=50, ncol=2))
  tResults_D2_D1_0x_vs_1x <- data.frame(matrix(nrow=50, ncol=2))
  tResults_D2_D1_0x_vs_2x <- data.frame(matrix(nrow=50, ncol=2))
  tResults_D2_D1_1x_vs_2x <- data.frame(matrix(nrow=50, ncol=2))
  colnames(tResults_D1_D2_0x_vs_1x) <- c('effectSize', 'pval')
  colnames(tResults_D1_D2_0x_vs_2x) <- c('effectSize', 'pval')
  colnames(tResults_D1_D2_1x_vs_2x) <- c('effectSize', 'pval')
  colnames(tResults_D2_D1_0x_vs_1x) <- c('effectSize', 'pval')
  colnames(tResults_D2_D1_0x_vs_2x) <- c('effectSize', 'pval')
  colnames(tResults_D2_D1_1x_vs_2x) <- c('effectSize', 'pval')
  for (i in 1:50){
    # do t tests
    t_D1_D2_0x_vs_1x <- t.test(plateD1_D2_0x_vs_1x[,i])
    t_D1_D2_0x_vs_2x <- t.test(plateD1_D2_0x_vs_2x[,i])
    t_D1_D2_1x_vs_2x <- t.test(plateD1_D2_1x_vs_2x[,i])
    t_D2_D1_0x_vs_1x <- t.test(plateD2_D1_0x_vs_1x[,i])
    t_D2_D1_0x_vs_2x <- t.test(plateD2_D1_0x_vs_2x[,i])
    t_D2_D1_1x_vs_2x <- t.test(plateD2_D1_1x_vs_2x[,i])
    # save results
    tResults_D1_D2_0x_vs_1x[i,] <- c(t_D1_D2_0x_vs_1x$estimate, t_D1_D2_0x_vs_1x$p.value)
    tResults_D1_D2_0x_vs_2x[i,] <- c(t_D1_D2_0x_vs_2x$estimate, t_D1_D2_0x_vs_2x$p.value)
    tResults_D1_D2_1x_vs_2x[i,] <- c(t_D1_D2_1x_vs_2x$estimate, t_D1_D2_1x_vs_2x$p.value)
    tResults_D2_D1_0x_vs_1x[i,] <- c(t_D2_D1_0x_vs_1x$estimate, t_D2_D1_0x_vs_1x$p.value)
    tResults_D2_D1_0x_vs_2x[i,] <- c(t_D2_D1_0x_vs_2x$estimate, t_D2_D1_0x_vs_2x$p.value)
    tResults_D2_D1_1x_vs_2x[i,] <- c(t_D2_D1_1x_vs_2x$estimate, t_D2_D1_1x_vs_2x$p.value)
  }
  # do FDR adjustment
  tResults_D1_D2_0x_vs_1x$pvalFDR <- p.adjust(tResults_D1_D2_0x_vs_1x$pval, method='fdr')
  tResults_D1_D2_0x_vs_2x$pvalFDR <- p.adjust(tResults_D1_D2_0x_vs_2x$pval, method='fdr')
  tResults_D1_D2_1x_vs_2x$pvalFDR <- p.adjust(tResults_D1_D2_1x_vs_2x$pval, method='fdr')
  tResults_D2_D1_0x_vs_1x$pvalFDR <- p.adjust(tResults_D2_D1_0x_vs_1x$pval, method='fdr')
  tResults_D2_D1_0x_vs_2x$pvalFDR <- p.adjust(tResults_D2_D1_0x_vs_2x$pval, method='fdr')
  tResults_D2_D1_1x_vs_2x$pvalFDR <- p.adjust(tResults_D2_D1_1x_vs_2x$pval, method='fdr')
  # define row names
  rownames(tResults_D1_D2_0x_vs_1x) <- colnames(plateD1)
  rownames(tResults_D1_D2_0x_vs_2x) <- colnames(plateD1)
  rownames(tResults_D1_D2_1x_vs_2x) <- colnames(plateD1)
  rownames(tResults_D2_D1_0x_vs_1x) <- colnames(plateD1)
  rownames(tResults_D2_D1_0x_vs_2x) <- colnames(plateD1)
  rownames(tResults_D2_D1_1x_vs_2x) <- colnames(plateD1)
  return(list(tResults_D1_D2_0x_vs_1x=tResults_D1_D2_0x_vs_1x,
              tResults_D1_D2_0x_vs_2x=tResults_D1_D2_0x_vs_2x,
              tResults_D1_D2_1x_vs_2x=tResults_D1_D2_1x_vs_2x,
              tResults_D2_D1_0x_vs_1x=tResults_D2_D1_0x_vs_1x,
              tResults_D2_D1_0x_vs_2x=tResults_D2_D1_0x_vs_2x,
              tResults_D2_D1_1x_vs_2x=tResults_D2_D1_1x_vs_2x))
}
betweenSampleTest <- function(interPlateNormData, 
                              effectSize, 
                              group1samples,
                              group2samples,
                              group1plate,
                              group2plate){
  techReps <- techRepsD1[techRepsD1!=100]
  group <- groupD[groupD!=100]
  group1plateData <- interPlateNormData[interPlateNormData$plate==group1plate,]
  group1data <- group1plateData[techReps %in% group1samples & group==0,4:53]
  group2plateData <- interPlateNormData[interPlateNormData$plate==group1plate,]
  group2data <- group2plateData[techReps %in% group2samples & group==effectSize,4:53]
  group1data <- log((group1data + 1), base=2)
  group2data <- log((group2data + 1), base=2)
  group1data$group <- rep(0, 15)
  group2data$group <- rep(1, 15)
  group1data$plate <- as.factor(rep(group1plate, 15))
  group2data$plate <- as.factor(rep(group2plate, 15))
  allData <- rbind(group1data, group2data)
  # if plates differ
  # fit linear model to each target
  if(group1plate != group2plate){
    lmResults <- data.frame(matrix(nrow=50, ncol=2))
    colnames(lmResults) <- c('groupCoef', 'groupPval')
    for (i in 1:50){
      fit <- lm(allData[,i] ~ as.factor(allData$group))
      lmResults[i,] <- c(coef(fit)[2],
                         summary(fit)$coefficients[2,4])
    }
    lmResults$groupPvalFDR <- p.adjust(lmResults$groupPval, method='fdr')
    rownames(lmResults) <- colnames(allData)[1:50]
    return(lmResults)
  } else if(group1plate==group2plate){
    tResults <- data.frame(matrix(nrow=50, ncol=2))
    colnames(tResults) <- c('groupDiff', 'groupPval')
    for (i in 1:50){
      t <- t.test(allData[,i] ~ allData$group)
      tResults[i,] <- c(t$estimate[2] - t$estimate[1],
                        t$p.value)
    }
    tResults$groupPvalFDR <- p.adjust(tResults$groupPval, method='fdr')
    rownames(tResults) <- colnames(allData)[1:50]
    return(tResults)
  }
}