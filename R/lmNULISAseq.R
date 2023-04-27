#' Linear regression model for NULISAseq differential expression test
#'
#' Fits linear regression model to each target in the NULISAseq data set. 
#' Outputs coefficients, t-statistics, unadjusted and adjusted p-values.
#'
#' @param data A matrix with targets in rows, samples in columns. 
#' Row names should be the target names, and column names are the sample names.
#' It is assumed that data has already been transformed
#' using \code{log2(x + 0.01)} for each NULISAseq normalized count value \code{x}.
#' @param sampleInfo A data frame with sample metadata. Rows are samples, 
#' columns are sample metadata variables. Differential expression analysis will 
#' only be done on the samples in sampleInfo, or a subset of those samples as 
#' specified using arguments \code{exclude_samples} or \code{sample_subset}. 
#' \code{sampleInfo} should have a column for each 
#' variable included in the linear regression models. String variables will 
#' be automatically treated as factors, and numeric variables will be 
#' treated as numeric.
#' @param sampleName_var The name of the column of sampleInfo that matches
#' the column names of \code{data}. This variable will be used to merge the 
#' target expression data with the sample metadata. 
#' @param modelFormula A string that represents the right hand side of the model 
#' formula (everything after the \code{~}) used for the linear model. For example \code{modelFormula = 
#' "disease + age + sex + plate"} test for differences in target expression 
#' by disease group, adjusted for age, sex, and plate. \code{modelFormula = 
#' "disease * age + sex + plate"} includes both main and interaction 
#' effects for disease and age. See \code{?lm()}.
#' @param reduced_modelFormula Optional reduced model formula 
#' that contains only a subset of the terms in modelFormula. 
#' The reduced model serves as null model for an F-test using \code{anova()}. 
#' This could be useful for testing the overall significance of factor 
#' variables with more than 2 levels. 
#' @param exclude_targets A vector of target names for targets that will be 
#' excluded from the differential expression analysis. Internal control targets, 
#' for example, should probably always be excluded.
#' @param exclude_samples A vector of sample names for samples that will be excluded
#' from differential expression analysis. External control wells (IPCs, NCs, SC,)
#' should usually be excluded.
#' @param target_subset Overrides exclude_targets. A vector of target names 
#' for targets that will be included in differential expression analysis.
#' @param sample_subset Overrides exclude_samples. A vector of sample names 
#' for samples that will be included in differential expression analysis.
#' @param return_model_fits Logical \code{TRUE} or \code{FALSE} (default).
#' Should a list of the model fits be returned? Might be useful for more 
#' detailed analyses and plotting. However, also requires using more memory.
#'
#' @return A list including the following:
#' \item{modelStats}{A data frame with rows corresponding to targets and columns 
#' corresponding to estimated model coefficients, unadjusted p-values, 
#' Bonferroni adjusted p-values, and Benjamini-Hochberg false discovery rate
#' adjusted p-values (see \code{?p.adjust()})}
#' \item{modelFits}{A list of length equal to number of targets containing
#' the model fit output from \code{lm()}. Only returned when 
#' \code{return_model_fits=TRUE}.}
#' \item{Fstats}{A data frame with rows corresponding to targets and columns }
#'
#' 
#'
#' @export
#'
lmNULISAseq <- function(data, 
                        sampleInfo,
                        sampleName_var,
                        modelFormula,
                        reduced_modelFormula=NULL,
                        exclude_targets=NULL,
                        exclude_samples=NULL,
                        target_subset=NULL,
                        sample_subset=NULL,
                        return_model_fits=FALSE){
  # get data target subset
  if(!is.null(exclude_targets) & is.null(target_subset)){
    data <- data[!(rownames(data) %in% exclude_targets),]
  } 
  if(!is.null(target_subset)){
    data <- data[target_subset,]
  }
  # get sample subset
  if(!is.null(exclude_samples) & is.null(sample_subset)){
    data <- data[,!(colnames(data) %in% exclude_samples)]
    sampleInfo <- sampleInfo[!(sampleInfo[,sampleName_var] %in% exclude_samples),]
  } 
  if(!is.null(sample_subset)){
    data <- data[,sample_subset]
    sampleInfo <- sampleInfo[(sampleInfo[,sampleName_var] %in% sample_subset),]
  }
  # define vector of targets
  targets <- rownames(data)
  # create empty objects to store results
  modelFits <- vector(mode='list', length=length(targets))
  names(modelFits) <- targets
  stats_list <- vector(mode='list', length=length(targets))
  names(stats_list) <- targets
  # loop over targets and fit model
  for(i in 1:length(targets)){
    target <- targets[i]
    target_data <- data.frame(sampleName=colnames(data),
                              target_data=data[target,])
    model_data <- merge(sampleInfo, target_data,
                        all.x=TRUE, all.y=FALSE,
                        by.x=sampleName_var, by.y='sampleName')
    model_formula <- as.formula(paste0('target_data ~ ', modelFormula))
    model_fit <- lm(model_formula, data=model_data)
    coef_table <- summary(model_fit)$coefficients
    coefs <- coef_table[2:nrow(coef_table),1]
    t_vals <- coef_table[2:nrow(coef_table),3]
    p_vals <- coef_table[2:nrow(coef_table),4]
    modelFits[[i]] <- model_fit
    stats_list[[i]] <- list(coefs=coefs,
                            t_vals=t_vals,
                            p_vals=p_vals)
  }
  # format output
  coef <- do.call(rbind, lapply(stats_list, function(x) x$coefs))
  t_val <- do.call(rbind, lapply(stats_list, function(x) x$t_vals))
  p_val <- do.call(rbind, lapply(stats_list, function(x) x$p_vals))
  p_val_FDR <- apply(p_val, 2, p.adjust, method='BH')
  p_val_bonf <- apply(p_val, 2, p.adjust, method='bonferroni')
  colnames(coef) <- paste0(colnames(coef), '_coef')
  colnames(t_val) <- paste0(colnames(t_val), '_tstat')
  colnames(p_val) <- paste0(colnames(p_val), '_pval_unadj')
  colnames(p_val_FDR) <- paste0(colnames(p_val_FDR), '_pval_FDR')
  colnames(p_val_bonf) <- paste0(colnames(p_val_bonf), '_pval_bonf')
  column_order <- c(rbind(colnames(coef), colnames(t_val), colnames(p_val), colnames(p_val_FDR), colnames(p_val_bonf)))
  modelStats <- cbind(coef, t_val, p_val, p_val_FDR, p_val_bonf)[,column_order]
  modelStats <- data.frame(target=rownames(modelStats), modelStats)
  # do Ftest if specified
  if(!is.null(reduced_modelFormula)){
    Fstats_list <- vector(mode='list', length=length(targets))
    names(Fstats_list) <- targets
    for(i in 1:length(targets)){
      target <- targets[i]
      target_data <- data.frame(sampleName=colnames(data),
                                target_data=data[target,])
      model_data <- merge(sampleInfo, target_data,
                          all.x=TRUE, all.y=FALSE,
                          by.x=sampleName_var, by.y='sampleName')
      model_formula <- as.formula(paste0('target_data ~ ', reduced_modelFormula))
      model_fit <- lm(model_formula, data=model_data)
      anova_test <- anova(model_fit, modelFits[[i]])
      Fstats_list[[i]] <- c(Fstat=anova_test$F[2], 
                            Df=anova_test$Df[2],
                            Ftest_pval=anova_test$`Pr(>F)`[2])
    }
    # format output
    Fstats <- do.call(rbind, Fstats_list)
    Fstats <- data.frame(target=rownames(Fstats), Fstats)
    Fstats$Ftest_pval_FDR <- p.adjust(Fstats$Ftest_pval, method='BH')
    Fstats$Ftest_pval_bonf <- p.adjust(Fstats$Ftest_pval, method='bonferroni')
  }
  if(return_model_fits==FALSE & is.null(reduced_modelFormula)){
    output <- list(modelStats=modelStats)
  } else if(return_model_fits==TRUE & is.null(reduced_modelFormula)){
    output <- list(modelStats=modelStats,
                   modelFits=modelFits)
  } else if(return_model_fits==FALSE & !is.null(reduced_modelFormula)){
    output <- list(modelStats=modelStats,
                   Fstats=Fstats)
  } else if(return_model_fits==TRUE & !is.null(reduced_modelFormula)){
    output <- list(modelStats=modelStats,
                   Fstats=Fstats,
                   modelFits=modelFits)
  }
  return(output)
}
