#' Linear Mixed Effect Model for NULISAseq Data - 
#' targets as outcome
#'
#' Fits linear mixed effects model to each target in the NULISAseq data set,
#' using univariate targets as outcome in the model. 
#' Outputs coefficients, t-statistics, unadjusted and adjusted p-values.
#' Uses \code{lme4} and \code{lmerTest} packages.
#'
#' @param data A matrix of normalized NULISAseq data
#' with targets in rows, samples in columns. 
#' Row names should be the target names, and column names are the sample names.
#' It is assumed that data has already been transformed
#' using \code{log2(x + 1)} for each NULISAseq normalized count value \code{x}.
#' @param sampleInfo A data frame with sample metadata. Rows are samples, 
#' columns are sample metadata variables. Linear mixed effect models will 
#' only be done on the samples in sampleInfo, or a subset of those samples as 
#' specified using arguments \code{exclude_samples} or \code{sample_subset}. 
#' \code{sampleInfo} should have a column for each 
#' variable included in the linear mixed effect models. String variables will 
#' be automatically treated as factors, and numeric variables will be 
#' treated as numeric.
#' @param sampleName_var The name of the column of sampleInfo that matches
#' the column names of \code{data}. This variable will be used to merge the 
#' target expression data with the sample metadata. 
#' @param modelFormula_fixed A string that represents the fixed effects part of the model 
#' formula on the used for the linear mixed effects model. 
#' For example \code{modelFormula_fixed = 
#' "disease + age + sex + plate"} tests for differences in target expression 
#' by disease group, adjusted for age, sex, and plate. \code{modelFormula = 
#' "disease * age + sex + plate"} includes both main and interaction 
#' effects for disease and age. See \code{?lmer()}.
#' @param modelFormula_random A string that represents the random effects part of the model 
#' formula on the used for the linear mixed effects model. 
#' For example \code{modelFormula_random = "(1|participant_ID)"} 
#' creates a subject specific random intercept, where the variable 
#' \code{participant_ID} (a column in \code{sampleInfo} data frame) denotes 
#' repeated measures on the same subject. For subject-specific random intercept and 
#' slopes (not recommended when time is categorical), 
#' use \code{modelFormula_random = "(1 + time|participant_ID)"}. For random 
#' subject nested within plate (which may be useful when analyzing
#' a large number of plates together), use \code{modelFormula_random = 
#' "(1|plate_ID:participant_ID)"}.
#' See \code{?lmer()}.
#' @param reduced_modelFormula_fixed Optional reduced model formula 
#' for fixed effects that contains only a subset of the terms in modelFormula. 
#' This could be an empty string if the full model contains only one term.
#' The reduced model serves as null model for a likelihood ratio test 
#' (LRT, which is a Chi-square test) using \code{anova()}. 
#' This could be useful for testing the overall significance of factor 
#' variables with more than 2 levels, for example, testing the overall significance 
#' of a categorical time effect. The reduced model uses the same random effects 
#' as specified in \code{modelFormula_random}.
#' @param reduced_modelFormula_random Optional reduced random effects formula. 
#' If not specified, the reduced model will use the same random effects structure 
#' as the full model. Specifying this allows testing the significance of random effects
#' components. For example, to test if participant random effects are needed, you could
#' specify \code{reduced_modelFormula_random = "(1|plate_ID)"} when the full model has
#' \code{modelFormula_random = "(1|plate_ID:participant_ID)"}.
#' @param exclude_targets A vector of target names for targets that will be 
#' excluded from the linear mixed effect models. Internal control targets, 
#' for example, should probably always be excluded.
#' @param exclude_samples A vector of sample names for samples that will be excluded
#' from the linear mixed effect models. External control wells (IPCs, NCs, SC,)
#' should usually be excluded.
#' @param target_subset Overrides exclude_targets. A vector of target names 
#' for targets that will be included in the linear mixed effect models.
#' @param sample_subset Overrides exclude_samples. A vector of sample names 
#' for samples that will be included in the linear mixed effect models.
#' @param return_model_fits Logical \code{TRUE} or \code{FALSE} (default).
#' Should a list of the model fits be returned? Might be useful for more 
#' detailed analyses and plotting. However, also requires using more memory.
#'
#' @return A list including the following:
#' \item{modelStats}{A data frame with rows corresponding to targets and columns 
#' corresponding to estimated model coefficients, unadjusted p-values, 
#' Bonferroni adjusted p-values, and Benjamini-Hochberg false discovery rate
#' adjusted p-values (see \code{?p.adjust()}).}
#' \item{modelFits}{A list of length equal to number of targets containing
#' the model fit output from \code{lm()}. Only returned when 
#' \code{return_model_fits=TRUE}.}
#' \item{LRTstats}{A data frame with rows corresponding to targets and columns.}
#'
#' 
#'
#' @export
#'
lmerNULISAseq <- function(data, 
                          sampleInfo,
                          sampleName_var,
                          modelFormula_fixed,
                          modelFormula_random,
                          reduced_modelFormula_fixed=NULL,
                          reduced_modelFormula_random=NULL,
                          exclude_targets=NULL,
                          exclude_samples=NULL,
                          target_subset=NULL,
                          sample_subset=NULL,
                          return_model_fits=FALSE){
  # get data target subset
  if(!is.null(exclude_targets) & is.null(target_subset)){
    data <- data[!(rownames(data) %in% exclude_targets), , drop = FALSE]
  }
  if(!is.null(target_subset)){
    data <- data[target_subset, , drop = FALSE]
  }
  # get sample subset
  if(!is.null(exclude_samples) & is.null(sample_subset)){
    data <- data[,!(colnames(data) %in% exclude_samples), drop = FALSE]
    sampleInfo <- sampleInfo[!(sampleInfo[,sampleName_var] %in% exclude_samples), , drop = FALSE]
  }
  if(!is.null(sample_subset)){
    data <- data[,sample_subset, drop = FALSE]
    sampleInfo <- sampleInfo[(sampleInfo[,sampleName_var] %in% sample_subset), , drop = FALSE]
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
    tryCatch({
      target <- targets[i]
      target_data <- data.frame(sampleName=colnames(data),
                                target_data=unlist(data[target,]))
      model_data <- merge(sampleInfo, target_data,
                          all.x=TRUE, all.y=FALSE,
                          by.x=sampleName_var, by.y='sampleName')
      model_formula <- as.formula(paste0('target_data ~ ', 
                                         modelFormula_fixed, ' + ', 
                                         modelFormula_random))
      model_data <- model_data[complete.cases(model_data[,all.vars(model_formula)]),]
      model_fit <- lmerTest::lmer(model_formula, data=model_data)
      
      #convert coef table to a tibble/dataframe with metric column having all relevant rownames
      coef_table <- summary(model_fit)$coefficients %>%
        tibble::as_tibble(rownames = 'metric')
      
      coef_table$metric <- ifelse(coef_table$metric == '(Intercept)', 'intercept', coef_table$metric)
      # rename the `Estimate, t value, Pr(>|t|)` columns to required values
      idx <- which(names(coef_table) %in% c('Estimate', "t value", "Pr(>|t|)"))
      names(coef_table)[idx] <- c('coef', 't_vals', 'p_vals')
      
      # create the named vector for coef, t_vals and p_vals these values also include Intercept values
      coefs <- coef_table[, c('metric', 'coef')] %>% tibble::deframe()
      t_vals <- coef_table[, c('metric', 't_vals')] %>% tibble::deframe()
      p_vals <- coef_table[, c('metric', 'p_vals')] %>% tibble::deframe()
      
      modelFits[[i]] <- model_fit
      stats_list[[i]] <- list(coefs=coefs,
                              t_vals=t_vals,
                              p_vals=p_vals)
    }, error = function(e){
      warning(sprintf("%s Index: %d Target: %s - %s",
                      format(Sys.time()), i, targets[i], conditionMessage(e)),
              call. = FALSE)
    })
  }

  # Check if all models failed
  if (all(sapply(stats_list, is.null))) {
    stop("All targets failed model fitting. Common causes:\n",
         "  - A covariate has only one level after removing samples with missing data\n",
         "  - Insufficient samples remain after filtering.\n",
         "Check the error messages above for details.")
  }

  all_predictors <- unique(unlist(lapply(stats_list, function(x) names(x$coefs))))
  # format output
  coef <- safe_extract_matrix(stats_list, "coefs", all_predictors)
  t_val <- safe_extract_matrix(stats_list, "t_vals", all_predictors)
  p_val <- safe_extract_matrix(stats_list, "p_vals", all_predictors)

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
  # do likelihood ratio test (LRT) if specified
  if(!is.null(reduced_modelFormula_fixed)){
    LRTstats_list <- vector(mode='list', length=length(targets))
    names(LRTstats_list) <- targets
    for(i in 1:length(targets)){
      target <- targets[i]
      tryCatch({
        target_data <- data.frame(sampleName=colnames(data),
                                  target_data=unlist(data[target,]))
        model_data <- merge(sampleInfo, target_data,
                            all.x=TRUE, all.y=FALSE,
                            by.x=sampleName_var, by.y='sampleName')
        model_formula <- as.formula(paste0('target_data ~ ', 
                                           reduced_modelFormula_fixed, ' + ', 
                                           modelFormula_random))
        model_data <- model_data[complete.cases(model_data[,all.vars(formula(modelFits[[i]]))]),]
        model_fit <- lmerTest::lmer(model_formula, data=model_data)
        anova_test <- suppressMessages(anova(model_fit, modelFits[[i]]))
        LRTstats_list[[i]] <- c(Chisq_stat=anova_test$Chisq[2], 
                                Df=anova_test$Df[2],
                                Chisq_test_pval=anova_test$`Pr(>Chisq)`[2])
      }, error = function(e){
        warning(sprintf("%s Index: %d Target: %s - %s",
                        format(Sys.time()), i, targets[i], conditionMessage(e)),
                call. = FALSE)
      })
    }
    # format output
    LRTstats <- do.call(rbind, LRTstats_list)
    LRTstats <- data.frame(target=rownames(LRTstats), LRTstats)
    LRTstats$Chisq_test_pval_FDR <- p.adjust(LRTstats$Chisq_test_pval, method='BH')
    LRTstats$Chisq_test_pval_bonf <- p.adjust(LRTstats$Chisq_test_pval, method='bonferroni')
  }
  if(return_model_fits==FALSE & is.null(reduced_modelFormula_fixed)){
    output <- list(modelStats=modelStats)
  } else if(return_model_fits==TRUE & is.null(reduced_modelFormula_fixed)){
    output <- list(modelStats=modelStats,
                   modelFits=modelFits)
  } else if(return_model_fits==FALSE & !is.null(reduced_modelFormula_fixed)){
    output <- list(modelStats=modelStats,
                   LRTstats=LRTstats)
  } else if(return_model_fits==TRUE & !is.null(reduced_modelFormula_fixed)){
    output <- list(modelStats=modelStats,
                   LRTstats=LRTstats,
                   modelFits=modelFits)
  }
  return(output)
}
