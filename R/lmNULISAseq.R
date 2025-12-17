#' Linear Regression Model for NULISAseq Data - 
#' targets as outcome
#'
#' Fits linear regression model to each target in the NULISAseq data set,
#' using univariate targets as outcome in the model. 
#' Outputs coefficients, t-statistics, unadjusted and adjusted p-values.
#'
#' @param data A matrix of normalized NULISAseq data with targets in rows, 
#' samples in columns. Row names should be the target names, and column names are the sample names. 
#' It is assumed that data has already been transformed
#' using \code{log2(x + 1)} for each NULISAseq normalized count value \code{x}.
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
#' from the differential expression analysis. External control wells (IPCs, NCs, SC,)
#' should usually be excluded.
#' @param target_subset Overrides exclude_targets. A vector of target names 
#' for targets that will be included in the differential expression analysis.
#' @param sample_subset Overrides exclude_samples. A vector of sample names 
#' for samples that will be included in the differential expression analysis.
#' @param return_model_fits Logical \code{TRUE} or \code{FALSE} (default).
#' Should a list of the model fits be returned? Might be useful for more
#' detailed analyses and plotting. However, also requires using more memory.
#' @param analysis_context Optional string to provide context in error messages.
#' Useful when calling lmNULISAseq from different analyses (e.g., "plate effect test").
#'
#' @return A list including the following:
#' \item{modelStats}{A data frame with rows corresponding to targets and columns 
#' corresponding to estimated model coefficients, unadjusted p-values, 
#' Bonferroni adjusted p-values, and Benjamini-Hochberg false discovery rate
#' adjusted p-values (see \code{?p.adjust()}).}
#' \item{modelFits}{A list of length equal to number of targets containing
#' the model fit output from \code{lm()}. Only returned when 
#' \code{return_model_fits=TRUE}.}
#' \item{Fstats}{A data frame with rows corresponding to targets and columns.}
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
                        return_model_fits=FALSE,
                        analysis_context=NULL){
  
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
      model_formula <- as.formula(paste0('target_data ~ ', modelFormula))
      model_data <- model_data[complete.cases(model_data[,all.vars(model_formula)]),]
      model_fit <- lm(model_formula, data=model_data)
      
      #convert coef table to a tibble/dataframe with metric column having all relevant rownames
      coef_table <- summary(model_fit)$coefficients %>%
        tibble::as_tibble(rownames = 'metric')
      
      # rename the `Estimate, t value, Pr(>|t|)` columns to required values
      idx <- which(names(coef_table) %in% c('Estimate', "t value", "Pr(>|t|)"))
      names(coef_table)[idx] <- c('coef', 't_vals', 'p_vals')
      
      coef_table$metric <- ifelse(coef_table$metric == '(Intercept)', 'intercept', coef_table$metric)
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
    context_msg <- if (!is.null(analysis_context)) {
      paste0(analysis_context, ": ")
    } else {
      ""
    }
    stop(context_msg,
         "All targets failed model fitting. Common causes:\n",
         "  - A covariate has only one level after removing samples with missing data\n",
         "  - Insufficient samples remain after filtering.\n",
         "Check the error messages above for details.",
         call. = FALSE)
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
  # do Ftest if specified
  if(!is.null(reduced_modelFormula)){
    Fstats_list <- vector(mode='list', length=length(targets))
    names(Fstats_list) <- targets
    for(i in 1:length(targets)){
      tryCatch({
      target <- targets[i]
      target_data <- data.frame(sampleName=colnames(data),
                                target_data=unlist(data[target,]))
      model_data <- merge(sampleInfo, target_data,
                          all.x=TRUE, all.y=FALSE,
                          by.x=sampleName_var, by.y='sampleName')
      model_formula <- as.formula(paste0('target_data ~ ', reduced_modelFormula))
      model_data <- model_data[complete.cases(model_data[,all.vars(formula(modelFits[[i]]))]),]
      model_fit <- lm(model_formula, data=model_data)
      anova_test <- anova(model_fit, modelFits[[i]])
      Fstats_list[[i]] <- c(Fstat=anova_test$F[2], 
                            Df=anova_test$Df[2],
                            Ftest_pval=anova_test$`Pr(>F)`[2])
      }, error = function(e){
        warning(sprintf("%s Index: %d Target: %s - %s",
                        format(Sys.time()), i, targets[i], conditionMessage(e)),
                call. = FALSE)
      })
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


#' Fill missing predictors with NA values
#'
#' Ensures a named vector contains all expected predictor names by filling in
#' missing predictors with NA. Used to standardize model coefficient vectors
#' across targets when some predictors may be absent due to singularities or
#' model fitting failures.
#'
#' @param x A named vector of model statistics (coefficients, t-values, or p-values)
#' @param all A character vector of all expected predictor names
#'
#' @return A named vector containing all predictors from 'all', with original
#'   values from 'x' where present and NA for missing predictors
#'
#' @export
fill_predictors <- function(x, all){
  fill <- sapply(setdiff(all, names(x)), function(x) NA)
  if(length(fill) !=0){
    c(x, fill)
  } else{
    x
  }
}

#' Safely extract and combine model statistics into a matrix
#'
#' Extracts a specific field (coefficients, t-values, or p-values) from a list
#' of model statistics and combines them into a matrix with targets as rows and
#' predictors as columns. Handles missing fields, NULL values, and targets with
#' different sets of predictors by filling missing values with NA.
#'
#' @param stats_list A named list where each element contains model statistics
#'   for one target. Each element should be a list with fields like 'coefs',
#'   't_vals', 'p_vals'
#' @param field_name A character string specifying which field to extract
#'   (e.g., 'coefs', 't_vals', 'p_vals')
#' @param all_predictors A character vector of all predictor names that should
#'   appear as columns in the output matrix
#'
#' @return A numeric matrix with targets as rows and predictors as columns,
#'   with NA values where predictors are missing for specific targets
#'
#' @export
safe_extract_matrix <- function(stats_list, field_name, all_predictors) {
  # Handle empty stats_list
  if (length(stats_list) == 0) {
    result <- matrix(numeric(0), nrow = 0, ncol = length(all_predictors))
    colnames(result) <- all_predictors
    return(result)
  }
  
  # Choose .id name that doesn't conflict with predictor names
  id_col_name <- if ("target" %in% all_predictors) ".target_id" else "target"
  
  result <- lapply(stats_list, function(x) {
    # Handle NULL or missing fields
    if (is.null(x) || is.null(x[[field_name]])) {
      return(setNames(rep(NA, length(all_predictors)), all_predictors))
    }
    # Apply fill_predictors
    fill_predictors(x[[field_name]], all_predictors)
  }) |>
    dplyr::bind_rows(.id = id_col_name) |>
    tibble::column_to_rownames(id_col_name) |>
    as.matrix()
  # Reorder columns to match all_predictors
  result[, all_predictors, drop = FALSE]
}