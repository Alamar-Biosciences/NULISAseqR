#' Linear Regression Model for NULISAseq Data - 
#' targets as predictors
#'
#' Fits linear regression model to each target in the NULISAseq data set,
#' using univariate targets as predictors in the model. 
#' Outputs coefficients, t-statistics, unadjusted and adjusted p-values.
#' This approach tests whether each 
#' target's expression is associated with the outcome variable, while 
#' adjusting for any specified covariates.
#' 
#' @param data A matrix of normalized NULISAseq data
#' with targets used as predictors in rows, samples in columns. 
#' Row names should be the target names, and column names are the sample names.
#' It is assumed that data has already been transformed
#' using \code{log2(x + 1)} for each NULISAseq normalized count value \code{x}.
#' @param sampleInfo A data frame with sample metadata including the response 
#' variable and covariates. Rows are samples, columns are sample metadata variables. 
#' Linear regression models will only be done on the samples in sampleInfo, or a subset of those samples as 
#' specified using arguments \code{exclude_samples} or \code{sample_subset}. 
#' \code{sampleInfo} should have a column for each 
#' variable included in the linear regression models. String variables will 
#' be automatically treated as factors, and numeric variables will be 
#' treated as numeric.
#' @param sampleName_var The name of the column of sampleInfo that matches
#' the column names of \code{data}. This variable will be used to merge the 
#' target expression data with the sample metadata.
#' @param response_var The name of the column of sampleInfo specifying the continuous numeric response 
#' variable.
#' @param modelFormula A string that represents the right hand side of the model formula.
#' The main effect of target expression 
#' will be automatically added as a predictor. Any interactions need to be specified 
#' in the model formula as "covariate * target". 
#' For example \code{modelFormula = "disease + age + sex + plate"} tests for variations of outcome explained by the 
#' target predictor, adjusted for disease group, age, sex, and plate. \code{modelFormula = 
#' "disease * target + age + sex + plate"} includes both main and interaction 
#' effects for disease and target expression. See \code{?lm()}.
#' @param reduced_modelFormula Optional reduced model formula 
#' that contains only a subset of the terms in modelFormula. 
#' The reduced model serves as null model for an F-test using \code{anova()}. 
#' This could be useful for testing the overall significance of factor 
#' variables with more than 2 levels. 
#' @param exclude_targets A vector of target names for targets that will be 
#' excluded from the linear regression models as predictors. Internal control targets, 
#' for example, should probably always be excluded.
#' @param exclude_samples A vector of sample names for samples that will be excluded
#' from the linear regression models. External control wells (IPCs, NCs, SC,)
#' should usually be excluded.
#' @param target_subset Overrides \code{exclude_targets}. A vector of target names 
#' for targets that will be included in the linear regression models as predictors.
#' @param sample_subset Overrides \code{exclude_samples}. A vector of sample names 
#' for samples that will be included in the linear regression models.
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
#' \item{Fstats}{A data frame with rows corresponding to targets and columns.}
#'
#' @importFrom magrittr %>%
#'
#' @export
#'
lmNULISAseq_predict <- function(data, 
                                sampleInfo,
                                sampleName_var,
                                response_var,
                                modelFormula,
                                reduced_modelFormula = NULL,
                                exclude_targets = NULL,
                                exclude_samples = NULL,
                                target_subset = NULL,
                                sample_subset = NULL,
                                return_model_fits = FALSE) {
  
  ##############################
  # Input Validation
  ##############################
  
  # check if response variable is numeric
  if (!is.numeric(sampleInfo[[response_var]])) {
    stop("Response variable must be numeric for linear regression models")
  }
  
  ##############################
  # Data Subsetting
  ##############################
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
  
  ##############################
  # Initialize Results Storage
  ##############################
  # define vector of targets
  targets <- rownames(data)
  
  # for reduced model if specified
  if (!is.null(reduced_modelFormula)) {
    reducedFits <- vector('list', length(targets))
    names(reducedFits) <- targets
    Fstats_list <- vector('list', length = length(targets))
    names(Fstats_list) <- targets
  }
  
  # collect all possible predictors across models
  all_predictors <- character()
  
  # create empty objects to store results
  modelFits <- vector(mode='list', length=length(targets))
  names(modelFits) <- targets
  stats_list <- vector(mode='list', length=length(targets))
  names(stats_list) <- targets
  
  ##############################
  # Main Analysis 
  ##############################
  
  for (i in seq_along(targets)) {
    target <- targets[i]
    
    tryCatch({
      # prepare data
      target_data <- data.frame(
        sampleName = colnames(data),
        target = unlist(data[target, ])
      )
      
      model_data <- merge(sampleInfo, target_data,
                          by.x = sampleName_var, 
                          by.y = 'sampleName')
      rownames(model_data) <- model_data[[sampleName_var]]
      
      # create formula
      model_formula <- as.formula(paste0(response_var, " ~ target + ", modelFormula))
      model_data <- model_data[complete.cases(model_data[, all.vars(model_formula)]), ]
      
      # fit full model
      model_fit <- lm(model_formula, data = model_data)
      
      # fit reduced model if specified
      if (!is.null(reduced_modelFormula)) {
        reduced_formula <- as.formula(paste0(response_var, " ~ target + ", reduced_modelFormula))
        
        reduced_fit <- lm(reduced_formula, data = model_data)
        
        # Perform anova Ftest
        anova_test <- suppressMessages(anova(reduced_fit, model_fit))
        Fstats_list[[i]] <- c(Fstat=anova_test$F[2], 
                              Df=anova_test$Df[2],
                              Ftest_pval=anova_test$`Pr(>F)`[2])
        
        if (return_model_fits) {
          reducedFits[[i]] <- reduced_fit
        }
      }
      
      # store fit if requested
      if (return_model_fits) modelFits[[i]] <- model_fit
      
      # convert coef table to a tibble/dataframe with metric column having all relevant rownames
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
      cat(format(Sys.time()), "Index:",i," Target:", targets[i],"\n")
      cat(format(Sys.time()), "Error:", conditionMessage(e),"\n")
    })
  }
  
  ##############################
  # Prepare Results Output
  ##############################
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
  
  ##############################
  # Prepare Final Output
  ##############################
  
  output <- list(modelStats = modelStats)
  
  # Add Ftest results if specified
  if (!is.null(reduced_modelFormula)) {
    Fstats <- do.call(rbind, Fstats_list)
    Fstats <- data.frame(target = targets, Fstats, stringsAsFactors = FALSE)
    Fstats$Ftest_pval_FDR <- p.adjust(Fstats$Ftest_pval, method = "BH")
    Fstats$Ftest_pval_bonf <- p.adjust(Fstats$Ftest_pval, method = "bonferroni")
    output$Fstats <- Fstats
  }
  
  # Add model fits if requested
  if (return_model_fits) {
    output$modelFits <- modelFits
    if (!is.null(reduced_modelFormula)) {
      output$reducedFits <- reducedFits
    }
  }
  
  return(output)
}