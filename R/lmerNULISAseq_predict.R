#' Linear Mixed Effect Model for NULISAseq Data - 
#' targets as predictor
#'
#' Fits linear mixed effects model to each target in the NULISAseq data set,
#' using univariate targets as predictors in the model
#' Outcome ~ univariate target predictor + ...
#' Outputs coefficients, t-statistics, unadjusted and adjusted p-values.
#' This approach tests whether each 
#' target's expression is associated with the outcome variable, while 
#' adjusting for any specified fixed effect covariates and accounting 
#' for random effects.
#' 
#' Uses \code{lme4} and \code{lmerTest} packages.
#'
#' @param data A matrix of normalized NULISAseq data
#' with targets used as predictors in rows, samples in columns. 
#' Row names should be the target names, and column names are the sample names.
#' It is assumed that data has already been transformed
#' using \code{log2(x + 1)} for each NULISAseq normalized count value \code{x}.
#' @param sampleInfo A data frame with sample metadata including the response 
#' variable and covariates. Rows are samples, columns are sample metadata variables. 
#' Linear mixed effect models will only be done on the samples in sampleInfo, or a subset of those samples as 
#' specified using arguments \code{exclude_samples} or \code{sample_subset}.
#' @param sampleName_var The name of the column of sampleInfo that matches
#' the column names of \code{data}. This variable will be used to merge the 
#' target expression data with the sample metadata. 
#' @param response_var The name of the column of sampleInfo specifying the continuous numeric response 
#' variable.
#' @param modelFormula_fixed A string that represents the fixed effects part of the model 
#' formula used for the linear mixed effects model. The main effect of target expression 
#' will be automatically added as a predictor. Any interactions need to be specified 
#' in the model formula as "covariate * target". For example, 
#' \code{"disease + age + sex + plate"} tests for variations of outcome explained by the 
#' target predictor, adjusted for disease group, age, sex, and plate. \code{modelFormula = 
#' "disease * target + age + sex + plate"} includes both main and interaction 
#' effects for disease and target expression. See \code{?lmer()}.
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
#' excluded from the linear mixed effect models as predictors. Internal control targets, 
#' for example, should probably always be excluded.
#' @param exclude_samples A vector of sample names for samples that will be excluded
#' from the linear mixed effect models. External control wells (IPCs, NCs, SC,) should usually be excluded.
#' @param target_subset Overrides exclude_targets. A vector of target names 
#' for targets that will be included in the linear mixed effect models as predictors.
#' @param sample_subset Overrides exclude_samples. A vector of sample names 
#' for samples that will be included in the linear mixed effect models.
#' @param return_model_fits Logical \code{TRUE} or \code{FALSE} (default).
#' Should a list of the model fits be returned? Might be useful for more 
#' detailed analyses and plotting. However, also requires using more memory.
#' @param control A list of control parameters for \code{lmer} model fitting, 
#' created by \code{lme4::lmerControl()}. Defaults to \code{lmerControl(optimizer = "bobyqa")}
#' which often helps with convergence issues. Other useful optimizers include "nloptwrap".
#' Additional control parameters can help with convergence:
#' \itemize{
#'   \item \code{optCtrl = list(maxfun = 100000)} - Increase maximum number of function evaluations
#'   \item \code{calc.derivs = FALSE} - Skip derivative calculations if having convergence issues
#'   \item \code{check.conv.grad = FALSE} - Skip gradient convergence checks if needed
#' }
#' See \code{?lme4::lmerControl} for all available options.
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
#' @importFrom magrittr %>%
#'
#' @export
#'
lmerNULISAseq_predict <- function(data, 
                                  sampleInfo,
                                  sampleName_var,
                                  response_var,
                                  modelFormula_fixed,
                                  modelFormula_random,
                                  reduced_modelFormula_fixed = NULL,
                                  reduced_modelFormula_random = NULL,
                                  exclude_targets = NULL,
                                  exclude_samples = NULL,
                                  target_subset = NULL,
                                  sample_subset = NULL,
                                  return_model_fits = FALSE,
                                  control = lme4::lmerControl(optimizer = "bobyqa")) {
  
  ##############################
  # Input Validation
  ##############################
  
  # check if response variable is numeric
  if (!is.numeric(sampleInfo[[response_var]])) {
    stop("Response variable must be numeric for linear mixed-effects models")
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
  if (!is.null(reduced_modelFormula_fixed)) {
    reducedFits <- vector('list', length(targets))
    names(reducedFits) <- targets
    LRTstats_list <- vector('list', length = length(targets))
    names(LRTstats_list) <- targets
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
      model_formula <- as.formula(paste0(response_var, " ~ target + ", 
                                         modelFormula_fixed, " + ", modelFormula_random))
      model_data <- model_data[complete.cases(model_data[, all.vars(model_formula)]), ]
      
      # fit full model
      model_fit <- lmerTest::lmer(model_formula, data = model_data, control = control)
      
      # fit reduced model if specified
      if (!is.null(reduced_modelFormula_fixed)) {
        reduced_formula <- as.formula(paste0(
          response_var, " ~ target + ", 
          reduced_modelFormula_fixed, " + ",
          ifelse(!is.null(reduced_modelFormula_random), 
                 reduced_modelFormula_random, 
                 modelFormula_random)
        ))
        
        reduced_fit <- lmerTest::lmer(reduced_formula, data = model_data, control = control)
        
        # Perform likelihood ratio test
        lrt <- suppressMessages(anova(reduced_fit, model_fit, test = "LRT"))
        LRTstats_list[[i]] <- c(Chisq_stat = lrt$Chisq[2],
                                Df = lrt$Df[2],
                                Chisq_test_pval = lrt$`Pr(>Chisq)`[2])
        
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
  
  # Add LRT results if specified
  if (!is.null(reduced_modelFormula_fixed)) {
    LRTstats <- do.call(rbind, LRTstats_list)
    LRTstats <- data.frame(target = targets, LRTstats, stringsAsFactors = FALSE)
    LRTstats$Chisq_test_pval_FDR <- p.adjust(LRTstats$Chisq_test_pval, method = "BH")
    LRTstats$Chisq_test_pval_bonf <- p.adjust(LRTstats$Chisq_test_pval, method = "bonferroni")
    output$LRTstats <- LRTstats
  }
  
  # Add model fits if requested
  if (return_model_fits) {
    output$modelFits <- modelFits
    if (!is.null(reduced_modelFormula_fixed)) {
      output$reducedFits <- reducedFits
    }
  }
  
  return(output)
}
