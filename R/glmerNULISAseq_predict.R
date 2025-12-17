#' Generalized Linear Mixed Effects Model for NULISAseq Data - 
#' targets as predictors
#'
#' Fits generalized linear mixed effects model to each target in the NULISAseq data set, 
#' using univariate targets as predictors in the model, supporting various response distributions through the family parameter. 
#' Outputs coefficients, odds ratios, z-statistics, unadjusted and adjusted p-values. Uses \code{lme4} package.
#'
#' @param data A matrix of normalized NULISAseq data
#' with targets in rows, samples in columns. 
#' Row names should be the target names, and column names are the sample names.
#' It is assumed that data has already been transformed
#' using \code{log2(x + 1)} for each NULISAseq normalized count value \code{x}.
#' @param sampleInfo A data frame with sample metadata including the response 
#' variable and covariates. Rows are samples, 
#' columns are sample metadata variables. Generalized linear mixed effects model will 
#' only be done on the samples in sampleInfo, or a subset of those samples as 
#' specified using arguments \code{exclude_samples} or \code{sample_subset}. 
#' \code{sampleInfo} should have a column for each 
#' variable included in the generalized linear mixed effect models. String variables will 
#' be automatically treated as factors, and numeric variables will be 
#' treated as numeric.
#' @param sampleName_var The name of the column of sampleInfo that matches
#' the column names of \code{data}. This variable will be used to merge the 
#' target expression data with the sample metadata. 
#' @param response_var The name of the column of sampleInfo specifying the response 
#' variable.
#' @param modelFormula_fixed A string that represents the fixed effects part of the model 
#' formula used for the generalized linear mixed effects model. 
#' The main effect of target expression will be automatically added as a predictor.
#' Any interactions need to be specified in the model formula as "covariate * target".
#' For example, when \code{response_var = disease}, \code{family = binomial()}, \code{modelFormula_fixed = 
#' "age + sex + plate"} tests for differences in the log-odds of the binary outcome
#' (disease status) as a function of age, sex, and plate. \code{modelFormula_fixed = 
#' "sex * target + age + plate"} includes both main and interaction 
#' effects for sex and target expression. See \code{?glmer()}.
#' @param modelFormula_random A string that represents the random effects part of the model 
#' formula on the used for the generalized linear mixed effects model. 
#' For example \code{modelFormula_random = "(1|participant_ID)"} 
#' creates a subject specific random intercept, where the variable 
#' \code{participant_ID} (a column in \code{sampleInfo} data frame) denotes 
#' repeated measures on the same subject. For subject-specific random intercept and 
#' slopes (not recommended when time is categorical), 
#' use \code{modelFormula_random = "(1 + time|participant_ID)"}. For random 
#' subject nested within plate (which may be useful when analyzing
#' a large number of plates together), use \code{modelFormula_random = 
#' "(1|plate_ID:participant_ID)"}.
#' See \code{?glmer()}.
#' @param reduced_modelFormula_fixed Optional reduced model formula 
#' for generalized linear mixed effects that contains only a subset of the terms in modelFormula. 
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
#' excluded from the generalized linear mixed effects models as predictors. Internal control targets, 
#' for example, should probably always be excluded.
#' @param exclude_samples A vector of sample names for samples that will be excluded
#' from the generalized linear mixed effects models. External control wells (IPCs, NCs, SC,) should usually be excluded.
#' @param target_subset Overrides exclude_targets. A vector of target names 
#' for targets that will be included in the generalized linear mixed effects models as predictors.
#' @param sample_subset Overrides exclude_samples. A vector of sample names 
#' for samples that will be included in the generalized linear mixed effects models.
#' @param family A family object for the \code{glm()} function specifying the distribution and link function.
#'     \describe{
#'     \item{\code{binomial(link = "logit")}}{For binary response data (default). Reports odds ratios (OR).}
#'     \item{\code{gaussian(link = "identity")}}{For continuous normally-distributed data. Reports coefficients.}
#'     \item{\code{poisson(link = "log")}}{For count data. Reports rate ratios (RR).}
#'     \item{\code{Gamma(link = "inverse")}}{For positive continuous data with constant coefficient of variation. Reports inverse ratios (IR).}
#'     \item{\code{inverse.gaussian(link = "1/mu^2")}}{For positive continuous data with variance increasing with mean^3. Reports inverse ratios (IR).}
#'     \item{\code{quasibinomial(link = "logit")}}{For overdispersed binary data. Reports odds ratios (OR).}
#'     \item{\code{quasipoisson(link = "log")}}{For overdispersed count data. Reports rate ratios (RR).}
#'     \item{\code{quasi(link = "identity", variance = "constant")}}{For custom quasi-likelihood models.}
#'   }
#'   The function automatically selects appropriate test statistics (z/t-values) and effect size measures based on the family.
#'   (see \code{?family()})
#' @param return_model_fits Logical \code{TRUE} or \code{FALSE} (default).
#' Should a list of the model fits be returned? Might be useful for more 
#' detailed analyses and plotting. However, also requires using more memory.
#' @param control A list of control parameters for \code{glmer} model fitting, 
#' created by \code{lme4::glmerControl()}. Defaults to \code{glmerControl(optimizer = "bobyqa")}
#' which often helps with convergence issues. Other useful optimizers include "nloptwrap".
#' Additional control parameters can help with convergence:
#' \itemize{
#'   \item \code{optCtrl = list(maxfun = 100000)} - Increase maximum number of function evaluations
#'   \item \code{calc.derivs = FALSE} - Skip derivative calculations if having convergence issues
#'   \item \code{check.conv.grad = FALSE} - Skip gradient convergence checks if needed
#' }
#' See \code{?lme4::glmerControl} for all available options.
#'
#' @return A list including the following:
#' \item{modelStats}{A data frame with rows corresponding to targets and columns 
#' corresponding to estimated model coefficients, effect sizes, standard errors, test statistics, unadjusted p-values, 
#' Bonferroni adjusted p-values, and Benjamini-Hochberg false discovery rate
#' adjusted p-values (see \code{?p.adjust()}).}
#' \item{LRTstats}{A data frame with rows corresponding to targets and columns with likelihood ratio 
#' test statistics including Chi-square statistic, degrees of freedom, unadjusted p-values, 
#' Bonferroni adjusted p-values, and Benjamini-Hochberg false discovery rate
#' adjusted p-values. (Only included when \code{reduced_modelFormula} is specified.)}
#' \item{modelFits}{A list of length equal to number of targets containing
#' the model fit output from \code{glm()}. Only returned when 
#' \code{return_model_fits=TRUE}.}
#' \item{reducedFits}{A list of length equal to number of targets containing
#' the reduced model fit output from \code{anova()}. Only returned when 
#' \code{return_model_fits=TRUE} and \code{reduced_modelFormula} is specified.}
#'
#' @importFrom lme4 glmer glmerControl
#'
#' @export
#' 
glmerNULISAseq_predict <- function(data, 
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
                                   family = binomial(),
                                   return_model_fits = FALSE,
                                   control = lme4::glmerControl(optimizer = "bobyqa")) {
  
  ##############################
  # Input Validation
  ##############################
  
  # Check if family is valid
  if (!inherits(family, "family")) {
    stop("family must be a valid family object (e.g., binomial(), poisson())")
  }
  
  # Validate response variable based on family
  if (family$family == "binomial" && length(unique(sampleInfo[[response_var]])) != 2) {
    stop("For binomial family, response variable must have exactly 2 levels")
  }
  if (family$family %in% c("poisson", "quasipoisson") && any(sampleInfo[[response_var]] < 0)) {
    stop("For Poisson families, response variable must be non-negative")
  }
  if (family$family %in% c("Gamma", "inverse.gaussian") && any(sampleInfo[[response_var]] <= 0)) {
    stop("For Gamma and inverse.gaussian families, response variable must be positive")
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
  
  # Define vector of targets
  targets <- rownames(data)
  
  ##############################
  # Initialize Results Storage
  ##############################
  
  modelFits <- vector('list', length(targets))
  names(modelFits) <- targets
  stats_list <- vector(mode = 'list', length = length(targets))
  names(stats_list) <- targets
  
  # For reduced model if specified
  if (!is.null(reduced_modelFormula_fixed)) {
    reducedFits <- vector('list', length(targets))
    names(reducedFits) <- targets
    LRTstats_list <- vector('list', length = length(targets))
    names(LRTstats_list) <- targets
  }
  
  ##############################
  # Main Analysis 
  ##############################
  
  for (i in seq_along(targets)) {
    target <- targets[i]
    
    tryCatch({
      # Prepare data
      target_data <- data.frame(
        sampleName = colnames(data),
        target = unlist(data[target, ])
      )
      
      model_data <- merge(sampleInfo, target_data,
                          by.x = sampleName_var, 
                          by.y = 'sampleName')
      rownames(model_data) <- model_data[[sampleName_var]]
      
      # Create formula
      model_formula <- as.formula(paste0(response_var, " ~ target + ", 
                                         modelFormula_fixed, " + ", modelFormula_random))
      model_data <- model_data[complete.cases(model_data[, all.vars(model_formula)]), ]
      
      # Fit full model
      model_fit <- lme4::glmer(model_formula, data = model_data, family = family, control = control)
      
      # Fit reduced model if specified
      if (!is.null(reduced_modelFormula_fixed)) {
        reduced_formula <- as.formula(paste0(
          response_var, " ~ target + ", 
          reduced_modelFormula_fixed, " + ",
          ifelse(!is.null(reduced_modelFormula_random), 
                 reduced_modelFormula_random, 
                 modelFormula_random)
        ))
        
        reduced_fit <- lme4::glmer(reduced_formula, data = model_data, family = family, control = control)
        
        # Perform likelihood ratio test
        lrt <- suppressMessages(anova(reduced_fit, model_fit, test = "LRT"))
        LRTstats_list[[i]] <- c(Chisq_stat = lrt$Chisq[2],
                                Df = lrt$Df[2],
                                Chisq_test_pval = lrt$`Pr(>Chisq)`[2])
        
        if (return_model_fits) {
          reducedFits[[i]] <- reduced_fit
        }
      }
      
      # Store fit if requested
      if (return_model_fits) modelFits[[i]] <- model_fit
      
      # Extract coefficients
      fit_summary <- summary(model_fit)
      coef_table <- fit_summary$coefficients
      
      # Convert to named vectors
      coefs <- setNames(coef_table[, "Estimate"], rownames(coef_table))
      ses <- setNames(coef_table[, "Std. Error"], rownames(coef_table))
      
      # Get appropriate test statistic
      stat_col <- if (grepl("z value", colnames(coef_table)[3], fixed = TRUE)) "z value" else "t value"
      stats <- setNames(coef_table[, stat_col], rownames(coef_table))
      
      pvals <- setNames(coef_table[, ncol(coef_table)], rownames(coef_table))
      
      # Calculate appropriate effect size
      effects <- get_effect_size(coef_table, family)
      
      # Clean names (replace (Intercept) with intercept)
      clean_names <- function(x) {
        names(x) <- ifelse(names(x) == "(Intercept)", "intercept", names(x))
        x
      }
      
      coefs <- clean_names(coefs)
      ses <- clean_names(ses)
      stats <- clean_names(stats)
      pvals <- clean_names(pvals)
      effects <- clean_names(effects)
      
      # Store results in stats_list
      stats_list[[i]] <- list(
        coefs = coefs,
        effects = effects,
        ses = ses,
        stats = stats,
        pvals = pvals
      )
      
    }, error = function(e) {
      message(format(Sys.time()), " Index: ", i, " Target: ", targets[i])
      message(format(Sys.time()), " Error: ", conditionMessage(e))
    })
  }
  
  # Check if all models failed
  if (all(sapply(stats_list, is.null))) {
    stop("All targets failed model fitting. Common causes:\n",
         "  - A covariate has only one level after removing samples with missing data\n",
         "  - Random effects grouping variable has insufficient levels\n",
         "  - Insufficient samples remain after filtering.\n",
         "Check the error messages above for details.")
  }
  
  ##############################
  # Prepare Results Output
  ##############################
  
  # Get all predictors across all models
  all_predictors <- unique(unlist(lapply(stats_list, function(x) names(x$coefs))))
  
  # Use safe_extract_matrix to create output matrices
  coef_mat <- safe_extract_matrix(stats_list, "coefs", all_predictors)
  effect_mat <- safe_extract_matrix(stats_list, "effects", all_predictors)
  se_mat <- safe_extract_matrix(stats_list, "ses", all_predictors)
  stat_mat <- safe_extract_matrix(stats_list, "stats", all_predictors)
  p_mat <- safe_extract_matrix(stats_list, "pvals", all_predictors)
  
  # Calculate adjusted p-values
  p_fdr_mat <- apply(p_mat, 2, p.adjust, method = "BH")
  p_bonf_mat <- apply(p_mat, 2, p.adjust, method = "bonferroni")
  
  # Create appropriate column names
  effect_suffix <- get_effect_suffix(family)
  stat_suffix <- get_stat_suffix(family)
  
  colnames(coef_mat) <- paste0(colnames(coef_mat), "_coef")
  colnames(effect_mat) <- paste0(colnames(effect_mat), effect_suffix)
  colnames(se_mat) <- paste0(colnames(se_mat), "_se")
  colnames(stat_mat) <- paste0(colnames(stat_mat), stat_suffix)
  colnames(p_mat) <- paste0(colnames(p_mat), "_pval_unadj")
  colnames(p_fdr_mat) <- paste0(colnames(p_fdr_mat), "_pval_FDR")
  colnames(p_bonf_mat) <- paste0(colnames(p_bonf_mat), "_pval_bonf")
  
  # Combine all results in logical order
  column_order <- c(rbind(
    colnames(coef_mat), 
    colnames(effect_mat),
    colnames(se_mat),
    colnames(stat_mat), 
    colnames(p_mat),
    colnames(p_fdr_mat),
    colnames(p_bonf_mat)
  ))
  
  modelStats <- cbind(coef_mat, effect_mat, se_mat, stat_mat, p_mat, p_fdr_mat, p_bonf_mat)[, column_order]
  modelStats <- data.frame(target = rownames(modelStats), modelStats, stringsAsFactors = FALSE)
  
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
