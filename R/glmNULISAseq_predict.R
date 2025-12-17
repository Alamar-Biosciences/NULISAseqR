#' Generalized Linear Model for NULISAseq Data - 
#' targets as predictors
#'
#' Fits generalized linear model to each target in the NULISAseq data set, 
#' using univariate targets as predictors in the model, supporting various response distributions through the family parameter. 
#' Outputs coefficients, odds ratios, z-statistics, unadjusted and adjusted p-values. 
#'
#' @param data A matrix of normalized NULISAseq data
#' with targets in rows, samples in columns. 
#' Row names should be the target names, and column names are the sample names.
#' It is assumed that data has already been transformed
#' using \code{log2(x + 1)} for each NULISAseq normalized count value \code{x}.
#' @param sampleInfo A data frame with sample metadata including the response 
#' variable and covariates. Rows are samples, 
#' columns are sample metadata variables. Generalized linear models will 
#' only be done on the samples in sampleInfo, or a subset of those samples as 
#' specified using arguments \code{exclude_samples} or \code{sample_subset}. 
#' \code{sampleInfo} should have a column for each 
#' variable included in the generalized linear models. String variables will 
#' be automatically treated as factors, and numeric variables will be 
#' treated as numeric.
#' @param sampleName_var The name of the column of sampleInfo that matches
#' the column names of \code{data}. This variable will be used to merge the 
#' target expression data with the sample metadata. 
#' @param response_var The name of the column of sampleInfo specifying the response 
#' variable.
#' @param modelFormula A string that represents the right hand side of the model 
#' formula (everything after the \code{~}) used for the generalized linear model. 
#' The main effect of target expression will be automatically added as a predictor.
#' Any interactions need to be specified in the model formula as "covariate * target".
#' For example when \code{response_var = disease}, \code{family = binomial()}, \code{modelFormula = 
#' "age + sex + plate"} tests for differences in the log-odds of the binary outcome
#' defined by the response variable (e.g., "disease"), adjusted for age, sex, and plate. \code{modelFormula = 
#' "sex * target + age + plate"} includes both main and interaction 
#' effects for sex and target expression. See \code{?glm()}.
#' @param reduced_modelFormula Optional reduced model formula 
#' that contains only a subset of the terms in modelFormula. 
#' The reduced model serves as null model for a likelihood ratio test 
#' (LRT, which is a Chi-square test) using \code{anova()}.
#' This could be useful for testing the overall significance of factor 
#' variables with more than 2 levels. 
#' @param exclude_targets A vector of target names for targets that will be 
#' excluded from the the generalized linear models as predictors. Internal control targets, 
#' for example, should probably always be excluded.
#' @param exclude_samples A vector of sample names for samples that will be excluded
#' from the generalized linear models. External control wells (IPCs, NCs, SC,)
#' should usually be excluded.
#' @param target_subset Overrides exclude_targets. A vector of target names 
#' for targets that will be included in the generalized linear models as predictors.
#' @param sample_subset Overrides exclude_samples. A vector of sample names 
#' for samples that will be included in the generalized linear models.
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
#' 
#'
#' @export
#' 
glmNULISAseq_predict <- function(data, 
                                 sampleInfo,
                                 sampleName_var,
                                 response_var,
                                 modelFormula,
                                 reduced_modelFormula = NULL,
                                 exclude_targets = NULL,
                                 exclude_samples = NULL,
                                 target_subset = NULL,
                                 sample_subset = NULL,
                                 family = binomial(),
                                 return_model_fits = FALSE) {
  
  ##############################
  # Input Validation
  ##############################
  
  # Check if family is valid
  if (!inherits(family, "family")) {
    stop("family must be a valid family object (e.g., binomial(), gaussian())")
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
  # Get data target subset
  if(!is.null(exclude_targets) & is.null(target_subset)){
    data <- data[!(rownames(data) %in% exclude_targets), , drop = FALSE]
  }
  if(!is.null(target_subset)){
    data <- data[target_subset, , drop = FALSE]
  }
  # Get sample subset
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
  if (!is.null(reduced_modelFormula)) {
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
      model_formula <- as.formula(paste0(response_var, " ~ target + ", modelFormula))
      model_data <- model_data[complete.cases(model_data[, all.vars(model_formula)]), ]
      
      # Fit full model
      model_fit <- glm(model_formula, data = model_data, family = family)
      
      # Fit reduced model if specified
      if (!is.null(reduced_modelFormula)) {
        reduced_formula <- as.formula(paste0(response_var, " ~ target + ", reduced_modelFormula))
        reduced_fit <- glm(reduced_formula, data = model_data, family = family)
        
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
  if (!is.null(reduced_modelFormula)) {
    LRTstats <- do.call(rbind, LRTstats_list)
    LRTstats <- data.frame(target = targets, LRTstats, stringsAsFactors = FALSE)
    LRTstats$Chisq_test_pval_FDR <- p.adjust(LRTstats$Chisq_test_pval, method = "BH")
    LRTstats$Chisq_test_pval_bonf <- p.adjust(LRTstats$Chisq_test_pval, method = "bonferroni")
    output$LRTstats <- LRTstats
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


#' Calculate effect size from GLM coefficients
#'
#' Computes the appropriate effect size measure based on the GLM family.
#' For binomial and quasi-binomial families, returns odds ratios (OR).
#' For Poisson and quasi-Poisson families, returns rate ratios (RR).
#' For Gamma and inverse Gaussian families, returns inverse ratios (IR).
#' For other families, returns raw coefficients.
#'
#' @param coef_table A coefficient table from \code{summary(glm_object)$coefficients}
#'   containing at minimum an "Estimate" column with regression coefficients
#' @param family A family object from the GLM specifying the distribution and link function
#'
#' @return A named numeric vector of effect sizes with names matching the
#'   predictor names from the coefficient table
#'
#' @details
#' The transformation applied depends on the family:
#' \itemize{
#'   \item binomial, quasibinomial: exp(estimate) for odds ratios
#'   \item poisson, quasipoisson: exp(estimate) for rate ratios
#'   \item Gamma, inverse.gaussian: 1/estimate for inverse ratios
#'   \item quasi: depends on variance function (mu(1-mu) -> OR, mu -> RR)
#'   \item gaussian and others: raw coefficient (no transformation)
#' }
#'
#' @keywords internal
get_effect_size <- function(coef_table, family) {
  estimate <- coef_table[, "Estimate"]
  
  effect <- switch(family$family,
                   "binomial" = exp(estimate),       # Odds Ratio
                   "poisson" = exp(estimate),       # Rate Ratio
                   "quasipoisson" = exp(estimate),  # Rate Ratio
                   "Gamma" = 1/estimate,            # Inverse (scale parameter)
                   "inverse.gaussian" = 1/estimate, # Inverse
                   "quasibinomial" = exp(estimate), # Odds Ratio
                   estimate)                         # Raw coefficient for others
  
  # Special handling for quasi families
  if (family$family == "quasi") {
    if (!is.null(family$varfun)) {
      if (family$varfun == "mu(1-mu)") effect <- exp(estimate) # quasibinomial
      else if (family$varfun == "mu") effect <- exp(estimate)  # quasipoisson
    }
  }
  
  return(effect)
}

#' Get effect size column suffix for GLM family
#'
#' Returns the appropriate suffix for effect size columns based on the GLM family.
#' This suffix is appended to predictor names when creating output column names
#' in the modelStats data frame.
#'
#' @param family A family object from the GLM specifying the distribution and link function
#'
#' @return A character string suffix:
#' \describe{
#'   \item{_OR}{Odds ratio (binomial, quasibinomial)}
#'   \item{_RR}{Rate ratio (poisson, quasipoisson)}
#'   \item{_IR}{Inverse ratio (Gamma, inverse.gaussian)}
#'   \item{_coef}{Raw coefficient (gaussian and other families)}
#' }
#'
#' @keywords internal
get_effect_suffix <- function(family) {
  suffix <- switch(family$family,
                   "binomial" = "_OR",
                   "poisson" = "_RR",
                   "quasipoisson" = "_RR",
                   "Gamma" = "_IR",
                   "inverse.gaussian" = "_IR",
                   "quasibinomial" = "_OR",
                   "_coef") # Default
  
  # Handle quasi family
  if (family$family == "quasi") {
    if (!is.null(family$varfun)) {
      if (family$varfun == "mu(1-mu)") suffix <- "_OR"
      else if (family$varfun == "mu") suffix <- "_RR"
    }
  }
  
  return(suffix)
}

#' Get test statistic column suffix for GLM family
#'
#' Returns the appropriate suffix for test statistic columns based on the GLM family.
#' GLMs use z-values for most families but t-values for quasi-likelihood families
#' with unknown dispersion. This suffix is appended to predictor names when creating
#' output column names in the modelStats data frame.
#'
#' @param family A family object from the GLM specifying the distribution and link function
#'
#' @return A character string suffix:
#' \describe{
#'   \item{_zval}{z-statistic (binomial, poisson, quasipoisson, Gamma, inverse.gaussian, quasibinomial)}
#'   \item{_tval}{t-statistic (gaussian, quasi with certain variance functions)}
#' }
#'
#' @details
#' Most GLM families use asymptotic z-tests because the dispersion parameter is known
#' or estimated from the family. The quasi-likelihood family may use t-tests when
#' the dispersion is estimated from the data, unless the variance function corresponds
#' to a standard family (binomial or Poisson).
#'
#' @keywords internal
get_stat_suffix <- function(family) {
  if (family$family %in% c("binomial", "poisson", "quasipoisson", 
                           "Gamma", "inverse.gaussian", "quasibinomial")) {
    return("_zval")
  } else if (family$family == "quasi") {
    if (!is.null(family$varfun)) {
      if (family$varfun %in% c("mu(1-mu)", "mu")) return("_zval")
    }
    return("_tval")
  } else {
    return("_tval")
  }
}