#' Linear regression model for differential expression test
#'
#' Fits linear regression model to each target in the data set. 
#' Outputs coefficients, unadjusted and adjusted p-values.
#'
#' @param data A matrix with targets in rows, samples in columns. 
#' Row names should be the target names, and column names are the sample names.
#' For NULISAseq data \code{x}, it is assumed that data has already been transformed
#' using \code{log2(x + 0.01)}.
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
#' @param modelFormula a string that represents the right hand side of the model 
#' formula (everything after the \code{~}) used for the linear model. For example \code{modelFormula = 
#' "disease + age + sex + plate"} test for differences in target expression 
#' by disease group, adjusted for age, sex, and plate. \code{modelFormula = 
#' "disease * age + sex + plate"} includes both main and interaction 
#' effects for disease and age. See \code{?lm()}.
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
#' detailed analyses and plotting. However also requires using more memory.
#'
#' @return A list including the following:
#' \item{stats}{A data frame with rows corresponding to targets and columns 
#' corresponding to estimated model coefficients, unadjusted p-values, 
#' Bonferroni adjusted p-values, and Benjamini-Hochberg false discovery rate
#' adjusted p-values (see \code{?p.adjust()})}
#' \item{modelFits}{A list of length equal to number of targets containing
#' the model fit output from \code{lm()}}
#'
#' 
#'
#' @export
#'
lmNULISAseq <- function(data, 
                        sampleInfo,
                        sampleName_var,
                        modelFormula,
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
    data <- data[,!(rownames(data) %in% exclude_samples)]
  } 
  if(!is.null(sample_subset)){
    data <- data[,sample_subset]
  }
  # define vector of targets
  targets <- rownames(data)
  # create empty objects to store results
  modelFits <- vector(mode=list, length=length(targets))
  # loop over targets and fit model
  for(i in 1:length(targets)){
    target <- targets[i]
    model_formula <- as.formula(paste0(target, '~', modelFormula))
    model_fit <- lm(model_formula, data=data)
  }
  
  return(cv_matrix)
}
