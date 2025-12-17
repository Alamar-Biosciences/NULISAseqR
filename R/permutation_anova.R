#' Permutation-Based ANOVA with Post-Hoc Pairwise Comparisons
#'
#' Performs a permutation test for a one-way ANOVA, followed by post-hoc pairwise
#' comparisons between group means. The function computes empirical p-values from
#' a shared permutation distribution and adjusts for multiple testing.
#'
#' @param response A numeric vector of response values.
#' @param group A factor vector indicating group membership for each observation.
#' @param B Number of permutations to perform. Default is 10,000.
#' @param seed Optional integer to set random seed for reproducibility.
#' @param plot Logical; if \code{TRUE}, plots the permutation distributions 
#' for the F-statistic and pairwise differences. Default is \code{TRUE}.
#' @param fast_pairwise Logical; if \code{TRUE}, will randomly sample 
#' only one pair mean difference for each permutation. If \code{FALSE} (default), 
#' will calculate every pairwise difference for each permutation.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{observed_F}}{Observed F-statistic from the linear model.}
#'   \item{\code{permutation_p}}{Permutation-based p-value for the global test.}
#'   \item{\code{perm_F_distribution}}{Vector of F-statistics from permuted datasets.}
#'   \item{\code{pairwise_results}}{Data frame with observed pairwise differences and p-values (raw and adjusted using Bonferroni, Holm, and BH).}
#'   \item{\code{perm_pairwise_distribution}}{Vector of all permuted pairwise mean differences pooled across all comparisons.}
#' }
#'
#' @examples
#' set.seed(123)
#' group <- factor(rep(c("A", "B", "C"), each = 10))
#' response <- rnorm(30) + rep(c(0, 0.5, 1), each = 10)
#' result <- permutation_anova(response, group, B = 1000, seed = 123)
#' result$pairwise_results
#'
#' @export

permutation_anova <- function(response, 
                              group, 
                              B = 10000, 
                              seed = NULL, 
                              plot = TRUE,
                              fast_pairwise = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  
  if (length(response) != length(group)) stop("Response and group must be of the same length.")
  if (!is.factor(group)) group <- factor(group)
  df <- data.frame(response = response, group = group)
  
  # function to calculate f stat (faster than lm and anova)
  custom_f_stat <- function(response, group) {
    group <- factor(group)
    group_means <- tapply(response, group, mean)
    group_vars <- tapply(response, group, var)
    group_ns <- table(group)
    grand_mean <- mean(response)
    
    ssb <- sum(group_ns * (group_means - grand_mean)^2)
    ssw <- sum((group_ns - 1) * group_vars)
    
    k <- length(group_means)
    N <- length(response)
    msb <- ssb / (k - 1)
    msw <- ssw / (N - k)
    
    F <- msb / msw
    return(F)
  }
  
  # Observed global F-statistic
  obs_F <- custom_f_stat(response = response, group = group)
  
  # Observed pairwise group mean differences
  group_levels <- levels(group)
  pairwise_obs <- combn(group_levels, 2, simplify = FALSE)
  obs_diffs <- sapply(pairwise_obs, function(pair) {
    m1 <- mean(response[group == pair[1]])
    m2 <- mean(response[group == pair[2]])
    m1 - m2
  })
  names(obs_diffs) <- sapply(pairwise_obs, paste, collapse = " - ")
  
  # Permutations
  perm_F <- numeric(B)
  if (fast_pairwise) {
    perm_diffs <- numeric(B)
  } else {
    perm_diffs <- matrix(NA, nrow = B, ncol = length(pairwise_obs))
    colnames(perm_diffs) <- names(obs_diffs)
  }
  
  for (b in 1:B) {
    perm_group <- sample(group)
    perm_F[b] <- custom_f_stat(response = response, group = perm_group)
    
    if (fast_pairwise){
      # randomly sample only one pairwise comparison per permutation
      pair <- sample(pairwise_obs, 1)[[1]]
      m1 <- mean(response[perm_group == pair[1]])
      m2 <- mean(response[perm_group == pair[2]])
      perm_diffs[b] <- m1 - m2
    } else {
      for (i in seq_along(pairwise_obs)) {
        pair <- pairwise_obs[[i]]
        m1 <- mean(response[perm_group == pair[1]])
        m2 <- mean(response[perm_group == pair[2]])
        perm_diffs[b, i] <- m1 - m2
      }
    }
    
  }
  
  # Global p-value
  p_F <- mean(perm_F >= obs_F)
  
  # Pairwise raw p-values
  p_pairwise_raw <- sapply(seq_along(obs_diffs), function(i) {
    mean(abs(c(perm_diffs)) >= abs(obs_diffs[i]))
  })
  names(p_pairwise_raw) <- names(obs_diffs)
  
  # Multiple testing corrections
  p_adj_bonf <- p.adjust(p_pairwise_raw, method = "bonferroni")
  p_adj_holm <- p.adjust(p_pairwise_raw, method = "holm")
  p_adj_bh   <- p.adjust(p_pairwise_raw, method = "BH")
  
  # Optional plot
  if (plot) {
    par(mfrow=c(1, 3))
    # boxplot
    boxplot(response ~ group, 
            main = 'Target NPQ distribution by plate',
            xlab= 'Plate',
            ylab = 'NPQ')
    
    # Fstat plot
    hist(perm_F, breaks = 50, col = "lightblue",
         main = "Permutation Distribution of F-statistic",
         xlab = "F-statistic")
    abline(v = obs_F, col = "red", lwd = 2)
    legend("topright", legend = paste("Obs F =", round(obs_F, 3)), col = "red", lwd = 2)
    
    # pairwise diff plot
    hist(c(perm_diffs), breaks = 50, col = "lightblue",
         main = "Permutation Distribution of Pairwise Mean Differences",
         xlab = "Pairwise Mean Difference")
    abline(v = obs_diffs, col = "red", lwd = 1)
  }
  
  pairwise_results <- data.frame(comparison=names(p_pairwise_raw),
                                 obs_diffs=obs_diffs,
                                 pval=p_pairwise_raw,
                                 p_bonf=p_adj_bonf,
                                 p_holm=p_adj_holm,
                                 p_FDR=p_adj_bh)
  
  return(list(
    observed_F = obs_F,
    permutation_F_p = p_F,
    perm_F_distribution = perm_F,
    pairwise_results=pairwise_results,
    perm_pairwise_distribution=c(perm_diffs)
  ))
}
