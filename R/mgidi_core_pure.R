#' MGIDI Core Algorithm
#'
#' @description
#' `r badge('experimental')`
#'
#' Pure R implementation of the MGIDI (Multi-trait Genotype-Ideotype Distance
#' Index) core algorithm. This implementation provides transparent, educational
#' code without relying on external statistical packages for core computations.
#'
#' The MGIDI index is computed as:
#' \loadmathjax
#' \mjsdeqn{MGIDI_i = \sqrt{\sum\limits_{j = 1}^f(F_{ij} - {F_j})^2}}
#'
#' @name mgidi_core
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @md

# Constant to avoid numerical issues in MGIDI factor contribution calculations
# When a genotype scores exactly the same as the ideotype on a factor, the
# difference is zero, which would cause division by zero when computing
# relative factor contributions (contrib_fac = |diff| / sum(|diff|) * 100).
# This small value preserves the ranking while avoiding computational errors.
# The value 1e-10 is chosen to be negligibly small compared to typical factor
# score differences while still being representable in double precision.
MGIDI_ZERO_TOLERANCE <- 1e-10

#' MGIDI Index Computation (Pure R)
#'
#' @description
#' Computes the Multi-trait Genotype-Ideotype Distance Index using pure R
#' implementation without external package dependencies for core computations.
#'
#' @param data A matrix or data frame with genotypes in rows and traits in columns.
#'   The first column can be genotype names (character/factor) or all columns can
#'   be numeric traits.
#' @param ideotype_directions A named vector or simple vector indicating direction
#'   for each trait. Use +1 or "h" for traits to maximize, -1 or "l" for traits
#'   to minimize. Names should match column names in data.
#' @param n_factors Number of factors to retain. If NULL (default), uses Kaiser
#'   criterion (eigenvalues >= 1).
#' @param mineval Minimum eigenvalue for factor retention. Default is 1.
#' @param weights Optional numeric vector of weights for each trait. Default is
#'   equal weights (all 1s).
#' @param use Method for handling missing values in correlation. Default is
#'   "complete.obs".
#' @param verbose Logical. If TRUE (default), print progress information.
#'
#' @return An object of class `mgidi_pure` with components:
#' \itemize{
#'   \item data: Original trait data
#'   \item rescaled_data: Rescaled trait data (0-100 scale)
#'   \item cormat: Correlation matrix
#'   \item PCA: PCA results (eigenvalues, variance explained)
#'   \item initial_loadings: Initial factor loadings
#'   \item rotated_loadings: Varimax-rotated factor loadings
#'   \item canonical_loadings: Canonical loadings for score computation
#'   \item scores_gen: Factor scores for all genotypes
#'   \item scores_ide: Factor scores for the ideotype
#'   \item MGIDI: Data frame with genotypes and MGIDI values
#'   \item contri_fac: Factor contribution to MGIDI for each genotype
#'   \item communalities: Communalities for each trait
#'   \item KMO: Kaiser-Meyer-Olkin measure
#'   \item MSA: Measure of Sampling Adequacy for each trait
#' }
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' # Create example data
#' set.seed(123)
#' data <- data.frame(
#'   GEN = paste0("G", 1:20),
#'   GY = rnorm(20, 5, 1),    # Grain yield (maximize)
#'   HI = rnorm(20, 0.4, 0.05), # Harvest index (maximize)
#'   PH = rnorm(20, 100, 10)  # Plant height (minimize)
#' )
#'
#' # Define ideotype directions
#' directions <- c(GY = 1, HI = 1, PH = -1)
#'
#' # Run MGIDI
#' result <- mgidi_index(data, directions)
#'
#' # View results
#' result$MGIDI
#' }
mgidi_index <- function(data,
                        ideotype_directions,
                        n_factors = NULL,
                        mineval = 1,
                        weights = NULL,
                        use = "complete.obs",
                        verbose = TRUE) {

  # Handle data input
  data <- as.data.frame(data)

  # Check for genotype column
  is_first_col_char <- !is.numeric(data[[1]])

  if (is_first_col_char) {
    gen_names <- as.character(data[[1]])
    trait_data <- data[, -1, drop = FALSE]
  } else {
    gen_names <- paste0("G", 1:nrow(data))
    trait_data <- data
  }

  # Convert to numeric matrix
  trait_data <- as.matrix(trait_data)
  rownames(trait_data) <- gen_names
  nvar <- ncol(trait_data)
  ngen <- nrow(trait_data)
  var_names <- colnames(trait_data)

  # Validate number of traits
  if (nvar < 2) {
    stop("MGIDI requires at least 2 traits.", call. = FALSE)
  }

  # Process ideotype directions
  directions <- validate_directions(ideotype_directions, nvar, var_names)

  # Set default weights
  if (is.null(weights)) {
    weights <- rep(1, nvar)
  } else {
    if (length(weights) != nvar) {
      stop("Length of weights must equal number of traits.", call. = FALSE)
    }
  }

  # Step 1: Rescale data to 0-100 scale based on directions
  rescaled <- matrix(0, nrow = ngen, ncol = nvar)
  rownames(rescaled) <- gen_names
  colnames(rescaled) <- var_names

  for (j in 1:nvar) {
    trait_vals <- trait_data[, j]
    trait_min <- min(trait_vals, na.rm = TRUE)
    trait_max <- max(trait_vals, na.rm = TRUE)

    if (directions[j] == 1) {
      # Maximize: higher values get higher scores (closer to 100)
      rescaled[, j] <- 100 * (trait_vals - trait_min) / (trait_max - trait_min)
    } else {
      # Minimize: lower values get higher scores (closer to 100)
      rescaled[, j] <- 100 * (trait_max - trait_vals) / (trait_max - trait_min)
    }
  }

  # Step 2: Compute correlation matrix
  cormat <- cor(rescaled, use = use)

  # Step 3: Eigendecomposition
  eig <- eigen(cormat)
  eigenvalues <- eig$values
  eigenvectors <- eig$vectors
  colnames(eigenvectors) <- paste0("PC", 1:nvar)
  rownames(eigenvectors) <- var_names

  # Step 4: Determine number of factors
  if (is.null(n_factors)) {
    n_factors <- sum(eigenvalues >= mineval)
    if (n_factors == 0) {
      n_factors <- 1
      if (verbose) {
        message("No eigenvalues >= ", mineval, ". Using 1 factor.")
      }
    }
  }

  # Step 5: Compute initial loadings
  if (n_factors == 1) {
    initial_loadings <- eigenvectors[, 1, drop = FALSE] * sqrt(eigenvalues[1])
    colnames(initial_loadings) <- "FA1"
  } else {
    sqrt_eigenvalues <- sqrt(eigenvalues[1:n_factors])
    initial_loadings <- eigenvectors[, 1:n_factors] %*% diag(sqrt_eigenvalues)
    colnames(initial_loadings) <- paste0("FA", 1:n_factors)
  }
  rownames(initial_loadings) <- var_names

  # Step 6: Apply Varimax rotation (if more than 1 factor)
  if (n_factors > 1) {
    rotated_result <- varimax(initial_loadings)
    rotated_loadings <- rotated_result$loadings[]
    class(rotated_loadings) <- "matrix"
  } else {
    rotated_loadings <- initial_loadings
  }
  colnames(rotated_loadings) <- paste0("FA", 1:n_factors)
  rownames(rotated_loadings) <- var_names

  # Step 7: Compute KMO and MSA
  R_inv <- solve(cormat)
  k <- nvar
  partial <- R_inv
  for (i in 1:k) {
    for (j in 1:k) {
      if (i != j) {
        partial[i, j] <- -R_inv[i, j] / sqrt(R_inv[i, i] * R_inv[j, j])
      }
    }
  }

  KMO <- sum((cormat[!diag(k)])^2) /
    (sum((cormat[!diag(k)])^2) + sum((partial[!diag(k)])^2))

  MSA <- sapply(1:k, function(i) {
    sum_r2 <- sum(cormat[i, -i]^2)
    sum_q2 <- sum(partial[i, -i]^2)
    sum_r2 / (sum_r2 + sum_q2)
  })
  names(MSA) <- var_names

  # Step 8: Compute communalities
  comm <- rowSums(rotated_loadings^2)
  unique <- 1 - comm

  # Step 9: Compute canonical loadings and factor scores
  # Standardize rescaled data by SD only (center = FALSE)
  sd_vals <- apply(rescaled, 2, sd, na.rm = TRUE)
  Z <- scale(rescaled, center = FALSE, scale = sd_vals)

  # Canonical loadings
  canonical_loadings <- t(t(rotated_loadings) %*% R_inv)
  colnames(canonical_loadings) <- paste0("FA", 1:n_factors)
  rownames(canonical_loadings) <- var_names

  # Factor scores
  scores <- Z %*% canonical_loadings
  rownames(scores) <- gen_names
  colnames(scores) <- paste0("FA", 1:n_factors)

  # Step 10: Construct ideotype
  # Ideotype has maximum rescaled values (all 100) weighted
  ideotype_vals <- rep(100, nvar) * weights
  ideotype_scaled <- ideotype_vals / sd_vals
  ideotype_scores <- matrix(ideotype_scaled %*% canonical_loadings, nrow = 1)
  rownames(ideotype_scores) <- "IDEOTYPE"
  colnames(ideotype_scores) <- paste0("FA", 1:n_factors)

  # Step 11: Compute MGIDI (Euclidean distance from ideotype)
  gen_ide <- sweep(scores, 2, ideotype_scores, "-")
  # Avoid exact zeros to prevent issues in factor contribution calculation
  gen_ide[gen_ide == 0] <- MGIDI_ZERO_TOLERANCE

  mgidi_vals <- sqrt(rowSums(gen_ide^2))
  mgidi_order <- order(mgidi_vals)
  mgidi_sorted <- mgidi_vals[mgidi_order]

  # Create MGIDI data frame
  MGIDI <- data.frame(
    Genotype = names(mgidi_sorted),
    MGIDI = as.numeric(mgidi_sorted),
    stringsAsFactors = FALSE
  )

  # Step 12: Compute factor contributions
  contri_fac <- abs(gen_ide) / rowSums(abs(gen_ide)) * 100
  contri_fac <- data.frame(
    GEN = rownames(contri_fac),
    contri_fac,
    stringsAsFactors = FALSE
  )

  # PCA summary
  pca_summary <- data.frame(
    PC = paste0("PC", 1:nvar),
    Eigenvalues = eigenvalues,
    Variance_pct = eigenvalues / sum(eigenvalues) * 100,
    Cum_variance_pct = cumsum(eigenvalues / sum(eigenvalues) * 100)
  )

  # Factor analysis summary
  fa_summary <- data.frame(
    VAR = var_names,
    rotated_loadings,
    Communality = comm,
    Uniqueness = unique,
    stringsAsFactors = FALSE
  )

  # Print summary if verbose
  if (verbose) {
    cat("\n-------------------------------------------------------------------------------\n")
    cat("Principal Component Analysis\n")
    cat("-------------------------------------------------------------------------------\n")
    print(pca_summary, row.names = FALSE)
    cat("\n-------------------------------------------------------------------------------\n")
    cat("Factor Analysis - Loadings after Varimax rotation\n")
    cat("-------------------------------------------------------------------------------\n")
    print(round(fa_summary[, -1], 4), row.names = FALSE)
    cat("\n-------------------------------------------------------------------------------\n")
    cat("Communality Mean:", round(mean(comm), 4), "\n")
    cat("KMO:", round(KMO, 4), "\n")
    cat("-------------------------------------------------------------------------------\n")
    cat("Top 10 genotypes by MGIDI (lower is better):\n")
    cat("-------------------------------------------------------------------------------\n")
    print(head(MGIDI, 10), row.names = FALSE)
    cat("-------------------------------------------------------------------------------\n")
  }

  # Return results
  result <- list(
    data = trait_data,
    rescaled_data = rescaled,
    cormat = cormat,
    PCA = pca_summary,
    FA = fa_summary,
    initial_loadings = initial_loadings,
    rotated_loadings = rotated_loadings,
    canonical_loadings = canonical_loadings,
    scores_gen = scores,
    scores_ide = ideotype_scores,
    gen_ide = gen_ide,
    MGIDI = MGIDI,
    contri_fac = contri_fac,
    communalities = comm,
    KMO = KMO,
    MSA = MSA,
    directions = directions,
    weights = weights,
    n_factors = n_factors
  )

  class(result) <- c("mgidi_pure", "list")
  return(result)
}


#' Validate Ideotype Directions
#'
#' @description
#' Validates and standardizes the ideotype directions vector.
#'
#' @param directions A vector of directions. Can be numeric (+1/-1) or
#'   character ("h"/"l"). Can be named or unnamed.
#' @param n_traits Number of traits expected
#' @param trait_names Names of traits (for validation if directions are named)
#'
#' @return A numeric vector of directions (+1 or -1) with trait names
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' # Numeric directions
#' validate_directions(c(1, 1, -1), 3, c("GY", "HI", "PH"))
#'
#' # Character directions
#' validate_directions(c("h", "h", "l"), 3, c("GY", "HI", "PH"))
#'
#' # Named directions
#' validate_directions(c(GY = 1, HI = 1, PH = -1), 3, c("GY", "HI", "PH"))
#' }
validate_directions <- function(directions, n_traits, trait_names = NULL) {
  # Handle character input (comma-separated string)
  if (length(directions) == 1 && is.character(directions) && grepl(",", directions)) {
    directions <- unlist(strsplit(directions, "\\s*(,|\\s)\\s*"))
    directions <- tolower(trimws(directions))
  }

  # Convert character directions to numeric
  if (is.character(directions)) {
    directions <- sapply(directions, function(x) {
      x <- tolower(trimws(x))
      if (x %in% c("h", "high", "max", "maximize", "increase")) {
        return(1)
      } else if (x %in% c("l", "low", "min", "minimize", "decrease")) {
        return(-1)
      } else if (!is.na(suppressWarnings(as.numeric(x)))) {
        return(as.numeric(x))
      } else {
        stop("Invalid direction: '", x, "'. Use 'h' or 'l' or numeric values.", call. = FALSE)
      }
    })
  }

  # Ensure numeric
  directions <- as.numeric(directions)

  # Validate length
  if (length(directions) != n_traits) {
    stop("Length of ideotype_directions (", length(directions),
         ") must equal number of traits (", n_traits, ").", call. = FALSE)
  }

  # Standardize to +1 or -1
  directions <- sign(directions)
  directions[directions == 0] <- 1  # Default to maximize if 0

  # Add names
  if (!is.null(trait_names)) {
    names(directions) <- trait_names
  }

  return(directions)
}


#' Ideotype Construction
#'
#' @description
#' Constructs the ideal genotype based on factor scores and trait directions.
#'
#' @param factor_scores Matrix of factor scores (genotypes x factors)
#' @param loadings Matrix of factor loadings (traits x factors)
#' @param directions Vector of trait directions (+1 for maximize, -1 for minimize)
#'
#' @return A numeric vector of ideotype factor scores
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' # Example with simulated data
#' scores <- matrix(rnorm(20), nrow = 10)
#' loadings <- matrix(c(0.8, 0.1, 0.1, 0.9), nrow = 2)
#' directions <- c(1, -1)
#' ideotype_construction(scores, loadings, directions)
#' }
ideotype_construction <- function(factor_scores, loadings, directions) {
  n_factors <- ncol(loadings)
  ideotype <- numeric(n_factors)

  for (j in 1:n_factors) {
    # Compute weighted sum of loadings by directions
    weighted_sum <- sum(loadings[, j] * directions)

    # If positive, ideotype has maximum; if negative, minimum
    if (weighted_sum > 0) {
      ideotype[j] <- max(factor_scores[, j])
    } else {
      ideotype[j] <- min(factor_scores[, j])
    }
  }

  names(ideotype) <- colnames(loadings)
  return(ideotype)
}


#' MGIDI Distance Calculation
#'
#' @description
#' Calculates the Euclidean distance from each genotype to the ideotype.
#'
#' @param factor_scores Matrix of factor scores (genotypes x factors)
#' @param ideotype Vector of ideotype factor scores
#'
#' @return A named vector of MGIDI distances
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' scores <- matrix(rnorm(20), nrow = 10)
#' rownames(scores) <- paste0("G", 1:10)
#' ideotype <- c(2, 2)
#' mgidi_distance(scores, ideotype)
#' }
mgidi_distance <- function(factor_scores, ideotype) {
  # Compute distance from each genotype to ideotype
  diff_matrix <- sweep(factor_scores, 2, ideotype, "-")
  distances <- sqrt(rowSums(diff_matrix^2))
  return(distances)
}


#' Selection Differential Calculation
#'
#' @description
#' Computes the selection differential (genetic gains) for selected genotypes.
#'
#' @param original_data Matrix of original trait values (genotypes x traits)
#' @param selected_genotypes Character vector of selected genotype names
#' @param trait_directions Vector of trait directions (+1/-1)
#'
#' @return A data frame with selection differential information
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' # Create example data
#' data <- matrix(rnorm(60), nrow = 20, ncol = 3)
#' rownames(data) <- paste0("G", 1:20)
#' colnames(data) <- c("GY", "HI", "PH")
#' selected <- c("G1", "G5", "G10")
#' directions <- c(1, 1, -1)
#' selection_differential(data, selected, directions)
#' }
selection_differential <- function(original_data, selected_genotypes, trait_directions) {
  original_data <- as.matrix(original_data)

  # Population means
  pop_means <- colMeans(original_data, na.rm = TRUE)

  # Selected means
  sel_data <- original_data[selected_genotypes, , drop = FALSE]
  sel_means <- colMeans(sel_data, na.rm = TRUE)

  # Selection differential
  SD <- sel_means - pop_means
  SD_pct <- (SD / abs(pop_means)) * 100

  # Determine if selection is in desired direction
  desired <- trait_directions
  goal_achieved <- ifelse(
    (desired == 1 & SD > 0) | (desired == -1 & SD < 0),
    "Yes", "No"
  )

  result <- data.frame(
    Trait = colnames(original_data),
    Direction = ifelse(desired == 1, "Increase", "Decrease"),
    Xo = pop_means,
    Xs = sel_means,
    SD = SD,
    SD_pct = SD_pct,
    Goal = goal_achieved,
    stringsAsFactors = FALSE
  )

  return(result)
}
