#' Factor Analysis from Scratch
#'
#' @description
#' `r badge('experimental')`
#'
#' A collection of pure R functions for factor analysis without relying on
#' external statistical packages. These functions provide transparent
#' implementations of factor analysis algorithms used in MGIDI calculations.
#'
#' @name factor_analysis_pure
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @md

#' Principal Components Analysis on Correlation Matrix
#'
#' @description
#' Performs principal component analysis on a correlation matrix using
#' eigendecomposition (pure R implementation).
#'
#' @param X A numeric matrix or data frame where rows are observations
#'   and columns are variables.
#' @param use Method for handling missing values in correlation computation.
#'   Default is "complete.obs".
#'
#' @return A list with components:
#' \itemize{
#'   \item eigenvalues: A vector of eigenvalues
#'   \item eigenvectors: A matrix of eigenvectors (loadings)
#'   \item variance_explained: Proportion of variance explained by each component
#'   \item cumulative_variance: Cumulative proportion of variance explained
#'   \item correlation_matrix: The correlation matrix used
#' }
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' X <- matrix(rnorm(50), nrow = 10)
#' colnames(X) <- paste0("V", 1:5)
#' pca_result <- principal_components(X)
#' pca_result$eigenvalues
#' pca_result$variance_explained
#' }
principal_components <- function(X, use = "complete.obs") {
  # Convert to matrix if needed
  X <- as.matrix(X)

  # Compute correlation matrix
  R <- cor(X, use = use)

  # Eigendecomposition
  eig <- eigen(R)

  # Extract eigenvalues and eigenvectors
  eigenvalues <- eig$values
  eigenvectors <- eig$vectors

  # Compute variance explained
  total_var <- sum(eigenvalues)
  var_explained <- eigenvalues / total_var
  cum_var <- cumsum(var_explained)

  # Set names
  p <- ncol(X)
  colnames(eigenvectors) <- paste0("PC", 1:p)
  if (!is.null(colnames(X))) {
    rownames(eigenvectors) <- colnames(X)
  }
  names(eigenvalues) <- paste0("PC", 1:p)
  names(var_explained) <- paste0("PC", 1:p)
  names(cum_var) <- paste0("PC", 1:p)

  return(list(
    eigenvalues = eigenvalues,
    eigenvectors = eigenvectors,
    variance_explained = var_explained,
    cumulative_variance = cum_var,
    correlation_matrix = R
  ))
}


#' Varimax Rotation for Factor Loadings
#'
#' @description
#' Applies Varimax rotation to factor loadings to achieve a simpler structure
#' (pure R implementation).
#'
#' @param loadings A matrix of factor loadings (variables x factors)
#' @param max_iter Maximum number of iterations. Default is 100.
#' @param tolerance Convergence tolerance. Default is 1e-5.
#' @param normalize Logical. If TRUE, use Kaiser normalization. Default is TRUE.
#'
#' @return A list with components:
#' \itemize{
#'   \item loadings: The rotated loading matrix
#'   \item rotation_matrix: The rotation matrix applied
#' }
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' # Create sample loadings
#' loadings <- matrix(c(0.8, 0.1, 0.7, 0.2, 0.1, 0.9, 0.2, 0.8), nrow = 4, byrow = TRUE)
#' rotated <- varimax_rotation(loadings)
#' rotated$loadings
#' }
varimax_rotation <- function(loadings, max_iter = 100, tolerance = 1e-5,
                             normalize = TRUE) {
  # Convert to matrix if needed
  loadings <- as.matrix(loadings)

  p <- nrow(loadings)  # Number of variables
  k <- ncol(loadings)  # Number of factors

  # If only one factor, no rotation needed
  if (k == 1) {
    return(list(loadings = loadings, rotation_matrix = matrix(1)))
  }

  # Kaiser normalization
  if (normalize) {
    h2 <- rowSums(loadings^2)
    h <- sqrt(h2)
    h[h == 0] <- 1  # Avoid division by zero
    loadings_norm <- loadings / h
  } else {
    loadings_norm <- loadings
  }

  # Initialize rotation matrix as identity
  rotation <- diag(k)

  # Iterative rotation
  for (iter in 1:max_iter) {
    max_change <- 0

    # Rotate pairs of factors
    for (i in 1:(k - 1)) {
      for (j in (i + 1):k) {
        # Extract columns
        xi <- loadings_norm[, i]
        xj <- loadings_norm[, j]

        # Compute rotation angle (simplified varimax criterion)
        u <- xi^2 - xj^2
        v <- 2 * xi * xj

        A <- sum(u)
        B <- sum(v)
        C <- sum(u^2 - v^2)
        D <- 2 * sum(u * v)

        # Optimal angle
        num <- D - 2 * A * B / p
        den <- C - (A^2 - B^2) / p

        if (abs(den) < tolerance) {
          phi <- 0
        } else {
          phi <- 0.25 * atan2(num, den)
        }

        # Check for convergence
        if (abs(phi) > max_change) {
          max_change <- abs(phi)
        }

        # Apply rotation
        cos_phi <- cos(phi)
        sin_phi <- sin(phi)

        new_xi <- cos_phi * xi + sin_phi * xj
        new_xj <- -sin_phi * xi + cos_phi * xj

        loadings_norm[, i] <- new_xi
        loadings_norm[, j] <- new_xj

        # Update rotation matrix
        rot_ij <- diag(k)
        rot_ij[i, i] <- cos_phi
        rot_ij[j, j] <- cos_phi
        rot_ij[i, j] <- sin_phi
        rot_ij[j, i] <- -sin_phi

        rotation <- rotation %*% rot_ij
      }
    }

    # Check convergence
    if (max_change < tolerance) {
      break
    }
  }

  # Undo Kaiser normalization
  if (normalize) {
    rotated_loadings <- loadings_norm * h
  } else {
    rotated_loadings <- loadings_norm
  }

  # Preserve names
  if (!is.null(rownames(loadings))) {
    rownames(rotated_loadings) <- rownames(loadings)
  }
  if (!is.null(colnames(loadings))) {
    colnames(rotated_loadings) <- colnames(loadings)
  }

  return(list(loadings = rotated_loadings, rotation_matrix = rotation))
}


#' Compute Factor Scores
#'
#' @description
#' Computes factor scores for observations using the regression method
#' (pure R implementation).
#'
#' @param X A numeric matrix or data frame of observations (n x p)
#' @param loadings A matrix of factor loadings (p x k)
#' @param correlation_matrix Optional correlation matrix. If NULL, computed from X.
#' @param standardize Logical. If TRUE, standardize X before computing scores.
#'   Default is TRUE.
#'
#' @return A matrix of factor scores (n x k)
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' X <- matrix(rnorm(50), nrow = 10)
#' colnames(X) <- paste0("V", 1:5)
#' pca <- principal_components(X)
#' # Use first 2 components as loadings
#' loadings <- pca$eigenvectors[, 1:2] * sqrt(pca$eigenvalues[1:2])
#' scores <- factor_scores(X, loadings)
#' scores
#' }
factor_scores <- function(X, loadings, correlation_matrix = NULL,
                          standardize = TRUE) {
  # Convert to matrix if needed
  X <- as.matrix(X)
  loadings <- as.matrix(loadings)

  # Standardize X if requested
  if (standardize) {
    Z <- scale(X, center = TRUE, scale = TRUE)
  } else {
    Z <- X
  }

  # Compute correlation matrix if not provided
  if (is.null(correlation_matrix)) {
    R <- cor(X, use = "complete.obs")
  } else {
    R <- correlation_matrix
  }

  # Compute canonical loadings using regression method
  # Scores = Z * R^-1 * L
  R_inv <- solve(R)
  canonical_loadings <- R_inv %*% loadings

  # Compute scores
  scores <- Z %*% canonical_loadings

  # Set column names
  k <- ncol(loadings)
  if (is.null(colnames(loadings))) {
    colnames(scores) <- paste0("FA", 1:k)
  } else {
    colnames(scores) <- colnames(loadings)
  }

  # Preserve row names
  if (!is.null(rownames(X))) {
    rownames(scores) <- rownames(X)
  }

  return(scores)
}


#' Calculate Communalities
#'
#' @description
#' Calculates communalities from factor loadings. Communality represents the
#' proportion of variance in each variable that is accounted for by the factors.
#'
#' @param loadings A matrix of factor loadings (p x k)
#'
#' @return A named vector of communalities
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' loadings <- matrix(c(0.8, 0.1, 0.7, 0.2, 0.1, 0.9, 0.2, 0.8), nrow = 4, byrow = TRUE)
#' rownames(loadings) <- paste0("V", 1:4)
#' communalities(loadings)
#' }
communalities <- function(loadings) {
  # Convert to matrix if needed
  loadings <- as.matrix(loadings)

  # Communality is the sum of squared loadings for each variable
  comm <- rowSums(loadings^2)

  # Preserve names
  if (!is.null(rownames(loadings))) {
    names(comm) <- rownames(loadings)
  }

  return(comm)
}


#' Calculate Uniquenesses
#'
#' @description
#' Calculates uniquenesses (specific variances) from factor loadings.
#' Uniqueness represents the proportion of variance in each variable
#' that is NOT accounted for by the factors.
#'
#' @param loadings A matrix of factor loadings (p x k)
#'
#' @return A named vector of uniquenesses
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' loadings <- matrix(c(0.8, 0.1, 0.7, 0.2, 0.1, 0.9, 0.2, 0.8), nrow = 4, byrow = TRUE)
#' rownames(loadings) <- paste0("V", 1:4)
#' uniquenesses(loadings)
#' }
uniquenesses <- function(loadings) {
  comm <- communalities(loadings)
  unique <- 1 - comm
  return(unique)
}


#' Extract Factors Based on Kaiser Criterion
#'
#' @description
#' Extracts the number of factors to retain based on the Kaiser criterion
#' (eigenvalues > mineval, default 1).
#'
#' @param eigenvalues A vector of eigenvalues
#' @param mineval Minimum eigenvalue for retention. Default is 1.
#'
#' @return Number of factors to retain
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' eigenvalues <- c(3.2, 1.5, 0.8, 0.3, 0.2)
#' kaiser_criterion(eigenvalues)
#' }
kaiser_criterion <- function(eigenvalues, mineval = 1) {
  n_factors <- sum(eigenvalues >= mineval)
  if (n_factors == 0) {
    warning("No eigenvalues >= ", mineval, ". Using at least 1 factor.", call. = FALSE)
    n_factors <- 1
  }
  return(n_factors)
}


#' Kaiser-Meyer-Olkin (KMO) Test
#'
#' @description
#' Computes the Kaiser-Meyer-Olkin measure of sampling adequacy
#' for factor analysis.
#'
#' @param R A correlation matrix
#'
#' @return KMO value (0 to 1). Values > 0.6 are generally acceptable for FA.
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' X <- matrix(rnorm(50), nrow = 10)
#' R <- cor(X)
#' kmo_test(R)
#' }
kmo_test <- function(R) {
  # Compute partial correlation matrix
  R_inv <- solve(R)
  k <- ncol(R)

  # Create partial correlation matrix
  partial <- R_inv
  for (i in 1:k) {
    for (j in 1:k) {
      if (i != j) {
        partial[i, j] <- -R_inv[i, j] / sqrt(R_inv[i, i] * R_inv[j, j])
      }
    }
  }

  # Compute KMO
  sum_r2 <- sum((R[!diag(k)])^2)
  sum_q2 <- sum((partial[!diag(k)])^2)

  kmo <- sum_r2 / (sum_r2 + sum_q2)

  return(kmo)
}


#' Measure of Sampling Adequacy (MSA)
#'
#' @description
#' Computes the Measure of Sampling Adequacy for each variable.
#'
#' @param R A correlation matrix
#'
#' @return A named vector of MSA values for each variable
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' X <- matrix(rnorm(50), nrow = 10)
#' colnames(X) <- paste0("V", 1:5)
#' R <- cor(X)
#' msa_individual(R)
#' }
msa_individual <- function(R) {
  # Compute partial correlation matrix
  R_inv <- solve(R)
  k <- ncol(R)

  # Create partial correlation matrix
  partial <- R_inv
  for (i in 1:k) {
    for (j in 1:k) {
      if (i != j) {
        partial[i, j] <- -R_inv[i, j] / sqrt(R_inv[i, i] * R_inv[j, j])
      }
    }
  }

  # Compute MSA for each variable
  msa <- sapply(1:k, function(i) {
    sum_r2 <- sum(R[i, -i]^2)
    sum_q2 <- sum(partial[i, -i]^2)
    sum_r2 / (sum_r2 + sum_q2)
  })

  # Set names
  if (!is.null(colnames(R))) {
    names(msa) <- colnames(R)
  }

  return(msa)
}
