#' Pure R Matrix Operations
#'
#' @description
#' `r badge('experimental')`
#'
#' A collection of pure R functions for fundamental matrix operations without
#' relying on external statistical packages for core computations. These functions
#' are designed for educational purposes and to provide transparent implementations
#' of matrix algorithms used in MGIDI calculations.
#'
#' @name matrix_operations
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @md

#' Matrix Multiplication Using Nested Loops
#'
#' @description
#' Performs matrix multiplication using nested loops (pure R implementation).
#'
#' @param A A numeric matrix (m x n)
#' @param B A numeric matrix (n x p)
#'
#' @return A numeric matrix (m x p) resulting from A %*% B
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' A <- matrix(c(1, 2, 3, 4), nrow = 2)
#' B <- matrix(c(5, 6, 7, 8), nrow = 2)
#' matrix_multiply(A, B)
#' }
matrix_multiply <- function(A, B) {
  # Convert to matrix if needed
  A <- as.matrix(A)
  B <- as.matrix(B)

  # Check dimensions compatibility
  if (ncol(A) != nrow(B)) {
    stop("Incompatible matrix dimensions for multiplication. ",
         "ncol(A) = ", ncol(A), " must equal nrow(B) = ", nrow(B), call. = FALSE)
  }

  # Get dimensions
  m <- nrow(A)
  n <- ncol(A)
  p <- ncol(B)

  # Initialize result matrix
  C <- matrix(0, nrow = m, ncol = p)

  # Perform multiplication using nested loops
  for (i in 1:m) {
    for (j in 1:p) {
      for (k in 1:n) {
        C[i, j] <- C[i, j] + A[i, k] * B[k, j]
      }
    }
  }

  # Preserve row and column names if available
  if (!is.null(rownames(A))) {
    rownames(C) <- rownames(A)
  }
  if (!is.null(colnames(B))) {
    colnames(C) <- colnames(B)
  }

  return(C)
}


#' Matrix Transposition
#'
#' @description
#' Transposes a matrix (pure R implementation).
#'
#' @param A A numeric matrix
#'
#' @return The transpose of matrix A
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' A <- matrix(1:6, nrow = 2)
#' matrix_transpose(A)
#' }
matrix_transpose <- function(A) {
  # Convert to matrix if needed
  A <- as.matrix(A)

  # Get dimensions
  m <- nrow(A)
  n <- ncol(A)

  # Initialize result matrix with swapped dimensions
  # Note: This explicit loop implementation is for educational purposes.
  # In practice, R's built-in t() function is more efficient.
  At <- matrix(0, nrow = n, ncol = m)

  # Transpose using explicit indexing
  for (i in 1:m) {
    for (j in 1:n) {
      At[j, i] <- A[i, j]
    }
  }

  # Swap row and column names
  if (!is.null(colnames(A))) {
    rownames(At) <- colnames(A)
  }
  if (!is.null(rownames(A))) {
    colnames(At) <- rownames(A)
  }

  return(At)
}


#' Matrix Inversion Using Gauss-Jordan Elimination
#'
#' @description
#' Computes the inverse of a square matrix using Gauss-Jordan elimination
#' (pure R implementation).
#'
#' @param A A square numeric matrix
#' @param tolerance Tolerance for detecting singular matrices. Default is 1e-10.
#'
#' @return The inverse of matrix A
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' A <- matrix(c(4, 7, 2, 6), nrow = 2)
#' A_inv <- matrix_inverse(A)
#' # Verify: A %*% A_inv should be identity
#' round(A %*% A_inv, 10)
#' }
matrix_inverse <- function(A, tolerance = 1e-10) {
  # Convert to matrix if needed
  A <- as.matrix(A)

  # Check if square
  n <- nrow(A)
  if (n != ncol(A)) {
    stop("Matrix must be square to compute inverse.", call. = FALSE)
  }

  # Create augmented matrix [A | I]
  augmented <- cbind(A, diag(n))

  # Gauss-Jordan elimination
  for (col in 1:n) {
    # Find pivot (largest absolute value in column)
    pivot_row <- which.max(abs(augmented[col:n, col])) + col - 1

    # Check for singular matrix
    if (abs(augmented[pivot_row, col]) < tolerance) {
      stop("Matrix is singular or nearly singular and cannot be inverted.", call. = FALSE)
    }

    # Swap rows if necessary
    if (pivot_row != col) {
      temp <- augmented[col, ]
      augmented[col, ] <- augmented[pivot_row, ]
      augmented[pivot_row, ] <- temp
    }

    # Scale pivot row
    pivot <- augmented[col, col]
    augmented[col, ] <- augmented[col, ] / pivot

    # Eliminate column in other rows
    for (row in 1:n) {
      if (row != col) {
        factor <- augmented[row, col]
        augmented[row, ] <- augmented[row, ] - factor * augmented[col, ]
      }
    }
  }

  # Extract inverse from augmented matrix
  A_inv <- augmented[, (n + 1):(2 * n)]

  # Preserve row and column names
  if (!is.null(rownames(A))) {
    rownames(A_inv) <- rownames(A)
  }
  if (!is.null(colnames(A))) {
    colnames(A_inv) <- colnames(A)
  }

  return(A_inv)
}


#' Matrix Standardization (Z-scores)
#'
#' @description
#' Standardizes a matrix by computing Z-scores for each column
#' (subtracting mean and dividing by standard deviation).
#'
#' @param X A numeric matrix or data frame
#' @param center Logical. If TRUE, center the data by subtracting column means.
#'   Default is TRUE.
#' @param scale_by_sd Logical. If TRUE, scale by dividing by column standard
#'   deviations. Default is TRUE.
#' @param na.rm Logical. If TRUE, NA values are removed when computing
#'   means and standard deviations. Default is TRUE.
#'
#' @return A standardized matrix with the same dimensions as X
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' X <- matrix(c(1, 2, 3, 10, 20, 30), nrow = 3)
#' colnames(X) <- c("A", "B")
#' matrix_standardize(X)
#' }
matrix_standardize <- function(X, center = TRUE, scale_by_sd = TRUE, na.rm = TRUE) {
  # Convert to matrix if needed
  X <- as.matrix(X)

  n <- nrow(X)
  p <- ncol(X)

  # Initialize result matrix
  Z <- matrix(0, nrow = n, ncol = p)

  for (j in 1:p) {
    col_data <- X[, j]

    # Calculate mean
    col_mean <- if (center) mean(col_data, na.rm = na.rm) else 0

    # Calculate standard deviation
    col_sd <- if (scale_by_sd) {
      sqrt(sum((col_data[!is.na(col_data)] - mean(col_data, na.rm = na.rm))^2,
               na.rm = na.rm) / (sum(!is.na(col_data)) - 1))
    } else {
      1
    }

    # Avoid division by zero
    if (col_sd < 1e-10) {
      col_sd <- 1
    }

    # Standardize
    Z[, j] <- (col_data - col_mean) / col_sd
  }

  # Preserve row and column names
  if (!is.null(rownames(X))) {
    rownames(Z) <- rownames(X)
  }
  if (!is.null(colnames(X))) {
    colnames(Z) <- colnames(X)
  }

  return(Z)
}


#' Correlation Matrix Computation
#'
#' @description
#' Computes the Pearson correlation matrix for a data matrix
#' (pure R implementation).
#'
#' @param X A numeric matrix or data frame where rows are observations
#'   and columns are variables.
#' @param use Method for handling missing values. Options are "complete.obs"
#'   (default) or "pairwise.complete.obs".
#'
#' @return A correlation matrix (p x p) where p is the number of columns in X
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' X <- matrix(rnorm(30), nrow = 10)
#' colnames(X) <- c("A", "B", "C")
#' correlation_matrix(X)
#' }
correlation_matrix <- function(X, use = "complete.obs") {
  # Convert to matrix if needed
  X <- as.matrix(X)

  # Handle missing values
  if (use == "complete.obs") {
    complete_rows <- complete.cases(X)
    X <- X[complete_rows, , drop = FALSE]
  }

  n <- nrow(X)
  p <- ncol(X)

  # Initialize correlation matrix
  R <- matrix(1, nrow = p, ncol = p)

  # Compute correlations
  for (i in 1:p) {
    for (j in 1:p) {
      if (i != j) {
        if (use == "pairwise.complete.obs") {
          # Use only complete pairs
          valid <- !is.na(X[, i]) & !is.na(X[, j])
          xi <- X[valid, i]
          xj <- X[valid, j]
        } else {
          xi <- X[, i]
          xj <- X[, j]
        }

        # Compute means
        mean_i <- mean(xi)
        mean_j <- mean(xj)

        # Compute correlation
        num <- sum((xi - mean_i) * (xj - mean_j))
        den_i <- sqrt(sum((xi - mean_i)^2))
        den_j <- sqrt(sum((xj - mean_j)^2))

        if (den_i > 0 && den_j > 0) {
          R[i, j] <- num / (den_i * den_j)
        } else {
          R[i, j] <- NA
        }
      }
    }
  }

  # Set names
  if (!is.null(colnames(X))) {
    rownames(R) <- colnames(X)
    colnames(R) <- colnames(X)
  }

  return(R)
}


#' Eigendecomposition Using Power Iteration Method
#'
#' @description
#' Computes eigenvalues and eigenvectors of a symmetric matrix using the
#' power iteration method with deflation (pure R implementation).
#'
#' @param A A symmetric numeric matrix
#' @param max_iter Maximum number of iterations for power iteration. Default is 1000.
#' @param tolerance Convergence tolerance. Default is 1e-10.
#'
#' @return A list with components:
#' \itemize{
#'   \item values: A vector of eigenvalues in decreasing order
#'   \item vectors: A matrix whose columns are the corresponding eigenvectors
#' }
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' A <- matrix(c(4, 1, 1, 3), nrow = 2)
#' result <- eigen_decomposition(A)
#' result$values
#' result$vectors
#' }
eigen_decomposition <- function(A, max_iter = 1000, tolerance = 1e-10) {
  # Convert to matrix if needed
  A <- as.matrix(A)

  # Check if square
  n <- nrow(A)
  if (n != ncol(A)) {
    stop("Matrix must be square for eigendecomposition.", call. = FALSE)
  }

  # Initialize storage for eigenvalues and eigenvectors
  eigenvalues <- numeric(n)
  eigenvectors <- matrix(0, nrow = n, ncol = n)

  # Working copy of matrix
  A_work <- A

  for (k in 1:n) {
    # Initialize random vector
    v <- rep(1, n)
    v <- v / sqrt(sum(v^2))

    # Power iteration
    for (iter in 1:max_iter) {
      # Matrix-vector multiplication
      v_new <- A_work %*% v

      # Compute eigenvalue (Rayleigh quotient)
      lambda <- sum(v * v_new)

      # Normalize
      norm_v <- sqrt(sum(v_new^2))
      if (norm_v < tolerance) {
        break
      }
      v_new <- v_new / norm_v

      # Check convergence
      if (max(abs(v_new - v)) < tolerance || max(abs(v_new + v)) < tolerance) {
        break
      }
      v <- v_new
    }

    # Store eigenvalue and eigenvector
    eigenvalues[k] <- lambda
    eigenvectors[, k] <- v

    # Deflation: remove the contribution of the found eigenpair
    A_work <- A_work - lambda * (v %*% t(v))
  }

  # Sort by decreasing eigenvalue
  ord <- order(eigenvalues, decreasing = TRUE)
  eigenvalues <- eigenvalues[ord]
  eigenvectors <- eigenvectors[, ord, drop = FALSE]

  # Set column names
  colnames(eigenvectors) <- paste0("PC", 1:n)
  if (!is.null(rownames(A))) {
    rownames(eigenvectors) <- rownames(A)
  }

  return(list(values = eigenvalues, vectors = eigenvectors))
}
