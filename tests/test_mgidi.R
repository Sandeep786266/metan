# Unit Tests for MGIDI Implementation
# ====================================
# Tests for the pure R implementation of MGIDI and related functions.

library(metan)

# Helper function for testing
test_that <- function(description, expr) {
  result <- tryCatch({
    eval(expr)
    TRUE
  }, error = function(e) {
    message("FAILED: ", description)
    message("  Error: ", e$message)
    FALSE
  })
  
  if (result) {
    message("PASSED: ", description)
  }
  
  return(result)
}

expect_equal <- function(x, y, tolerance = 1e-6) {
  if (is.numeric(x) && is.numeric(y)) {
    if (max(abs(x - y)) > tolerance) {
      stop("Values not equal within tolerance")
    }
  } else {
    if (!identical(x, y)) {
      stop("Values not identical")
    }
  }
}

expect_true <- function(x) {
  if (!isTRUE(x)) {
    stop("Expected TRUE but got FALSE")
  }
}

expect_false <- function(x) {
  if (!isFALSE(x)) {
    stop("Expected FALSE but got TRUE")
  }
}

expect_error <- function(expr) {
  result <- tryCatch({
    eval(expr)
    FALSE
  }, error = function(e) TRUE)
  
  if (!result) {
    stop("Expected an error but none occurred")
  }
}

# Run all tests
run_tests <- function() {
  cat("=================================================================\n")
  cat("Running MGIDI Unit Tests\n")
  cat("=================================================================\n\n")
  
  tests_passed <- 0
  tests_failed <- 0
  
  # Test 1: Matrix Multiplication
  cat("--- Matrix Operations Tests ---\n")
  
  if (test_that("matrix_multiply works correctly", {
    A <- matrix(c(1, 2, 3, 4), nrow = 2)
    B <- matrix(c(5, 6, 7, 8), nrow = 2)
    result <- matrix_multiply(A, B)
    expected <- A %*% B
    expect_equal(result, expected)
  })) tests_passed <- tests_passed + 1 else tests_failed <- tests_failed + 1
  
  # Test 2: Matrix Transpose
  if (test_that("matrix_transpose works correctly", {
    A <- matrix(1:6, nrow = 2)
    result <- matrix_transpose(A)
    expected <- t(A)
    expect_equal(result, expected)
  })) tests_passed <- tests_passed + 1 else tests_failed <- tests_failed + 1
  
  # Test 3: Matrix Inverse
  if (test_that("matrix_inverse works correctly", {
    A <- matrix(c(4, 7, 2, 6), nrow = 2)
    result <- matrix_inverse(A)
    expected <- solve(A)
    expect_equal(result, expected, tolerance = 1e-8)
  })) tests_passed <- tests_passed + 1 else tests_failed <- tests_failed + 1
  
  # Test 4: Matrix Standardization
  if (test_that("matrix_standardize works correctly", {
    X <- matrix(c(1, 2, 3, 10, 20, 30), nrow = 3)
    result <- matrix_standardize(X)
    expected <- scale(X)
    expect_equal(as.numeric(result), as.numeric(expected), tolerance = 1e-10)
  })) tests_passed <- tests_passed + 1 else tests_failed <- tests_failed + 1
  
  # Test 5: Correlation Matrix
  if (test_that("correlation_matrix works correctly", {
    set.seed(123)
    X <- matrix(rnorm(30), nrow = 10)
    result <- correlation_matrix(X)
    expected <- cor(X)
    expect_equal(result, expected, tolerance = 1e-10)
  })) tests_passed <- tests_passed + 1 else tests_failed <- tests_failed + 1
  
  # Test 6: Factor Analysis Tests
  cat("\n--- Factor Analysis Tests ---\n")
  
  if (test_that("principal_components returns correct structure", {
    set.seed(123)
    X <- matrix(rnorm(50), nrow = 10)
    colnames(X) <- paste0("V", 1:5)
    result <- principal_components(X)
    expect_true(length(result$eigenvalues) == 5)
    expect_true(all(result$eigenvalues >= 0))
    expect_true(abs(sum(result$variance_explained) - 1) < 1e-10)
  })) tests_passed <- tests_passed + 1 else tests_failed <- tests_failed + 1
  
  if (test_that("communalities are between 0 and 1", {
    loadings <- matrix(c(0.8, 0.1, 0.7, 0.2, 0.1, 0.9, 0.2, 0.8), nrow = 4, byrow = TRUE)
    comm <- communalities(loadings)
    expect_true(all(comm >= 0))
    expect_true(all(comm <= 1.01))  # Small tolerance for numerical precision
  })) tests_passed <- tests_passed + 1 else tests_failed <- tests_failed + 1
  
  if (test_that("varimax_rotation preserves communalities", {
    loadings <- matrix(c(0.8, 0.1, 0.7, 0.2, 0.1, 0.9, 0.2, 0.8), nrow = 4, byrow = TRUE)
    rotated <- varimax_rotation(loadings)
    comm_before <- communalities(loadings)
    comm_after <- communalities(rotated$loadings)
    expect_equal(comm_before, comm_after, tolerance = 1e-4)
  })) tests_passed <- tests_passed + 1 else tests_failed <- tests_failed + 1
  
  if (test_that("kaiser_criterion returns correct number", {
    eigenvalues <- c(3.2, 1.5, 0.8, 0.3, 0.2)
    result <- kaiser_criterion(eigenvalues)
    expect_equal(result, 2)
  })) tests_passed <- tests_passed + 1 else tests_failed <- tests_failed + 1
  
  # Test 7: MGIDI Core Tests
  cat("\n--- MGIDI Core Tests ---\n")
  
  if (test_that("validate_directions works with numeric input", {
    result <- validate_directions(c(1, 1, -1), 3, c("A", "B", "C"))
    expect_equal(result, c(A = 1, B = 1, C = -1))
  })) tests_passed <- tests_passed + 1 else tests_failed <- tests_failed + 1
  
  if (test_that("validate_directions works with character input", {
    result <- validate_directions(c("h", "h", "l"), 3, c("A", "B", "C"))
    expect_equal(result, c(A = 1, B = 1, C = -1))
  })) tests_passed <- tests_passed + 1 else tests_failed <- tests_failed + 1
  
  if (test_that("validate_directions errors on wrong length", {
    expect_error(validate_directions(c(1, 1), 3, c("A", "B", "C")))
  })) tests_passed <- tests_passed + 1 else tests_failed <- tests_failed + 1
  
  if (test_that("mgidi_index runs without errors", {
    set.seed(123)
    data <- data.frame(
      GEN = paste0("G", 1:10),
      V1 = rnorm(10, 5, 1),
      V2 = rnorm(10, 10, 2),
      V3 = rnorm(10, 50, 10)
    )
    result <- mgidi_index(data, c(1, 1, -1), verbose = FALSE)
    expect_true("mgidi_pure" %in% class(result))
  })) tests_passed <- tests_passed + 1 else tests_failed <- tests_failed + 1
  
  if (test_that("mgidi_index returns correct structure", {
    set.seed(123)
    data <- data.frame(
      GEN = paste0("G", 1:10),
      V1 = rnorm(10, 5, 1),
      V2 = rnorm(10, 10, 2),
      V3 = rnorm(10, 50, 10)
    )
    result <- mgidi_index(data, c(1, 1, -1), verbose = FALSE)
    expect_true(!is.null(result$MGIDI))
    expect_true(!is.null(result$data))
    expect_true(!is.null(result$scores_gen))
    expect_true(nrow(result$MGIDI) == 10)
  })) tests_passed <- tests_passed + 1 else tests_failed <- tests_failed + 1
  
  if (test_that("MGIDI values are non-negative", {
    set.seed(123)
    data <- data.frame(
      GEN = paste0("G", 1:10),
      V1 = rnorm(10, 5, 1),
      V2 = rnorm(10, 10, 2),
      V3 = rnorm(10, 50, 10)
    )
    result <- mgidi_index(data, c(1, 1, -1), verbose = FALSE)
    expect_true(all(result$MGIDI$MGIDI >= 0))
  })) tests_passed <- tests_passed + 1 else tests_failed <- tests_failed + 1
  
  if (test_that("MGIDI results are sorted (ascending)", {
    set.seed(123)
    data <- data.frame(
      GEN = paste0("G", 1:10),
      V1 = rnorm(10, 5, 1),
      V2 = rnorm(10, 10, 2),
      V3 = rnorm(10, 50, 10)
    )
    result <- mgidi_index(data, c(1, 1, -1), verbose = FALSE)
    mgidi_vals <- result$MGIDI$MGIDI
    expect_true(all(diff(mgidi_vals) >= 0))
  })) tests_passed <- tests_passed + 1 else tests_failed <- tests_failed + 1
  
  # Test 8: Selection Functions Tests
  cat("\n--- Selection Functions Tests ---\n")
  
  if (test_that("select_genotypes returns correct number", {
    set.seed(123)
    data <- data.frame(
      GEN = paste0("G", 1:20),
      V1 = rnorm(20, 5, 1),
      V2 = rnorm(20, 10, 2)
    )
    result <- mgidi_index(data, c(1, 1), verbose = FALSE)
    selected <- select_genotypes(result, intensity = 0.20)
    expect_equal(length(selected), 4)  # 20% of 20 = 4
  })) tests_passed <- tests_passed + 1 else tests_failed <- tests_failed + 1
  
  if (test_that("select_genotypes with n parameter works", {
    set.seed(123)
    data <- data.frame(
      GEN = paste0("G", 1:20),
      V1 = rnorm(20, 5, 1),
      V2 = rnorm(20, 10, 2)
    )
    result <- mgidi_index(data, c(1, 1), verbose = FALSE)
    selected <- select_genotypes(result, n = 5)
    expect_equal(length(selected), 5)
  })) tests_passed <- tests_passed + 1 else tests_failed <- tests_failed + 1
  
  if (test_that("rank_by_mgidi returns correct structure", {
    distances <- c(G1 = 2.5, G2 = 1.8, G3 = 3.2, G4 = 0.9)
    result <- rank_by_mgidi(distances)
    expect_true("Rank" %in% names(result))
    expect_true("Genotype" %in% names(result))
    expect_equal(result$Genotype[1], "G4")  # Lowest MGIDI
  })) tests_passed <- tests_passed + 1 else tests_failed <- tests_failed + 1
  
  if (test_that("genetic_gain calculates correctly", {
    result <- genetic_gain(5.0, 6.0, heritability = 0.5)
    expect_equal(result$selection_differential, 1.0)
    expect_equal(result$response_to_selection, 0.5)
  })) tests_passed <- tests_passed + 1 else tests_failed <- tests_failed + 1
  
  if (test_that("selection_differential calculates correctly", {
    data <- matrix(1:12, nrow = 4, ncol = 3)
    rownames(data) <- paste0("G", 1:4)
    colnames(data) <- c("V1", "V2", "V3")
    result <- selection_differential(data, c("G3", "G4"), c(1, 1, -1))
    expect_true("SD" %in% names(result))
    expect_true("SD_pct" %in% names(result))
  })) tests_passed <- tests_passed + 1 else tests_failed <- tests_failed + 1
  
  # Test 9: Edge Cases
  cat("\n--- Edge Case Tests ---\n")
  
  if (test_that("mgidi_index handles minimum traits (2)", {
    set.seed(123)
    data <- data.frame(
      GEN = paste0("G", 1:10),
      V1 = rnorm(10),
      V2 = rnorm(10)
    )
    result <- mgidi_index(data, c(1, -1), verbose = FALSE)
    expect_true(!is.null(result$MGIDI))
  })) tests_passed <- tests_passed + 1 else tests_failed <- tests_failed + 1
  
  if (test_that("mgidi_index handles numeric-only data", {
    set.seed(123)
    data <- matrix(rnorm(30), nrow = 10, ncol = 3)
    colnames(data) <- c("V1", "V2", "V3")
    result <- mgidi_index(data, c(1, 1, -1), verbose = FALSE)
    expect_true(!is.null(result$MGIDI))
  })) tests_passed <- tests_passed + 1 else tests_failed <- tests_failed + 1
  
  # Summary
  cat("\n=================================================================\n")
  cat("TEST SUMMARY\n")
  cat("=================================================================\n")
  cat("Tests passed:", tests_passed, "\n")
  cat("Tests failed:", tests_failed, "\n")
  cat("Total tests:", tests_passed + tests_failed, "\n")
  
  if (tests_failed == 0) {
    cat("\nAll tests PASSED!\n")
  } else {
    cat("\nSome tests FAILED. Please review the output above.\n")
  }
  
  return(invisible(list(passed = tests_passed, failed = tests_failed)))
}

# Run tests
if (interactive() || !exists("skip_tests")) {
  run_tests()
}
