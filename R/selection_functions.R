#' Selection and Ranking Functions
#'
#' @description
#' `r badge('experimental')`
#'
#' Functions for selecting and ranking genotypes based on MGIDI results.
#'
#' @name selection_functions
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @md

#' Select Genotypes Based on MGIDI
#'
#' @description
#' Selects the top genotypes based on MGIDI values (lower is better).
#'
#' @param mgidi_result An object of class `mgidi_pure` or `mgidi`
#' @param intensity Selection intensity as a proportion (0-1). Default is 0.20
#'   (select top 20%).
#' @param n Number of genotypes to select. If provided, overrides `intensity`.
#'
#' @return A character vector of selected genotype names
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' # Create example data
#' set.seed(123)
#' data <- data.frame(
#'   GEN = paste0("G", 1:20),
#'   GY = rnorm(20, 5, 1),
#'   HI = rnorm(20, 0.4, 0.05),
#'   PH = rnorm(20, 100, 10)
#' )
#' directions <- c(GY = 1, HI = 1, PH = -1)
#' result <- mgidi_index(data, directions, verbose = FALSE)
#' select_genotypes(result, intensity = 0.25)
#' }
select_genotypes <- function(mgidi_result, intensity = 0.20, n = NULL) {
  # Extract MGIDI data frame
  if (has_class(mgidi_result, "mgidi_pure")) {
    mgidi_df <- mgidi_result$MGIDI
    gen_col <- "Genotype"
  } else if (has_class(mgidi_result, "mgidi")) {
    mgidi_df <- mgidi_result$MGIDI
    gen_col <- "Genotype"
  } else {
    stop("Input must be an mgidi_pure or mgidi object.", call. = FALSE)
  }

  total_gen <- nrow(mgidi_df)

  # Determine number to select
  if (!is.null(n)) {
    n_select <- min(n, total_gen)
  } else {
    n_select <- max(1, round(total_gen * intensity))
  }

  # Select top genotypes (MGIDI is already sorted, lower is better)
  selected <- mgidi_df[[gen_col]][1:n_select]

  return(selected)
}


#' Rank Genotypes by MGIDI
#'
#' @description
#' Ranks genotypes by their MGIDI values.
#'
#' @param mgidi_result An object of class `mgidi_pure` or `mgidi`, or a named
#'   numeric vector of MGIDI distances
#'
#' @return A data frame with genotypes ranked by MGIDI (rank 1 = best)
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' # Using a named vector
#' distances <- c(G1 = 2.5, G2 = 1.8, G3 = 3.2, G4 = 0.9, G5 = 2.1)
#' rank_by_mgidi(distances)
#' }
rank_by_mgidi <- function(mgidi_result) {
  # Handle different input types
  if (has_class(mgidi_result, c("mgidi_pure", "mgidi"))) {
    mgidi_df <- mgidi_result$MGIDI
    if ("Genotype" %in% names(mgidi_df)) {
      genotypes <- mgidi_df$Genotype
      values <- mgidi_df$MGIDI
    } else {
      genotypes <- mgidi_df[[1]]
      values <- mgidi_df[[2]]
    }
  } else if (is.numeric(mgidi_result)) {
    genotypes <- names(mgidi_result)
    if (is.null(genotypes)) {
      genotypes <- paste0("G", seq_along(mgidi_result))
    }
    values <- as.numeric(mgidi_result)
  } else {
    stop("Input must be an mgidi object or a numeric vector.", call. = FALSE)
  }

  # Create ranked data frame
  ord <- order(values)
  result <- data.frame(
    Rank = 1:length(values),
    Genotype = genotypes[ord],
    MGIDI = values[ord],
    stringsAsFactors = FALSE
  )

  return(result)
}


#' Calculate Genetic Gain
#'
#' @description
#' Calculates breeding progress (genetic gain) from selection.
#'
#' @param population_mean Mean of the original population
#' @param selected_mean Mean of the selected individuals
#' @param heritability Optional heritability estimate. If provided, computes
#'   expected response to selection.
#'
#' @return A list with selection differential and genetic gain
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' genetic_gain(population_mean = 5.0, selected_mean = 6.2, heritability = 0.6)
#' }
genetic_gain <- function(population_mean, selected_mean, heritability = NULL) {
  # Selection differential
  S <- selected_mean - population_mean

  # Percentage
  S_pct <- (S / abs(population_mean)) * 100

  result <- list(
    population_mean = population_mean,
    selected_mean = selected_mean,
    selection_differential = S,
    selection_differential_pct = S_pct
  )

  # Expected response to selection (if heritability provided)
  if (!is.null(heritability)) {
    R <- heritability * S
    R_pct <- (R / abs(population_mean)) * 100
    result$heritability <- heritability
    result$response_to_selection <- R
    result$response_to_selection_pct <- R_pct
  }

  return(result)
}


#' Factor Contributions to MGIDI
#'
#' @description
#' Computes the relative contribution of each factor to the MGIDI index
#' for each genotype.
#'
#' @param mgidi_result An object of class `mgidi_pure` or `mgidi`
#' @param genotypes Optional character vector of specific genotypes to analyze.
#'   If NULL, all genotypes are analyzed.
#'
#' @return A data frame with factor contributions for each genotype
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' set.seed(123)
#' data <- data.frame(
#'   GEN = paste0("G", 1:20),
#'   GY = rnorm(20, 5, 1),
#'   HI = rnorm(20, 0.4, 0.05),
#'   PH = rnorm(20, 100, 10)
#' )
#' directions <- c(GY = 1, HI = 1, PH = -1)
#' result <- mgidi_index(data, directions, verbose = FALSE)
#' contribution_factors(result)
#' }
contribution_factors <- function(mgidi_result, genotypes = NULL) {
  if (!has_class(mgidi_result, c("mgidi_pure", "mgidi"))) {
    stop("Input must be an mgidi_pure or mgidi object.", call. = FALSE)
  }

  contri_fac <- mgidi_result$contri_fac

  # Filter genotypes if specified
  if (!is.null(genotypes)) {
    gen_col <- names(contri_fac)[1]  # Usually "GEN" or "Genotype"
    contri_fac <- contri_fac[contri_fac[[gen_col]] %in% genotypes, , drop = FALSE]
  }

  return(contri_fac)
}


#' Summary Statistics for Selected Genotypes
#'
#' @description
#' Computes summary statistics comparing selected vs unselected genotypes.
#'
#' @param mgidi_result An object of class `mgidi_pure` or `mgidi`
#' @param selected_genotypes Character vector of selected genotype names
#'
#' @return A data frame with trait summaries for selected and unselected groups
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' set.seed(123)
#' data <- data.frame(
#'   GEN = paste0("G", 1:20),
#'   GY = rnorm(20, 5, 1),
#'   HI = rnorm(20, 0.4, 0.05),
#'   PH = rnorm(20, 100, 10)
#' )
#' directions <- c(GY = 1, HI = 1, PH = -1)
#' result <- mgidi_index(data, directions, verbose = FALSE)
#' selected <- select_genotypes(result, intensity = 0.25)
#' summary_selected(result, selected)
#' }
summary_selected <- function(mgidi_result, selected_genotypes) {
  # Get original data
  data <- mgidi_result$data

  # Split into selected and unselected
  all_gen <- rownames(data)
  unselected <- setdiff(all_gen, selected_genotypes)

  sel_data <- data[selected_genotypes, , drop = FALSE]
  unsel_data <- data[unselected, , drop = FALSE]

  # Compute summary statistics
  traits <- colnames(data)
  n_traits <- length(traits)

  result <- data.frame(
    Trait = traits,
    Overall_Mean = colMeans(data, na.rm = TRUE),
    Overall_SD = apply(data, 2, sd, na.rm = TRUE),
    Selected_Mean = colMeans(sel_data, na.rm = TRUE),
    Selected_SD = apply(sel_data, 2, sd, na.rm = TRUE),
    Unselected_Mean = colMeans(unsel_data, na.rm = TRUE),
    Unselected_SD = apply(unsel_data, 2, sd, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  # Add selection differential
  result$SD <- result$Selected_Mean - result$Overall_Mean
  result$SD_pct <- result$SD / abs(result$Overall_Mean) * 100

  return(result)
}


#' Get Selected Genotype Details
#'
#' @description
#' Returns detailed information about selected genotypes including their
#' trait values, factor scores, and MGIDI values.
#'
#' @param mgidi_result An object of class `mgidi_pure` or `mgidi`
#' @param selected_genotypes Character vector of selected genotype names
#'
#' @return A list with:
#' \itemize{
#'   \item trait_values: Original trait values for selected genotypes
#'   \item rescaled_values: Rescaled trait values
#'   \item factor_scores: Factor scores
#'   \item mgidi_values: MGIDI values
#'   \item factor_contributions: Factor contributions to MGIDI
#' }
#' @export
#' @examples
#' \donttest{
#' library(metan)
#' set.seed(123)
#' data <- data.frame(
#'   GEN = paste0("G", 1:20),
#'   GY = rnorm(20, 5, 1),
#'   HI = rnorm(20, 0.4, 0.05),
#'   PH = rnorm(20, 100, 10)
#' )
#' directions <- c(GY = 1, HI = 1, PH = -1)
#' result <- mgidi_index(data, directions, verbose = FALSE)
#' selected <- select_genotypes(result, intensity = 0.25)
#' get_selected_details(result, selected)
#' }
get_selected_details <- function(mgidi_result, selected_genotypes) {
  result <- list()

  # Trait values
  result$trait_values <- mgidi_result$data[selected_genotypes, , drop = FALSE]

  # Rescaled values
  if (!is.null(mgidi_result$rescaled_data)) {
    result$rescaled_values <- mgidi_result$rescaled_data[selected_genotypes, , drop = FALSE]
  }

  # Factor scores
  if (!is.null(mgidi_result$scores_gen)) {
    result$factor_scores <- mgidi_result$scores_gen[selected_genotypes, , drop = FALSE]
  }

  # MGIDI values
  mgidi_df <- mgidi_result$MGIDI
  gen_col <- names(mgidi_df)[1]
  result$mgidi_values <- mgidi_df[mgidi_df[[gen_col]] %in% selected_genotypes, ]

  # Factor contributions
  if (!is.null(mgidi_result$contri_fac)) {
    contri_df <- mgidi_result$contri_fac
    gen_col <- names(contri_df)[1]
    result$factor_contributions <- contri_df[contri_df[[gen_col]] %in% selected_genotypes, ]
  }

  return(result)
}
