#' MGIDI Visualization Functions
#'
#' @description
#' `r badge('experimental')`
#'
#' Visualization functions for MGIDI analysis results.
#'
#' @name mgidi_plots
#' @author Tiago Olivoto \email{tiagoolivoto@@gmail.com}
#' @md

#' Plot MGIDI Ranking
#'
#' @description
#' Creates a bar plot showing MGIDI values for all genotypes.
#'
#' @param mgidi_result An object of class `mgidi_pure` or `mgidi`
#' @param n Number of genotypes to display. Default is NULL (all genotypes).
#' @param selected_only Logical. If TRUE, only show selected genotypes.
#' @param selection_intensity Selection intensity (0-1) for highlighting
#'   selected genotypes. Default is 0.20.
#' @param col_selected Color for selected genotypes. Default is "red".
#' @param col_nonselected Color for non-selected genotypes. Default is "gray".
#' @param title Plot title. Default is "MGIDI Ranking".
#' @param x_lab X-axis label.
#' @param y_lab Y-axis label.
#' @param rotate_x Logical. If TRUE, rotate x-axis labels. Default is TRUE.
#'
#' @return A ggplot object
#' @export
#' @import ggplot2
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
#' plot_mgidi_ranking(result)
#' }
plot_mgidi_ranking <- function(mgidi_result,
                               n = NULL,
                               selected_only = FALSE,
                               selection_intensity = 0.20,
                               col_selected = "red",
                               col_nonselected = "gray",
                               title = "MGIDI Ranking",
                               x_lab = "Genotypes",
                               y_lab = "MGIDI",
                               rotate_x = TRUE) {

  # Extract MGIDI data
  mgidi_df <- mgidi_result$MGIDI
  gen_col <- names(mgidi_df)[1]
  mgidi_col <- names(mgidi_df)[2]

  total_gen <- nrow(mgidi_df)
  n_selected <- round(total_gen * selection_intensity)

  # Add selection status
  mgidi_df$Status <- ifelse(1:total_gen <= n_selected, "Selected", "Not Selected")

  # Filter if needed
  if (selected_only) {
    mgidi_df <- mgidi_df[mgidi_df$Status == "Selected", ]
  } else if (!is.null(n)) {
    mgidi_df <- mgidi_df[1:min(n, total_gen), ]
  }

  # Reorder factor for plotting
  mgidi_df[[gen_col]] <- factor(mgidi_df[[gen_col]],
                                levels = rev(mgidi_df[[gen_col]]))

  # Create plot
  p <- ggplot(mgidi_df, aes(x = .data[[gen_col]], y = .data[[mgidi_col]], fill = .data[["Status"]])) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c("Selected" = col_selected,
                                 "Not Selected" = col_nonselected)) +
    labs(title = title, x = x_lab, y = y_lab) +
    coord_flip() +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )

  if (rotate_x && !selected_only && is.null(n)) {
    p <- p + coord_flip()
  }

  return(p)
}


#' Plot Factor Contributions
#'
#' @description
#' Creates a plot showing the contribution of each factor to the MGIDI
#' for selected genotypes.
#'
#' @param mgidi_result An object of class `mgidi_pure` or `mgidi`
#' @param genotypes Character vector of genotypes to display. Default is "selected"
#'   which shows only selected genotypes based on selection_intensity.
#' @param selection_intensity Selection intensity for determining selected
#'   genotypes. Default is 0.20.
#' @param type Plot type: "stacked" for stacked bar chart, "radar" for radar plot.
#'   Default is "stacked".
#' @param position Bar position: "fill" (relative proportions) or "stack"
#'   (absolute values). Default is "fill".
#' @param title Plot title.
#' @param x_lab X-axis label.
#' @param y_lab Y-axis label.
#'
#' @return A ggplot object
#' @export
#' @import ggplot2
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
#' plot_factor_contributions(result)
#' }
plot_factor_contributions <- function(mgidi_result,
                                      genotypes = "selected",
                                      selection_intensity = 0.20,
                                      type = "stacked",
                                      position = "fill",
                                      title = "Factor Contributions to MGIDI",
                                      x_lab = "Genotypes",
                                      y_lab = "Contribution") {

  # Extract contributions
  contri_df <- mgidi_result$contri_fac
  gen_col <- names(contri_df)[1]
  factor_cols <- names(contri_df)[-1]

  # Determine which genotypes to show
  if (identical(genotypes, "selected")) {
    total_gen <- nrow(mgidi_result$MGIDI)
    n_selected <- max(1, round(total_gen * selection_intensity))
    selected_gen <- mgidi_result$MGIDI[[1]][1:n_selected]
    contri_df <- contri_df[contri_df[[gen_col]] %in% selected_gen, ]
  } else if (identical(genotypes, "all")) {
    # Use all genotypes
  } else if (is.character(genotypes)) {
    contri_df <- contri_df[contri_df[[gen_col]] %in% genotypes, ]
  }

  # Reshape for ggplot - using base R to avoid external dependencies
  # This is a simple reshape operation
  n_rows <- nrow(contri_df)
  n_factors <- length(factor_cols)
  contri_long <- data.frame(
    stringsAsFactors = FALSE
  )
  for (i in 1:n_rows) {
    for (j in 1:n_factors) {
      contri_long <- rbind(contri_long, data.frame(
        GEN_temp = contri_df[[gen_col]][i],
        Factor = factor_cols[j],
        Contribution = contri_df[[factor_cols[j]]][i],
        stringsAsFactors = FALSE
      ))
    }
  }
  names(contri_long)[1] <- gen_col

  # Create plot
  if (type == "stacked") {
    p <- ggplot(contri_long, aes(x = .data[[gen_col]], y = .data[["Contribution"]],
                                        fill = .data[["Factor"]])) +
      geom_bar(stat = "identity", position = position,
               color = "black", linewidth = 0.2) +
      labs(title = title, x = x_lab, y = y_lab) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      scale_y_continuous(expand = c(0, 0))
  } else if (type == "radar") {
    # Radar plot (simplified version)
    p <- ggplot(contri_long, aes(x = .data[[gen_col]], y = .data[["Contribution"]],
                                        group = .data[["Factor"]], color = .data[["Factor"]])) +
      geom_polygon(aes(fill = .data[["Factor"]]), alpha = 0.2) +
      geom_line(linewidth = 0.8) +
      coord_polar() +
      labs(title = title) +
      theme_minimal() +
      theme(legend.position = "bottom")
  }

  return(p)
}


#' Plot Selection Gains
#'
#' @description
#' Creates a plot showing genetic gains for each trait from selection.
#'
#' @param mgidi_result An object of class `mgidi_pure` or `mgidi`
#' @param selected_genotypes Character vector of selected genotypes. If NULL,
#'   uses selection_intensity.
#' @param selection_intensity Selection intensity if selected_genotypes is NULL.
#'   Default is 0.20.
#' @param show_percentage Logical. If TRUE, show percentage gains. If FALSE,
#'   show absolute gains. Default is TRUE.
#' @param title Plot title.
#' @param x_lab X-axis label.
#' @param y_lab Y-axis label.
#' @param col_positive Color for positive (desired) gains.
#' @param col_negative Color for negative (undesired) gains.
#'
#' @return A ggplot object
#' @export
#' @import ggplot2
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
#' plot_selection_gains(result)
#' }
plot_selection_gains <- function(mgidi_result,
                                 selected_genotypes = NULL,
                                 selection_intensity = 0.20,
                                 show_percentage = TRUE,
                                 title = "Selection Gains by Trait",
                                 x_lab = "Traits",
                                 y_lab = NULL,
                                 col_positive = "forestgreen",
                                 col_negative = "tomato") {

  # Determine selected genotypes
  if (is.null(selected_genotypes)) {
    total_gen <- nrow(mgidi_result$MGIDI)
    n_selected <- max(1, round(total_gen * selection_intensity))
    selected_genotypes <- mgidi_result$MGIDI[[1]][1:n_selected]
  }

  # Get trait data and directions
  trait_data <- mgidi_result$data
  directions <- mgidi_result$directions

  # Compute selection differential
  pop_means <- colMeans(trait_data, na.rm = TRUE)
  sel_means <- colMeans(trait_data[selected_genotypes, , drop = FALSE], na.rm = TRUE)
  SD <- sel_means - pop_means
  SD_pct <- (SD / abs(pop_means)) * 100

  # Create data frame for plotting
  plot_df <- data.frame(
    Trait = names(SD),
    SD = SD,
    SD_pct = SD_pct,
    Direction = directions,
    stringsAsFactors = FALSE
  )

  # Determine if gain is in desired direction
  plot_df$Desired <- ifelse(
    (plot_df$Direction == 1 & plot_df$SD > 0) |
      (plot_df$Direction == -1 & plot_df$SD < 0),
    "Desired", "Undesired"
  )

  # Select which value to plot
  if (show_percentage) {
    y_var <- "SD_pct"
    if (is.null(y_lab)) y_lab <- "Selection Gain (%)"
  } else {
    y_var <- "SD"
    if (is.null(y_lab)) y_lab <- "Selection Gain (absolute)"
  }

  # Create plot
  p <- ggplot(plot_df, aes(x = .data[["Trait"]], y = .data[[y_var]], fill = .data[["Desired"]])) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_fill_manual(values = c("Desired" = col_positive,
                                 "Undesired" = col_negative)) +
    labs(title = title, x = x_lab, y = y_lab) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  return(p)
}


#' Biplot for MGIDI Factor Scores
#'
#' @description
#' Creates a biplot showing genotypes in factor score space.
#'
#' @param mgidi_result An object of class `mgidi_pure` or `mgidi`
#' @param factors Which factors to plot. Default is c(1, 2).
#' @param selection_intensity Selection intensity for highlighting. Default is 0.20.
#' @param show_loadings Logical. If TRUE, show loading vectors. Default is TRUE.
#' @param show_ideotype Logical. If TRUE, show ideotype position. Default is TRUE.
#' @param col_selected Color for selected genotypes.
#' @param col_nonselected Color for non-selected genotypes.
#' @param col_ideotype Color for ideotype point.
#' @param point_size Size of genotype points.
#' @param label_size Size of labels.
#' @param title Plot title.
#'
#' @return A ggplot object
#' @export
#' @import ggplot2
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
#' biplot_mgidi(result)
#' }
biplot_mgidi <- function(mgidi_result,
                         factors = c(1, 2),
                         selection_intensity = 0.20,
                         show_loadings = TRUE,
                         show_ideotype = TRUE,
                         col_selected = "red",
                         col_nonselected = "gray60",
                         col_ideotype = "blue",
                         point_size = 3,
                         label_size = 3,
                         title = "MGIDI Biplot") {

  # Check number of factors
  n_factors <- ncol(mgidi_result$scores_gen)
  if (n_factors < 2) {
    stop("Biplot requires at least 2 factors. Only ", n_factors, " factor(s) available.",
         call. = FALSE)
  }

  # Extract scores
  scores <- mgidi_result$scores_gen
  loadings <- mgidi_result$rotated_loadings

  # Determine selected genotypes
  total_gen <- nrow(scores)
  n_selected <- max(1, round(total_gen * selection_intensity))
  selected_gen <- mgidi_result$MGIDI[[1]][1:n_selected]

  # Create data frame for genotypes
  gen_df <- data.frame(
    Genotype = rownames(scores),
    Factor1 = scores[, factors[1]],
    Factor2 = scores[, factors[2]],
    stringsAsFactors = FALSE
  )
  gen_df$Status <- ifelse(gen_df$Genotype %in% selected_gen, "Selected", "Not Selected")

  # Create plot
  p <- ggplot() +
    # Add genotype points
    geom_point(data = gen_df,
               aes(x = .data[["Factor1"]], y = .data[["Factor2"]], color = .data[["Status"]]),
               size = point_size, alpha = 0.8) +
    scale_color_manual(values = c("Selected" = col_selected,
                                  "Not Selected" = col_nonselected)) +
    # Add labels - use ggrepel if available, otherwise use geom_text
    {
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        ggrepel::geom_text_repel(data = gen_df,
                                 aes(x = .data[["Factor1"]], y = .data[["Factor2"]], label = .data[["Genotype"]]),
                                 size = label_size, max.overlaps = 15)
      } else {
        geom_text(data = gen_df,
                  aes(x = .data[["Factor1"]], y = .data[["Factor2"]], label = .data[["Genotype"]]),
                  size = label_size, vjust = -0.5)
      }
    } +
    # Add axis lines
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = title,
      x = paste0("Factor ", factors[1]),
      y = paste0("Factor ", factors[2])
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      aspect.ratio = 1
    )

  # Add ideotype
  if (show_ideotype && !is.null(mgidi_result$scores_ide)) {
    ide_scores <- mgidi_result$scores_ide
    ide_df <- data.frame(
      x = ide_scores[, factors[1]],
      y = ide_scores[, factors[2]]
    )
    p <- p +
      geom_point(data = ide_df, aes(x = x, y = y),
                 shape = 17, size = point_size * 1.5, color = col_ideotype) +
      annotate("text", x = ide_df$x, y = ide_df$y,
               label = "IDEOTYPE", vjust = -1.5, color = col_ideotype,
               size = label_size, fontface = "bold")
  }

  # Add loadings as arrows
  if (show_loadings) {
    # Scale loadings for visibility
    scale_factor <- max(abs(gen_df$Factor1), abs(gen_df$Factor2)) /
      max(abs(loadings[, factors]))
    load_df <- data.frame(
      Trait = rownames(loadings),
      x = loadings[, factors[1]] * scale_factor * 0.8,
      y = loadings[, factors[2]] * scale_factor * 0.8
    )

    p <- p +
      geom_segment(data = load_df,
                   aes(x = 0, y = 0, xend = x, yend = y),
                   arrow = arrow(length = unit(0.2, "cm")),
                   color = "darkgreen", linewidth = 0.5) +
      geom_text(data = load_df,
                aes(x = x * 1.1, y = y * 1.1, label = Trait),
                color = "darkgreen", size = label_size, fontface = "bold")
  }

  return(p)
}


#' Print Method for mgidi_pure Objects
#'
#' @description
#' Print method for objects of class `mgidi_pure`.
#'
#' @param x An object of class `mgidi_pure`
#' @param n Number of genotypes to display. Default is 10.
#' @param ... Additional arguments (not used)
#'
#' @return Invisible x
#' @export
#' @method print mgidi_pure
print.mgidi_pure <- function(x, n = 10, ...) {
  cat("MGIDI Analysis Results\n")
  cat("======================\n\n")

  cat("Number of genotypes:", nrow(x$data), "\n")
  cat("Number of traits:", ncol(x$data), "\n")
  cat("Number of factors:", x$n_factors, "\n")
  cat("KMO:", round(x$KMO, 4), "\n\n")

  cat("Traits and directions:\n")
  dir_labels <- ifelse(x$directions == 1, "Maximize", "Minimize")
  for (i in seq_along(x$directions)) {
    cat("  ", names(x$directions)[i], ":", dir_labels[i], "\n")
  }

  cat("\nTop", min(n, nrow(x$MGIDI)), "genotypes by MGIDI:\n")
  print(head(x$MGIDI, n), row.names = FALSE)

  invisible(x)
}
