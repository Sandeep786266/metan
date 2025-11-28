# Example: Wheat Breeding MGIDI Analysis
# =======================================
# This example demonstrates the use of MGIDI for selecting superior wheat
# genotypes based on multiple traits.

library(metan)

# Create simulated wheat breeding data
# ------------------------------------
set.seed(2024)
n_genotypes <- 30

wheat_data <- data.frame(
  GEN = paste0("W", sprintf("%02d", 1:n_genotypes)),
  
  # Grain yield (kg/ha) - MAXIMIZE
  GY = rnorm(n_genotypes, mean = 5500, sd = 800),
  
  # Harvest index (ratio) - MAXIMIZE
  HI = rnorm(n_genotypes, mean = 0.42, sd = 0.05),
  
  # Thousand kernel weight (g) - MAXIMIZE
  TKW = rnorm(n_genotypes, mean = 38, sd = 5),
  
  # Protein content (%) - MAXIMIZE
  PROT = rnorm(n_genotypes, mean = 12.5, sd = 1.2),
  
  # Plant height (cm) - MINIMIZE (for lodging resistance)
  PH = rnorm(n_genotypes, mean = 95, sd = 12),
  
  # Days to heading (days) - MINIMIZE (for early maturity)
  DTH = rnorm(n_genotypes, mean = 75, sd = 5)
)

# Display the first few rows
cat("Wheat Breeding Data (first 10 genotypes):\n")
print(head(wheat_data, 10))

# Define breeding objectives (ideotype directions)
# ------------------------------------------------
# +1 = Maximize (higher values are better)
# -1 = Minimize (lower values are better)

ideotype <- c(
  GY = 1,    # Maximize grain yield
  HI = 1,    # Maximize harvest index
  TKW = 1,   # Maximize kernel weight
  PROT = 1,  # Maximize protein content
  PH = -1,   # Minimize plant height
  DTH = -1   # Minimize days to heading
)

cat("\nBreeding Objectives:\n")
for (trait in names(ideotype)) {
  direction <- if (ideotype[trait] == 1) "Maximize" else "Minimize"
  cat(sprintf("  %s: %s\n", trait, direction))
}

# Run MGIDI analysis
# ------------------
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("Running MGIDI Analysis\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

result <- mgidi_index(
  data = wheat_data,
  ideotype_directions = ideotype,
  n_factors = NULL,  # Auto-determine using Kaiser criterion
  verbose = TRUE
)

# Select top 20% of genotypes
# ---------------------------
selected <- select_genotypes(result, intensity = 0.20)
cat("\nSelected genotypes (top 20%):", paste(selected, collapse = ", "), "\n")

# Compute selection differential
# ------------------------------
sel_diff <- selection_differential(
  original_data = result$data,
  selected_genotypes = selected,
  trait_directions = ideotype
)

cat("\nSelection Differential:\n")
print(round(sel_diff, 3))

# Summary statistics
# ------------------
summary_stats <- summary_selected(result, selected)
cat("\nSummary Statistics:\n")
print(round(summary_stats, 3))

# Factor contributions for selected genotypes
# -------------------------------------------
cat("\nFactor Contributions for Selected Genotypes:\n")
contributions <- contribution_factors(result, selected)
print(round(contributions, 2))

# Visualization (if running interactively)
# ----------------------------------------
if (interactive()) {
  # MGIDI ranking plot
  p1 <- plot_mgidi_ranking(result, selection_intensity = 0.20)
  print(p1)
  
  # Factor contribution plot
  p2 <- plot_factor_contributions(result, selection_intensity = 0.20)
  print(p2)
  
  # Selection gains plot
  p3 <- plot_selection_gains(result, selection_intensity = 0.20)
  print(p3)
  
  # Biplot
  p4 <- biplot_mgidi(result, selection_intensity = 0.20)
  print(p4)
}

# Print final recommendations
# ---------------------------
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("RECOMMENDATIONS\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

cat("\nBased on MGIDI analysis, the following genotypes are recommended\n")
cat("for advancement in the breeding program:\n\n")

for (i in seq_along(selected)) {
  mgidi_val <- result$MGIDI[result$MGIDI$Genotype == selected[i], "MGIDI"]
  cat(sprintf("  %d. %s (MGIDI = %.3f)\n", i, selected[i], mgidi_val))
}

cat("\nThese genotypes show the best overall performance across all\n")
cat("breeding objectives, balancing high yield, quality traits, and\n")
cat("agronomic characteristics.\n")

# Detailed performance of top genotype
# ------------------------------------
top_gen <- selected[1]
cat("\n", paste(rep("-", 50), collapse = ""), "\n")
cat("Detailed Performance of Top Genotype:", top_gen, "\n")
cat(paste(rep("-", 50), collapse = ""), "\n")

top_data <- wheat_data[wheat_data$GEN == top_gen, -1]
pop_means <- colMeans(wheat_data[, -1])
deviations <- as.numeric(top_data) - pop_means
pct_dev <- (deviations / abs(pop_means)) * 100

detail_df <- data.frame(
  Trait = names(top_data),
  Value = as.numeric(top_data),
  Population_Mean = pop_means,
  Deviation = deviations,
  Deviation_Pct = pct_dev,
  Direction = ifelse(ideotype == 1, "Maximize", "Minimize")
)
print(round(detail_df, 2))
