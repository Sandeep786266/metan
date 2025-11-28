# Example: Simple MGIDI Demonstration
# ====================================
# A simple example to demonstrate the basic usage of the MGIDI algorithm.
#
# This example uses the pure R implementation of MGIDI (mgidi_index function)
# which provides educational, transparent code without external dependencies
# for core statistical computations.

library(metan)

# Create simple example data
# --------------------------
set.seed(42)

# 10 genotypes with 3 traits
simple_data <- data.frame(
  GEN = paste0("G", 1:10),
  Yield = c(5.2, 6.1, 4.8, 7.2, 5.5, 6.8, 4.5, 5.9, 6.3, 5.0),
  Quality = c(8.5, 7.2, 9.1, 6.8, 8.0, 7.5, 9.5, 7.8, 7.0, 8.8),
  Height = c(95, 110, 85, 120, 100, 105, 80, 98, 115, 90)
)

cat("Example Data:\n")
print(simple_data)

# Define ideotype
# ---------------
# We want to:
#   - MAXIMIZE Yield (higher is better)
#   - MAXIMIZE Quality (higher is better)  
#   - MINIMIZE Height (shorter plants are better)

directions <- c(Yield = 1, Quality = 1, Height = -1)

cat("\nBreeding Objectives:\n")
cat("  Yield: Maximize\n")
cat("  Quality: Maximize\n")
cat("  Height: Minimize\n")

# Run MGIDI
# ---------
cat("\n--- Running MGIDI Analysis ---\n")

result <- mgidi_index(
  data = simple_data,
  ideotype_directions = directions,
  verbose = TRUE
)

# Display full MGIDI ranking
cat("\nComplete MGIDI Ranking:\n")
print(result$MGIDI)

# Select top 30% (3 genotypes)
selected <- select_genotypes(result, intensity = 0.30)
cat("\nSelected genotypes (top 30%):", paste(selected, collapse = ", "), "\n")

# Show trait values for selected genotypes
cat("\nTrait values for selected genotypes:\n")
print(simple_data[simple_data$GEN %in% selected, ])

# Compare with simple ranking
# ---------------------------
cat("\n--- Comparison with Simple Trait Ranking ---\n")

# Standardize traits according to direction
std_data <- simple_data
std_data$Yield_std <- (std_data$Yield - min(std_data$Yield)) / 
                       (max(std_data$Yield) - min(std_data$Yield))
std_data$Quality_std <- (std_data$Quality - min(std_data$Quality)) / 
                         (max(std_data$Quality) - min(std_data$Quality))
std_data$Height_std <- (max(std_data$Height) - std_data$Height) / 
                        (max(std_data$Height) - min(std_data$Height))

# Simple sum of standardized values
std_data$Simple_Index <- std_data$Yield_std + std_data$Quality_std + std_data$Height_std
std_data <- std_data[order(-std_data$Simple_Index), ]

cat("\nSimple Index Ranking:\n")
print(std_data[, c("GEN", "Yield_std", "Quality_std", "Height_std", "Simple_Index")])

# Compare selections
simple_top <- std_data$GEN[1:3]
cat("\nComparison of selections:\n")
cat("  MGIDI selected:", paste(selected, collapse = ", "), "\n")
cat("  Simple index selected:", paste(simple_top, collapse = ", "), "\n")

# Interpretation
# --------------
cat("\n--- Interpretation ---\n")
cat("
MGIDI provides a more sophisticated selection approach because it:

1. Uses factor analysis to identify underlying trait patterns
2. Accounts for correlations between traits
3. Creates an ideotype based on factor scores, not raw values
4. Considers the multivariate nature of the breeding objectives

The simple sum index treats all traits as independent, which may
lead to suboptimal selection when traits are correlated.

Factor Analysis Summary:
")
print(result$FA)

cat("\nCommunalities (proportion of trait variance explained by factors):\n")
print(round(result$communalities, 3))
