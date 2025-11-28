# MGIDI Tutorial: Multi-trait Selection in Plant Breeding

## Introduction

This tutorial demonstrates how to use the MGIDI (Multi-trait Genotype-Ideotype Distance Index) implementation for selecting superior genotypes in plant breeding programs.

## Getting Started

### Loading the Package

```r
library(metan)
```

### Preparing Your Data

MGIDI requires data in one of these formats:

1. **Data frame with genotype column**: First column is genotype names, remaining columns are traits
2. **Numeric matrix**: All columns are traits (genotype names will be auto-generated)

```r
# Example: Data frame with genotype column
data <- data.frame(
  GEN = paste0("G", 1:20),    # Genotype names
  GY = rnorm(20, 5, 1),       # Grain yield
  HI = rnorm(20, 0.4, 0.05),  # Harvest index
  PH = rnorm(20, 100, 10)     # Plant height
)

# View the data
head(data)
```

## Basic Usage

### Step 1: Define Breeding Objectives

Create a vector specifying the direction for each trait:
- `+1` or `"h"`: Maximize (higher values are better)
- `-1` or `"l"`: Minimize (lower values are better)

```r
# Numeric format
directions <- c(GY = 1, HI = 1, PH = -1)

# Or character format
directions <- c("h", "h", "l")
```

### Step 2: Run MGIDI Analysis

```r
result <- mgidi_index(
  data = data,
  ideotype_directions = directions,
  n_factors = NULL,    # Auto-determine number of factors
  verbose = TRUE       # Print summary output
)
```

### Step 3: View Results

```r
# MGIDI ranking (lower = better)
result$MGIDI

# Factor loadings
result$rotated_loadings

# Factor scores for genotypes
result$scores_gen

# Communalities
result$communalities
```

### Step 4: Select Top Genotypes

```r
# Select top 20% of genotypes
selected <- select_genotypes(result, intensity = 0.20)

# Or select specific number
selected <- select_genotypes(result, n = 5)

print(selected)
```

## Understanding the Output

### MGIDI Data Frame

The `result$MGIDI` data frame contains:
- `Genotype`: Genotype name
- `MGIDI`: MGIDI value (lower is better)

### Factor Analysis Results

- `result$PCA`: Principal component analysis summary
- `result$FA`: Factor analysis with loadings and communalities
- `result$KMO`: Kaiser-Meyer-Olkin measure (should be > 0.6)
- `result$MSA`: Measure of sampling adequacy per trait

### Factor Contributions

The `result$contri_fac` data frame shows the contribution of each factor to the MGIDI for each genotype. High contribution indicates weakness in that factor.

## Visualization

### MGIDI Ranking Plot

```r
# Bar plot of MGIDI values
p1 <- plot_mgidi_ranking(result, selection_intensity = 0.20)
print(p1)
```

### Factor Contributions Plot

```r
# Stacked bar chart of factor contributions
p2 <- plot_factor_contributions(result, selection_intensity = 0.20)
print(p2)
```

### Selection Gains Plot

```r
# Show genetic gains per trait
p3 <- plot_selection_gains(result, selection_intensity = 0.20)
print(p3)
```

### Biplot

```r
# PCA-style biplot of factor scores
p4 <- biplot_mgidi(result, selection_intensity = 0.20)
print(p4)
```

## Advanced Usage

### Using Weights

Assign different weights to traits to prioritize certain characteristics:

```r
# Higher weight for yield (GY)
result <- mgidi_index(
  data = data,
  ideotype_directions = c(1, 1, -1),
  weights = c(2, 1, 1)  # GY has double weight
)
```

### Specifying Number of Factors

Override the Kaiser criterion:

```r
# Force 2 factors
result <- mgidi_index(
  data = data,
  ideotype_directions = c(1, 1, -1),
  n_factors = 2
)
```

### Computing Selection Differential

```r
selected <- select_genotypes(result, intensity = 0.20)

sel_diff <- selection_differential(
  original_data = result$data,
  selected_genotypes = selected,
  trait_directions = result$directions
)

print(sel_diff)
```

### Getting Detailed Information on Selected Genotypes

```r
selected <- select_genotypes(result, intensity = 0.20)

details <- get_selected_details(result, selected)

# Trait values
details$trait_values

# Factor scores
details$factor_scores

# MGIDI values
details$mgidi_values
```

## Interpreting Results

### KMO Value
- \> 0.9: Excellent
- 0.8-0.9: Good
- 0.7-0.8: Acceptable
- 0.6-0.7: Mediocre
- < 0.6: Poor (factor analysis may not be appropriate)

### Communalities
- Values close to 1 indicate traits are well-represented by the factors
- Low values (< 0.5) suggest the trait may not fit the factor model well

### Factor Contributions
- Low contribution = strength (genotype is close to ideotype for that factor)
- High contribution = weakness (genotype is far from ideotype for that factor)

### Selection Differential
- Positive SD for "maximize" traits = good selection
- Negative SD for "minimize" traits = good selection
- The `Goal` column indicates whether selection achieved the desired direction

## Best Practices

1. **Standardize data**: MGIDI automatically rescales data, but ensure your input is clean
2. **Check KMO**: Values should be > 0.6 for reliable factor analysis
3. **Examine communalities**: Low values may indicate problematic traits
4. **Consider correlations**: Highly correlated traits will load on the same factor
5. **Use appropriate selection intensity**: Typically 10-30% depending on program needs

## Common Issues and Solutions

### Issue: KMO is too low
**Solution**: Consider removing traits with low MSA or increasing sample size

### Issue: Only one factor extracted
**Solution**: May indicate traits are highly correlated; consider reducing trait set

### Issue: MGIDI values are very similar
**Solution**: Genotypes may be genuinely similar; consider more diverse germplasm

### Issue: Selection differential is in wrong direction
**Solution**: Check that ideotype directions are correctly specified

## Complete Example

```r
library(metan)

# Create data
set.seed(123)
data <- data.frame(
  GEN = paste0("G", 1:30),
  Yield = rnorm(30, 5, 1),
  Quality = rnorm(30, 80, 10),
  Height = rnorm(30, 100, 15),
  Maturity = rnorm(30, 120, 10)
)

# Define breeding objectives
# Maximize: Yield, Quality
# Minimize: Height, Maturity
directions <- c(Yield = 1, Quality = 1, Height = -1, Maturity = -1)

# Run MGIDI
result <- mgidi_index(data, directions, verbose = TRUE)

# Select top 20%
selected <- select_genotypes(result, intensity = 0.20)
cat("Selected genotypes:", paste(selected, collapse = ", "), "\n")

# View selection gains
sel_diff <- selection_differential(result$data, selected, directions)
print(sel_diff)

# Create visualizations
p1 <- plot_mgidi_ranking(result)
p2 <- plot_selection_gains(result)

# Combine plots (requires patchwork)
print(p1)
print(p2)
```

## References

Olivoto, T., and Nardino, M. (2020). MGIDI: toward an effective multivariate selection in biological experiments. Bioinformatics, 37(10), 1383-1389.

## See Also

- `?mgidi_index` - Main MGIDI function
- `?select_genotypes` - Selection function
- `?plot_mgidi_ranking` - Visualization
- `?mgidi` - Original metan implementation
