# MGIDI Algorithm Documentation

## Multi-trait Genotype-Ideotype Distance Index

### Introduction

The MGIDI (Multi-trait Genotype-Ideotype Distance Index) is a multivariate selection index designed for plant breeding programs. It was proposed by Olivoto and Nardino (2020) as a tool for selecting superior genotypes based on multiple traits simultaneously.

### Mathematical Foundation

#### MGIDI Formula

The MGIDI for genotype *i* is calculated as the Euclidean distance from the genotype's factor scores to the ideotype's factor scores:

$$MGIDI_i = \sqrt{\sum_{j=1}^{f}(F_{ij} - F_j^*)^2}$$

Where:
- $MGIDI_i$ = Multi-trait genotype-ideotype distance index for genotype *i*
- $F_{ij}$ = Factor score of genotype *i* on factor *j*
- $F_j^*$ = Ideal score for factor *j* (ideotype)
- $f$ = Number of factors retained

The genotype with the lowest MGIDI is closest to the ideotype and should be preferred.

### Algorithm Steps

#### Step 1: Data Preparation and Rescaling

The original trait values are rescaled to a 0-100 scale based on the breeding objectives:

For traits to **maximize**:
$$X_{rescaled} = 100 \times \frac{X - X_{min}}{X_{max} - X_{min}}$$

For traits to **minimize**:
$$X_{rescaled} = 100 \times \frac{X_{max} - X}{X_{max} - X_{min}}$$

#### Step 2: Correlation Matrix

Compute the Pearson correlation matrix among rescaled traits:

$$R = \frac{1}{n-1}(Z - \bar{Z})'(Z - \bar{Z})$$

#### Step 3: Eigendecomposition

Perform eigendecomposition of the correlation matrix:

$$R = V \Lambda V'$$

Where:
- $V$ = Matrix of eigenvectors
- $\Lambda$ = Diagonal matrix of eigenvalues

#### Step 4: Factor Extraction (Kaiser Criterion)

Retain factors with eigenvalues ≥ 1 (or custom threshold):

$$k = \sum_{j=1}^{p} I(\lambda_j \geq 1)$$

#### Step 5: Initial Loadings

Compute initial factor loadings:

$$L_{initial} = V_k \sqrt{\Lambda_k}$$

Where $V_k$ and $\Lambda_k$ are the retained eigenvectors and eigenvalues.

#### Step 6: Varimax Rotation

Apply Varimax rotation to achieve a simpler factor structure:

$$L_{rotated} = L_{initial} \times T$$

Where $T$ is the orthogonal rotation matrix that maximizes the variance of squared loadings.

#### Step 7: Factor Scores

Compute factor scores using the regression method:

$$F = Z \times R^{-1} \times L_{rotated}$$

#### Step 8: Ideotype Construction

For each factor *j*, the ideotype score is determined by the trait directions:

$$F_j^* = \begin{cases} 
\max(F_j) & \text{if } \sum_k (L_{kj} \times d_k) > 0 \\
\min(F_j) & \text{otherwise}
\end{cases}$$

Where $d_k$ is the direction for trait *k* (+1 for maximize, -1 for minimize).

#### Step 9: MGIDI Calculation

Compute the Euclidean distance from each genotype to the ideotype.

#### Step 10: Ranking and Selection

Rank genotypes by MGIDI (ascending - lower is better) and select the top percentage.

### Factor Analysis in MGIDI

Factor analysis reduces the dimensionality of the trait data while preserving the essential information:

1. **Communality**: The proportion of variance in each trait explained by the factors
   $$h^2_k = \sum_{j=1}^{f} L_{kj}^2$$

2. **Uniqueness**: The proportion of variance not explained
   $$u^2_k = 1 - h^2_k$$

3. **KMO Test**: Kaiser-Meyer-Olkin measure of sampling adequacy
   $$KMO = \frac{\sum \sum_{i \neq j} r_{ij}^2}{\sum \sum_{i \neq j} r_{ij}^2 + \sum \sum_{i \neq j} q_{ij}^2}$$

   Where $r_{ij}$ are correlations and $q_{ij}$ are partial correlations.

### Selection Differential

The selection differential measures the difference between the population mean and the selected mean:

$$SD = \bar{X}_s - \bar{X}_o$$

$$SD\% = \frac{\bar{X}_s - \bar{X}_o}{|\bar{X}_o|} \times 100$$

If heritability ($h^2$) is known, the expected response to selection is:

$$R = h^2 \times SD$$

### Advantages of MGIDI

1. **Multivariate Approach**: Considers all traits simultaneously
2. **Accounts for Correlations**: Uses factor analysis to handle correlated traits
3. **Flexible Ideotype**: Different directions can be specified for each trait
4. **Weighted Selection**: Traits can be weighted according to importance
5. **Interpretable**: Factor contributions show strengths and weaknesses

### Interpretation Guidelines

- **MGIDI Values**: Lower values indicate genotypes closer to the ideotype
- **Factor Contributions**: High contribution indicates weakness in that factor
- **Selection Intensity**: Typically 10-30% of genotypes are selected
- **KMO**: Values > 0.6 indicate adequate sampling for factor analysis

### References

Olivoto, T., and Nardino, M. (2020). MGIDI: toward an effective multivariate selection in biological experiments. Bioinformatics, 37(10), 1383-1389. https://doi.org/10.1093/bioinformatics/btaa981

### Implementation Notes

The pure R implementation in this package:

1. Uses standard R matrix operations for core computations
2. Implements Varimax rotation without external packages
3. Provides detailed intermediate results for educational purposes
4. Is fully compatible with the original metan package implementation
