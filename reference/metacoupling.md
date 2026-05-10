# Metacoupling Analysis

Metacoupling Analysis

## Usage

``` r
metacoupling(
  data,
  swm_peri = NULL,
  swm_tele = NULL,
  weight = NULL,
  method = c("standard", "wang", "fan"),
  threads = 1
)
```

## Arguments

- data:

  A numeric matrix or data.frame. Rows are observations, columns are
  indicators.

- swm_peri:

  A numeric matrix representing the **peri (local) spatial weight
  matrix**. Must be square with dimension equal to `nrow(data)`. If
  `NULL`, a zero matrix is used.

- swm_tele:

  A numeric matrix representing the **tele (long-distance) spatial
  weight matrix**. Must be square with dimension equal to `nrow(data)`.
  If `NULL`, a zero matrix is used.

- weight:

  Numeric vector of indicator weights. Must have length equal to
  `ncol(data)`. If `NULL`, equal weights are used.

- method:

  Coupling model. One of `"standard"`, `"wang"`, or `"fan"`.

- threads:

  Number of threads used in computation.

## Value

A data.frame with:

- `Intra_C`: intra-system coupling degree

- `Intra_D`: intra-system coordination degree

- `Peri_C`: peri-coupling degree

- `Peri_D`: peri coordination degree

- `Tele_C`: tele-coupling degree

- `Tele_D`: tele coordination degree

## Details

Full model definitions and formulas are available at:
<https://github.com/stscl/coupling/discussions/8>

## Note

Input values should be normalized to `[0, 1]`. Spatial weight matrices
are typically symmetric.

## Examples

``` r
set.seed(42)
mat = matrix(runif(20), nrow = 5)
swm1 = apply(matrix(runif(25), 5, 5), 1, \(.x) .x / sum(.x))
swm2 = apply(matrix(runif(25), 5, 5), 1, \(.x) .x / sum(.x))
coupling::metacoupling(mat, swm1, swm2)
#>     Intra_C   Intra_D    Peri_C    Peri_D    Tele_C    Tele_D
#> 1 0.9497342 0.8199577 0.9090374 0.7621545 1.3099835 1.1089423
#> 2 0.9905121 0.9136491 1.1690717 1.0278882 0.5735267 0.5064461
#> 3 0.6926110 0.5050229 0.7869004 0.6461619 0.6749018 0.5429212
#> 4 0.9148146 0.7122036 0.8078109 0.6369757 1.2478349 0.9992342
#> 5 0.9877654 0.7649259 0.8117074 0.6539408 0.7006156 0.5491195
```
