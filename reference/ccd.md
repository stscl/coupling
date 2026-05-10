# Coupling Coordination Degree (CCD)

Coupling Coordination Degree (CCD)

## Usage

``` r
ccd(data, weight = NULL, method = c("standard", "wang", "fan"), threads = 1)
```

## Arguments

- data:

  A numeric matrix or data.frame. Rows are observations, columns are
  indicators.

- weight:

  Numeric vector of indicator weights. Must have length equal to
  `ncol(data)`. If `NULL`, equal weights are used.

- method:

  Coupling model. One of `"standard"`, `"wang"`, or `"fan"`.

- threads:

  Number of threads used in computation.

## Value

A data.frame with:

- `C`: coupling degree

- `D`: coordination degree

## Details

Full model definitions and formulas are available at:
<https://github.com/stscl/coupling/discussions/3>

## Note

Input values should be normalized to `[0, 1]`.

## Examples

``` r
set.seed(42)
mat = matrix(runif(20), nrow = 5)
coupling::ccd(mat)
#>           C         D
#> 1 0.9497342 0.8199577
#> 2 0.9905121 0.9136491
#> 3 0.6926110 0.5050229
#> 4 0.9148146 0.7122036
#> 5 0.9877654 0.7649259
```
