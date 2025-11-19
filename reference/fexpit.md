# Helper function: expanded expit

Helper function: expanded expit

## Usage

``` r
fexpit(x, sqz = 0.000001)
```

## Arguments

- x:

  a vector of values between -Inf and +Inf

- sqz:

  the amount by which we 'squoze', default is .000001

## Value

       a vector of values between 0 and 1 inclusive

## Examples

``` r
x <- rnorm(n = 1000)
summary(x)
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> -3.40416 -0.71838 -0.03249 -0.02624  0.68695  3.40723 

sqz <- 1 / (10**6)
p <- fexpit(x, sqz = sqz)
summary(p)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.03217 0.32775 0.49188 0.49502 0.66529 0.96793 

all((abs(x - flogit(p)) / x) < sqz)
#> [1] TRUE
all(abs(x - flogit(fexpit(x))) < sqz)
#> [1] TRUE
```
