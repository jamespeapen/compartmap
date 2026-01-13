# Transform/normalize counts using TF-IDF

Transform/normalize counts using TF-IDF

## Usage

``` r
transformTFIDF(mat, scale.factor = 100000, count.min = 0, count.max = 1)
```

## Arguments

- mat:

  n x p input matrix (n = samples/cells; p = rna counts)

- scale.factor:

  Scaling factor for the term-frequency (TF)

- count.min:

  The minimum expression count used for TF-IDF. Binarizes when
  `count.min` = 0 and `count.max` = 1.

- count.max:

  The maximum expression count used for TF-IDF. Binarizes when
  `count.min` = 0 and `count.max` = 1. binarizes the matrix. A `cap`
  value greater than 1 will cap counts at that value.

## Value

A TF-IDF transformed matrix of the same dimensions as the input

## Details

This function and its helpers were modeled after or taken from:

- http://andrewjohnhill.com/images/posts/2019-5-6-dimensionality-reduction-for-scatac-data/analysis.html

- https://divingintogeneticsandgenomics.rbind.io/post/clustering-scatacseq-data-the-tf-idf-way/

## Examples

``` r
m <- 1000
n <- 100
mat <- round(matrix(runif(m * n), m, n))
# Input needs to be a tall matrix
tfidf <- transformTFIDF(mat)
```
