MultiBD
======

`MultiBD` is an `R` package for direct likelihood-based inference of multivariate birth-death processes. 

## Installation

1. Install (if necessary) package dependencies and helpers:
```{r}
install.packages(c("Rcpp", "RcppParallel", "BH", "devtools"))
```

2. Install `MultiBD` from `github` (until package becomes available via `CRAN`):
```{r}
devtools::install_github("msuchard/MultiBD")
```

## Short example

```{r}
execute_short_example() # TODO
```


## Vignettes

1. **NAME**

## License
`MultiBD` is licensed under Apache License 2.0

## Development status

[![Build Status](https://travis-ci.org/msuchard/MultiBD.svg?branch=master)](https://travis-ci.org/msuchard/MultiBD)

Beta

## Acknowledgements
- This project is supported in part through the National Science Foundation grant DMS 1264153 and National Institutes of Health grant R01 AI107034.

## References

1. Ho LST, Xu J, Crawford FW, Minin VN, Suchard MA.
[Birth(death)/birth-death processes and their computable transition probabilities with statistical applications](http://arxiv.org).
*arXiv preprint arXiv*:0000.00000, 2016.
