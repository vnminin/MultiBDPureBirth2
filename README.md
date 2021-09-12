MultiBDPureBirth2
======

`MultiBDPureBirth2` is a slight extention of the `MultiBD`  `R` package for direct likelihood-based inference of multivariate birth-death processes.  

## Installation

```
1. Install  `MultiBDPureBirth2` from `github`:
```{r}
devtools::install_github("vnminin/MultiBDPureBirth2")
```

## Short example

```{r}
library(MultiBDPureBirth2)

power_list = list(0.9,1,1)
names(power) = c("powS", "powI_inf",  "powI_rem")

SIR_prob_pure_birth(0.5, 0.1, 0.2, S0=1000, I0=10, nSI=4, nIR=0, direction = c("Forward","Backward"), power = power, nblocks = 20, tol = 1e-12, computeMode = 0, nThreads = 1)
```


## Vignettes

1. [Simple MCMC under SIR](https://github.com/msuchard/MultiBD/blob/master/inst/doc/SIR-MCMC.pdf)
2. [SIR model and proposed branching approximation](https://github.com/msuchard/MultiBD/blob/master/inst/doc/SIRtrans.pdf)

## License
`MultiBDPureBirth2` is licensed under Apache License 2.0

## Acknowledgements
- This project is supported in part through the National Science Foundation grant DMS 1264153 and National Institutes of Health grant R01 AI107034.

## References

1. Ho LST, Xu J, Crawford FW, Minin VN, Suchard MA (2018).
[Birth/birth-death processes and their computable transition probabilities with biological applications](https://link.springer.com/article/10.1007/s00285-017-1160-3).
Journal of Mathematical Biology 76(4) 911-944.
2. Ho LST, Crawford FW, Suchard MA (2018).
[Direct likelihood-based inference for discretely observed stochastic compartmental models of infectious disease](https://arxiv.org/abs/1608.06769).
Annals of Applied Statistics. In press.
