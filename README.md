
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lrspline

<!-- badges: start -->
<!-- badges: end -->

The R package lrspline provides functions for fitting smoothing splines
to massive data.

## Installation

You can install the *development* version of lrspline from
[Github](https://github.com/danqingxu/lrspline):

    # install.packages("devtools")
    devtools::install_github("danqingxu/lrspline")

<!-- ## Usage -->
<!-- ```{r example, eval = FALSE} -->
<!-- library(QDRS) -->
<!-- # EHR demo data -->
<!-- # sample.set: matrix of binary features. -->
<!-- # sample.group: a grouping factor with two levels "Case" and "Control".  -->
<!-- # training: an index vector that indicates training rows. -->
<!-- sample.set = EHR$sample.set -->
<!-- sample.group = EHR$sample.group -->
<!-- training = EHR$training -->
<!-- # Compute QDRSs -->
<!-- PheRS.res = PheRS(X = sample.set, group = sample.group) -->
<!-- Eigen.res = eigen.score(X = sample.set, training = training, scale = TRUE) -->
<!-- PC.res = PC(X = sample.set, group = sample.group, training = training, scale = TRUE, pc.num = 1:2) -->
<!-- LPC.res = LPC(X = sample.set, group = sample.group, training = training, scale = TRUE) -->
<!-- LVS.res = LVS.score(X = sample.set, Y = NULL, family = "binomial", starting.choice = "random", p.seed = 124) -->
<!-- NMF1.res = NMF1(X = sample.set) -->
<!-- # Assess the performance -->
<!-- pairwise.wilcox(x = LPC.res$scores, g = sample.group) -->
<!-- pairwise.auc(x = LPC.res$scores, g = sample.group) -->
<!-- pairwise.upr(x = LPC.res$scores, g = sample.group) -->
<!-- ``` -->

## License

This package is free and open source software, licensed under GPL
(&gt;=2).
