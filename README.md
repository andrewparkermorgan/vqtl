vqtl
====

This package implements several methods for mapping variance-controlling loci (vQTL) in experimental crosses, whose input and output are compatible with Karl Broman's `R/qtl` (http://www.rqtl.org/).  Also included are functions for producing (static) QTL charts with `ggplot2`.

Methods implemented
----
* Double generalized linear model (DGLM) cf. Ronnegard & Valdar 2011 (http://dx.doi.org/10.1534/genetics.111.127068)

TODO: methods to add
----
* Ridge regression with marker heteroscedasticity (heteroscedastic effects model, HEM) cf. Shen _et al._ 2013 (http://dx.doi.org/10.1534/genetics.112.146720)
* Ridge regression with linear predictor for variance components (hierarchical generalized linear model, HGLM) using Ronnegard's `hglm` package (http://cran.r-project.org/web/packages/hglm/)