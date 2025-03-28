## Rare-variant association studies: When are aggregation tests more powerful than single-variant tests? 

This repository contains code to perform analytic calculations and simulations from our paper. 

#### Analytic Calculations

Code to perform analytic calculations for single-variant (SV), burden, and SKAT tests can be found in `scripts/Analytic_Calculations.R`. Power for these three tests can be computed in the general case where we ONLY assume independent variants, and in the simple case where we further assume all variants have equal MAF, all causal variants have equal effect sizes, and all variants have equal weights in the mask. Also see our Shiny app: [Analytic Calculations App](https://debrajbose.shinyapps.io/analytic_calculations/)
