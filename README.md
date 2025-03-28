## Rare-variant association studies: When are aggregation tests more powerful than single-variant tests? 

This repository contains code to perform analytic calculations and simulations from our paper. 

### Analytic Calculations

Code to perform analytic calculations for single-variant (SV), burden, and SKAT tests can be found in `scripts/Analytic_Calculations.R`. Power for these three tests can be computed in the general case where we ONLY assume independent variants, and in the simple case where we further assume all variants have equal MAF, all causal variants have equal effect sizes, and all variants have equal weights in the mask. Also see our [R Shiny app](https://debrajbose.shinyapps.io/analytic_calculations/) for a more user-friendly way to explore how power for these tests are affected by the different parameters.

### Simulations

Code to perform simulations for (SV), burden, SKAT, and SKAT-O tests can be found in `scripts/Simulations.R`. Genotypes of 100,000 individuals for 100 variants in a dummy gene is provided in `data/dummy_gene.txt` along with their annotations `data/annotations.txt` (the three columns are indicators of whether a variant is a PTV, deleterious missense, or other missense variant respectively). Please remember to change the locations of these files while running `scripts/Simulations.R`. Results obtained by running the script with default parameters are provided in the `results` folder (the three files correspond to the three masks used in our analyses). 
