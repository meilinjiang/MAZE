# MAZE
Mediation Analysis for ZEro-inflated mediators

## Installing the R package
```r
require(devtools)
devtools::install_github("https://github.com/meilinjiang/MAZE", build_vignettes = TRUE)
```

## Detailed vignettes
```r
browseVignettes("MAZE")
```

## Example
```r
library(MAZE)
# load the example dataset "zinb10"
data(zinb10)
# call MAZE() to perform mediation analysis
maze_out <- MAZE(data = zinb10,
                 distM = c('zilonm', 'zinbm', 'zipm'),  K = 1,
                 selection = 'AIC',
                 X = 'X', M = 'Mobs', Y = 'Y', Z = NULL,
                 XMint = c(TRUE, FALSE),
                 x1 = 0, x2 = 1, zval = NULL, mval = 0,
                 B = 20, seed = 1)
## results of selected mediation model
maze_out$results_effects # indirect and direct effects
maze_out$selected_model_name # distribution of the mediator and number of components K in the selected mediation model
maze_out$results_parameters # model parameters
maze_out$BIC; maze_out$AIC # BIC and AIC of the selected mediation model
```
