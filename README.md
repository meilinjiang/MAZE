# MAZE
Mediation Analysis for ZEro-inflated mediators

## Installing the R package
require(devtools)

devtools::install_github("https://github.com/meilinjiang/MAZE", 
                         build_vignettes = TRUE)

## Detailed vignettes
browseVignettes("MAZE")


## Example
library(MAZE) 

data(zinb10) # load the example dataset "zinb10"

maze_out <- MAZE(data=zinb10, 
                 distM=c('zilonm', 'zinbm', 'zipm'), 
                 K = 1,
                 selection = "AIC",
                 X='X', M='Mobs', Y='Y', Z=NULL, 
                 XMint = c(TRUE, FALSE),
                 x1=0, x2=1, zval = NULL, mval = 0,
                 B=20, seed=1)

maze_out$results_effects # indirect and direct effects

maze_out$selected_disM # distribution of the mediator in the selected mediation model
