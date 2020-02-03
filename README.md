# R codes and data accompanying the paper "The role of intrinsic dimension in high-resolution player tracking data - Insights in basketball"

__Authors:__

*Edgar Santos-Fernandez*, PhD. School of Mathematical Sciences. Y Block, Floor 8, Gardens Point Campus
Queensland University of Technology. GPO Box 2434. Brisbane, QLD 4001. Australia

*Francesco Denti*. Department of Economics, Università della Svizzera italiana, Lugano, Switzerland
Department of Statistics and Quantitative Methods, University of Milan - Bicocca, Milan, Italy


*Kerrie Mengersen*, Dist. Professor of Statistics. Deputy Director, ARC Centre for Mathematical and Statistical Frontiers;
Director, QUT Centre for Data Science. School of Mathematical Sciences. Y Block, Floor 8, Gardens Point Campus.
Queensland University of Technology. GPO Box 2434. Brisbane, QLD 4001. Australia


*Antonietta Mira*, Professor of Statistics. Director
Data Science Lab. Institute of Computational Science
Università della Svizzera italiana, USI, Lugano.
University of Insubria, Italy


# Installation:
Computations in the article are based on the R package ```Rhidalgo``` which stands for Heterogeneous Intrinsic
Dimension Algorithm.

```
library(devtools)
install_github("https://github.com/EdgarSantos-Fernandez/id_basketball/Rhidalgo")
```

```
library("Rhidalgo")
```

# Analysis of movement data:

## Data: 

* game.12.25.2015.CLE.at.GSW.RDS: Shot chart variables from the game GSW vs CLE played on 25-12-2015.
Locations of the 10 players in the court when the shot was taken are the columns x_ti and y_ti. Where x and y are the locations in the horizontal and the vertical axis respectively. t is the team (a = away) a (h = home). i is the player number (1-5).

* movement_data_play15.RDS: Player tracking data from event id = 15.  https://youtu.be/jb57MFQLoRo?t=17  

![Alt text](https://github.com/EdgarSantos-Fernandez/id_basketball/blob/master/event15-4.gif?raw=true "Title")

