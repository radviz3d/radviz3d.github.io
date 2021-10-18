## ---- echo=FALSE-----------------------------------------------------------------------
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, error = FALSE)


## --------------------------------------------------------------------------------------
library(MixSim)
library(RColorBrewer)
library(MASS)
library(freqparcoord)
library(expm)
library(grDevices)
library(radviz3d)
library(rgl)



## --------------------------------------------------------------------------------------
set.seed(2020)
res_radviz <- list()
sim.palette <- brewer.pal(8,"Dark2")[c(5,8,4,7,3)]
# load dataset:sim_data from the radviz3d package
for (i in 1:3){
  df <- sim_data[[i]]
  radialvis3d(data = df[,-1], cl = factor(df$class), domrp = F, doGtrans = F, lwd = 2, 
            alpha = 0.025, point.cex = 0.15, color = sim.palette, class.labels = NULL)
  rgl::rgl.viewpoint(zoom = 0.6)
  res_radviz[[i]] <- rgl::rglwidget()
}



## --------------------------------------------------------------------------------------
res_radviz[[1]]


## --------------------------------------------------------------------------------------
res_radviz[[2]]


## --------------------------------------------------------------------------------------
res_radviz[[3]]

