estimate.pars <- function(x, cl) {
    p <- ncol(x) 
    id <- as.integer(cl)
    K <- max(id)
    ## estimate mixture parameters for calculation of overlap
    Pi <- prop.table(tabulate(id))
    Mu <- t(sapply(1:K, function(k){ colMeans(x[id == k,]) }))
    S <- sapply(1:K, function(k){ var(x[id == k, ]) })
    dim(S) <- c(p, p, K)
    list(Pi = Pi, Mu = Mu, S = S)
}

prmtd.pars <- function(idx, parlist, cl) {
    Mu <- parlist$Mu[,idx]
    S <- parlist$S[idx, idx, ]
    list(Mu = Mu, S = S) 
}

pars.2d <- function(parlist, projmat) {
    Mu <- parlist$Mu %*% projmat
    S <- array(dim = c(2, 2, dim(parlist$S)[3]))
    for (i in 1:(dim(S)[3]))
        S[,,i] <- t(projmat) %*% parlist$S[,,i] %*% projmat
    list(Mu = Mu, S = S)
}

directionless_circular_permutations <- function( n ) {
    library(gtools)
    n1 <- n - 1L
    v <- seq.int( n1 )
    ix <- combinations( n1, 2L )
    jx <- permutations( n-3L, n-3L )
    jxrows <- nrow( jx )
    jxoffsets <- seq.int( jxrows )
    result <- matrix( n, nrow = factorial( n1 )/2L, ncol = n )
    k <- seq( 2L, n - 2L )
    for ( i in seq.int( nrow( ix ) ) ) {
        result[ ( i - 1L ) * jxrows + jxoffsets, k ] <-
            matrix( v[ -ix[ i, ] ][ jx ], nrow = jxrows )
    }
    result[ , 1L ] <- rep( ix[ , 1L ], each = jxrows )
    result[ , n1 ] <- rep( ix[ , 2L ], each = jxrows )
    result
}

frobenius.distance <- function(mat1, mat2) (sqrt(sum(diag((mat1 - mat2) %*% t(mat1-mat2)))))    

optimal_2d_axis_order <- function(x, cl, proj.mat) {
    library(MixSim)
    ll <- estimate.pars(x, cl)
    target.overlap <- overlap(Pi = ll$Pi, Mu = ll$Mu, S = ll$S)
    
    ## get all directionless permutations of coordinates
    ##    #print(target.overlap)
    p <- ncol(x)
    list.idx <- directionless_circular_permutations(p)
    curr.idx <- NULL
    curr.mse <- Inf

    for (i in 1:nrow(list.idx)) {
        idx <- list.idx[i,]
        ll.idx <- prmtd.pars(idx = idx, parlist = ll, cl = cl)
        ll.2d <- pars.2d(parlist = ll.idx, projmat = proj.mat) 
        idx.overlap <- overlap(Pi = ll$Pi, Mu = ll.2d$Mu, S = ll.2d$S)
        ##        #print(idx.overlap)
        ## idx.mse <- max((idx.overlap$OmegaMap - target.overlap$OmegaMap)^2)

        idx.mse <- frobenius.distance(idx.overlap$OmegaMap, target.overlap$OmegaMap)
        ##print(i)
        #print(idx)
        #print(idx.mse)
        if (idx.mse < curr.mse) {
            curr.idx <- idx
            curr.mse <- idx.mse
        }
        #print(curr.idx)
        #print(curr.mse)
    }
    #print(curr.idx)
    curr.idx
}


draw_circle <- function(nlines = 1000, radius = 1, ...) (lines(radius * cbind(cos((0:nlines)*2*pi/1000), sin((0:nlines)*2*pi/nlines)),...))


farthest.point.from.k.centers <- function(x, ctr) {

    x.arr <- array(x, dim = c(dim(x), nrow(ctr)))
    ctr.arr <- array(rep(ctr, each = nrow(x)), dim = dim(x.arr))

    each.dist <- apply(X = (x.arr - ctr.arr)^2, MARGIN = c(1, 3), FUN = sum)

    ##
    ## find closest distance of each point to a center
    ## 

    x.min <- apply(X = each.dist, MARGIN = 1, FUN = min)

    ##
    ## find point that is the fathest from any center
    ##

    which.max(x.min)
}


separated.class.points <- function(x, cl)
{
    ##
    ## find most separated points that belong to the different groups
    ##

    tmp.arr <- apply(X = x, MARGIN = 1, FUN = function(x) (sum(x^2)))

    id <- which.max(tmp.arr)

    x.arr <- x
    inc.arr <- matrix(x[id, ], ncol = ncol(x))

    cls <- cl[id]
    tmp.cl <- cl

    for (i in 1:length(unique(cl))) {
        x.arr <- x.arr[tmp.cl != tmp.cl[id],]
        tmp.cl <- tmp.cl[tmp.cl != tmp.cl[id]]
        id <- farthest.point.from.k.centers(x = x.arr, ctr = inc.arr) 
        cls <- c(cls, tmp.cl[id])
        inc.arr <- rbind(inc.arr, x.arr[id, ])
    }
    list(cls, inc.arr)
}

##
## x = data matrix (of n rows and p columns)
## angles = projections of each axis on to the circle (uniformly distributed along the circle by default)
## modify = TRUE if points be plotted by square root, FALSE by default
## cl = class indicator (no class by default)
## palette = the desired palette to be used, colorblind friendly by default.
## coord.labels = labels of the coordinate axis
## coord.font = font for coordinate axes labels, bold face by default
## coord.cex = character size of coordinate labels
## class.labels = class labels in case labels are present, no effect if cl is NULL 
## class.labels.locations = locations where class labels should be put, default is to try and put them close to class points that are as far are possible
## class.cex = character size of class labels
## font.family = font family to be used for the labels
## 
## Written: Ranjan Maitra, Ames, Iowa 50011
## Last modified: 20018/03/31
##

modradviz2d.fn <- function(x, angles = 2*pi*(0:(ncol(x)-1))/ncol(x), modify = FALSE, cl = NULL, pch = 16, ret.proj = F,
                           palette = rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")), coord.labels = colnames(x), coord.font = 2, coord.cex = 1.1, class.cex = 1.1,  class.labels = levels(factor(cl)), class.labels.locations = NULL, ...) {
    ##
    ## scale the dataset to be between 0 and 1 in each coordinate
    ##
#    y <- apply(X = x, MARGIN = 2, FUN = function(z)((z - min(z))/(max(z) - min(z))))
    y <- x
    if (modify)
        y <- apply(X = y, MARGIN = 2, FUN = function(x)(sqrt(x)/max(sqrt(x))))
    ##
    ## build projection matrix
    ##

    pr <- cbind(cos(angles), sin(angles))

    xproj <- y %*% pr

    max.rad <- sqrt(max(rowSums(xproj^2)))
    
    par(mar = rep(0,4))
    plot(c(0,0), type = "n", xlim = c(-1.2, 1.2)*max.rad, ylim = c(-1.2, 1.2)*max.rad, xlab = "", ylab = "", axes = FALSE)

    ##
    ## draw the circle
    ## 
    
    draw_circle(col="gray40", radius = max.rad)

    ##
    ## draw the major guidelines
    ##

    for (i in 1:3) {
        if (modify) {
            draw_circle(radius = max.rad*sqrt(i)/2, col = "gray20", lwd = 0.25)
        } else
            draw_circle(radius = max.rad*i/4, col = "gray20", lwd = 0.25)
    }

    ##
    ## draw the minor guidelines
    ##

   
    for (i in 1:12) {
        if (modify) {
            draw_circle(radius = max.rad*sqrt(i)/sqrt(12), col = "gray60", lwd = 0.25)
        } else
            draw_circle(radius = max.rad*i/12, col = "gray60", lwd = 0.25)
    }


    axis.proj <- diag(ncol(x)) %*% pr *max.rad

    for (i in 1:ncol(x))
        lines(x = c(0,axis.proj[i,1]), y = c(0,axis.proj[i,2]), col = 'gray40')

    text(x = axis.proj, labels = coord.labels, font = coord.font, cex = coord.cex, ...)
    if (is.null(cl)) {
        points(x = xproj, col = palette[1], pch = 20, fg = "gray80")
    } else {
        points(x = xproj, col = palette[as.numeric(factor(cl))], pch = pch, fg = "gray80")

        if (is.null(class.labels.locations)) {
            ll <- separated.class.points(x = xproj, cl = as.numeric(factor(cl)))
            text(ll[[2]] + c(-0.01, -0.01)*max.rad, col = palette[ll[[1]]], font = 4, cex = class.cex,  labels = class.labels[ll[[1]]], ...)
        }
        else
            text(class.labels.locations, col = palette, font = 4, cex = class.cex, labels = class.labels, ...)
    }
    if (ret.proj)
    xproj
}


##
## use with labels included
##

## modradviz2d(x = iris[,1:4], cl = iris[,5], modify = F, coord.labels = c("Sepal Length", "Sepal Width", "Petal Length", "Petal Length"), class.labels.locations = matrix(c(0.4, 0.6, -0.42, -0.1, 0.1, -0.67), ncol = 2, byrow = T))



##
## use without labels
##

## modradviz2d(x = iris[,1:4], modify = F, coord.labels = c("Sepal Length", "Sepal Width", "Petal Length", "Petal Length"))

##
## used with the 'xkcd' font family
##

## library(xkcd)


## modradviz2d(x = iris[,1:4], cl = iris[,5], modify = F, coord.labels = c("Sepal Length", "Sepal Width", "Petal Length", "Petal Length"), class.labels.locations = matrix(c(0.45, 0.6, -0.42, -0.1, 0.12, -0.67), ncol = 2, byrow = T), font.family = 'xkcd')

modradviz2d <- function(x, angles = 2*pi*(0:(ncol(x)-1))/ncol(x), modify = FALSE, cl = NULL, palette = rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")), coord.labels = colnames(x), coord.font = 2, coord.cex = 1.1, 
                        class.labels = levels(factor(cl)), class.labels.locations = NULL, opt.axis.order = TRUE, ret.proj = F, ...) {
  ##
  ## change all axis to be between 0 and 1
  ##
  
  y <- apply(X = x, MARGIN = 2, FUN = function(z)((z - min(z))/(max(z) - min(z))))
  x <- apply(X = y, MARGIN = 1, FUN = function(z)(if (sum(z) == 0)
    z else
      z/sum(z)))
  x <- t(x)
  
  idx.opt <- 1:ncol(x)
  if (!is.null(cl) & (opt.axis.order)) 
  idx.opt <- optimal_2d_axis_order(x = x, cl = cl, proj.mat =  cbind(cos(angles), sin(angles)))
  modradviz2d.fn(x = x[,idx.opt], angles = angles, modify = modify, cl = cl, 
                 palette = palette, coord.labels = coord.labels[idx.opt], 
                 coord.font = coord.font, coord.cex = coord.cex, 
                 class.labels = class.labels, ret.proj = ret.proj,
                 class.labels.locations = class.labels.locations, ...) -> xproj
  if (ret.proj) 
  xproj
}

#modradviz2d(x = iris[,1:4], cl = iris[,5], modify = F, coord.labels = c("Sepal Length", "Sepal Width", "Petal Length", "Petal Length"), class.labels.locations = matrix(c(0.4, 0.6, -0.42, -0.1, 0.1, -0.67), ncol = 2, byrow = T))

                                        # modradviz2d(x = iris[,1:4], cl = iris[,5], modify = F, coord.labels = c("Sepal Length", "Sepal Width", "Petal Length", "Petal Width"))

if (FALSE) {
    library(MixSim)
    
    Q <- MixGOM(goMega = 0.01, K = 4, p = 9)
    
    A <- simdataset(n = 500, Pi = Q$Pi, Mu = Q$Mu, S = Q$S)
    data = A$X
    colnames(data) = paste("x",1:9,sep='')
    modradviz2d(x = A$X, cl = factor(A$id), class.labels = paste("Group ", rep(1:5)), palette = RColorBrewer::brewer.pal(n = 8, name = "Dark2"))
}


#load("/home/maitra/WRITES/ericksad/General\ shaped/GRB-Z.rda")


#modradviz2d(x = iris[,1:4], cl = iris[,5], modify = F, coord.labels = c("Sepal Length", "Sepal Width", "Petal Length", "Petal Length"), class.labels.locations = matrix(c(0.45, 0.6, -0.42, -0.1, 0.12, -0.67), ncol = 2, byrow = T))



#modradviz2d(x = Z[,-1], cl = Z[,1], modify = F,
#            coord.labels=c(expression(T[50]),  expression(T[90]), expression(F[1]),   expression(F[2]),   expression(F[3]), 
# expression(F[4]),   expression(P[64]),  expression(P[256]),expression(P[1024])),            coord.cex = 2.25, class.cex = 2.5,
#            palette = RColorBrewer::brewer.pal(n = 8, name = "Dark2"),
#            class.labels.locations = matrix(
#            c( -0.07601210, 0.1085,  0.04003691, 0.08065406,  0.06672818, 
#              0.27016917, 0.202, 0.10471225, -0.17529177, -0.02834751),
#            ncol = 2, byrow = F)
#            )

#dev.copy2pdf(file = "GRB-radviz2d.pdf", height = 12, width = 12)
