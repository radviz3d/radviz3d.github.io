stars = function(x,cl,labels=NULL,col.palette){
  x.range <- apply(x, 2, range)
  z <- t((t(x) - x.range[1,]) / (x.range[2,] - x.range[1,]))
  d <- dim(z)[2] 
  prj <- t(sapply((1:d)/d, function(i) c(cos(2*pi*i), sin(2*pi*i))))
  star <- z %*% prj
  plot(rbind(apply(star, 2, range), apply(prj*1.25, 2, range)), 
       type="n", bty="n", xaxt="n", yaxt="n",
       main="", xlab="", ylab="")
  tmp <- apply(prj, 1, function(v) lines(rbind(c(0,0), v)))
  if(is.null(labels)) labels=paste0("X",1:ncol(z))
  text(prj * 1.1, labels=labels, cex=0.8, col="black")
  for (i in 1:nlevels(cl)) {
    points(star[which(cl==levels(cl)[i]),], pch=20, col=col.palette[i]);
  }
}