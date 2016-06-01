###############################################################################
# (Functions)
# TODO: Add comment
# 
# Author: yenlow
# 22-Feb-16 refactoring, Nathanael Romano
###############################################################################

plotOR <- function(data, filename=NULL, title=NULL, 
                   legendgrps=c("ORmatched", "ORadj"), width=7, height=11,
                   bg="transparent", legend=TRUE) {
	names <- rownames(data)
	names <- c(names, "")
	col <- c("black","darkgray")
	pch <- c(18,18)
  cex <- c(1.2,1.2)
	lty <- c(1,1)
  xmin <- 0.9
	
	or.grp1 <- c(data$ORmatched, NA)  
	lower.grp1 <- c(data$ORlow_matched, NA)
	upper.grp1 <- c(data$ORupp_matched, NA)
	
	or.grp2 <- c(data$ORadj, NA)
	lower.grp2 <- c(data$ORlow_adj, NA)
	upper.grp2 <- c(data$ORupp_adj, NA)
	
	# start plotting
	if (!is.null(filename)){
	  png(filename, width=width, height=height, units="in", bg=bg, res=300)
	}
	par(mar=c(2, 5, 3, 1))
	
	y.axis <- (length(legendgrps)+1)*seq(length(or.grp1),1)
	y.axis.adj <- y.axis-1
	y.axis.adjspace <- y.axis-2
	x.lim=c(max(min(floor(10*min(lower.grp1, lower.grp2, na.rm=T))/10,xmin),0.01),
	        ceiling(10*max(upper.grp1, upper.grp2, na.rm=T))/10)

	# grp1
	plot(or.grp1, y.axis, xlim=x.lim,axes=F, type = "p",ylab="",
       pch=pch[1], col=col[1], main="", log="x",cex=cex[1]) 
	segments(lower.grp1, y.axis, upper.grp1, y.axis, lty=lty[1],lwd=1.5, 
	         col=col[1])
	
	# grp2
	points(or.grp2, y.axis.adj, pch=pch[2], col=col[2],cex=cex[2])
	segments(lower.grp2, y.axis.adj, upper.grp2, y.axis.adj, lty=lty[2],
	         lwd=1.5, col=col[2])
	# coef +/-1.96*se = 95% interval, lwd adjusts line thickness
	
	# plot gridlines
  if (x.lim[2]-x.lim[1]<2) {
    xticks=round(seq(x.lim[1],x.lim[2],by=0.1),1)
  } else if(x.lim[2]-x.lim[1]<5) {
    xticks=round(seq(x.lim[1],x.lim[2],by=0.5),1)
    if(x.lim[1]==xmin) xticks=c(xmin,xticks)
  } else {
    xticks=round(seq(x.lim[1],x.lim[2],by=1))
    if(x.lim[1]==xmin) xticks=c(xmin,xticks)
  }
	abline(v=xticks, lty=1, col="gray", lwd=0.3)  # plot vertical gridlines
	abline(v=1, lty=1, col="black", lwd=1.5) # plot OR=1 reference line
	
	# plot axes
	axis(1, at=c(xticks,xmin,1), labels =c(xticks,xmin,1), tick = T,
	     # draw x-axis and labels with tick marks
			cex.axis = 0.8, mgp = c(2,.7,0))
	    # reduce label size, moves labels closer to tick marks    
	axis(2, at = y.axis, label = names, las = 2, tick = F, 
	     mgp = c(2,.4,0), cex.axis = 0.8)
	# draw y-axis with tick marks, make labels perpendicular 
	# to axis and closer to axis
	
	#add title
	title(title, cex.main=1)
	if(legend==TRUE) legend("topleft", legend=legendgrps,
	                        pch=pch, col=col, lty=lty, lwd=1.5, cex=0.8)
	# legend(0, 1.3, legend=c("no matching", "psm matching"), fill=c(col1, col2))
	# legend(2, length(y.axis)-1, legend=c("no matching", "psm matching"), 
	# fill=c(col1, col2))
	if(!is.null(filename)) dev.off() 
	# turn off pdf device; graph is created in working directory.
	
}

############# SPLIT VIOLION PLOTS #############
# https://gist.github.com/mbjoseph/5852613
vioplot2log <- function(datas, range=1.5, h = NULL, ylim=NULL, names = NULL,
                      horizontal = FALSE, col = "gray", border = "black", 
                      lty = 1, lwd=1, rectCol = "black", colMed = "black", 
                      pchMed = "|", at, add = FALSE, wex=1, drawRect = TRUE, 
                      side="both",xaxt=NULL, title="",cex.main=1,
                      mar=c(2,5,2,0.5), boxsize=0.05, log=TRUE){
  n <- length(datas)
  if (missing(at)) at <- 1:n
  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)
  q1 <- vector(mode = "numeric", length = n)
  q2 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  med <- vector(mode = "numeric", length = n)
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  baserange <- c(10^-5, -Inf)
  args <- list(display = "none")
  radj <- ifelse(side == "right", 0, 1)
  ladj <- ifelse(side == "left", 0, 1)
  if (!(is.null(h))) args <- c(args, h = h)
  med.dens <- rep(NA, n)
  for (i in 1:n) {
    data <- datas[[i]]
    data.min <- min(data)
    data.max <- max(data)
    q1[i] <- quantile(data, 0.25)
    q2[i] <- quantile(data, 0.5)
    q3[i] <- quantile(data, 0.75)
    med[i] <- median(data)
    iqd <- q3[i] - q1[i]
    upper[i] <- min(q3[i] + range * iqd, data.max)
    lower[i] <- max(q1[i] - range * iqd, data.min)
    est.xlim <- c(min(lower[i], data.min), max(upper[i],data.max))
    smout <- do.call("sm.density", c(list(data, xlim = est.xlim),args))
    med.dat <- do.call("sm.density", c(list(data, xlim=est.xlim,
                              eval.points=med[i], display = "none")))
    med.dens[i] <- med.dat$estimate
    hscale <- 0.4/max(smout$estimate) * wex
    base[[i]] <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    med.dens[i] <- med.dens[i] * hscale
    t <- range(base[[i]])
    baserange[1] <- max(min(baserange[1], t[1],na.rm=T),10^-6)
    baserange[2] <- min(max(baserange[2], t[2],na.rm=T),10^5)
  }
  if (!add) {
    xlim <- if (n == 1)
      at + c(-0.5, 0.5)
    else range(at) + min(diff(at))/2 * c(-1, 1)
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  }
  else {
    label <- names
  }
  boxwidth <- boxsize * wex
  if (!add)
    plot.new()
    par(mar=mar)
  if (!horizontal) {
    if (!add) {
      if (log==TRUE) {
        plot.window(xlim = xlim, ylim = ylim, log="y") 
      } else { 
        plot.window(xlim = xlim, ylim = ylim)
      }
      title(main=title,cex.main=cex.main)
      axis(2)
      axis(1, at=at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(x = c(at[i] - radj*height[[i]], rev(at[i] + ladj*height[[i]])),
              y = c(base[[i]], rev(base[[i]])),
              col = col, border = border,
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd,
              lty = lty)
        rect(at[i] - radj*boxwidth/2,
             q1[i],
             at[i] + ladj*boxwidth/2,
             q3[i], col = rectCol)
        # median line segment
        lines(x = c(at[i] - radj*med.dens[i],
                    at[i],
                    at[i] + ladj*med.dens[i]),
              y = rep(med[i],3))
      }
    }
  }
  else {  # if horizontal
    if (!add) {
      if (log==TRUE) {
        plot.window(xlim = ylim, ylim = xlim, log="x")
      } else {
        plot.window(xlim = ylim, ylim = xlim)
      }
      title(main=title,cex.main=cex.main)
      if (xaxt!="n") axis(1)
      axis(2, at = at, label = label, las=2)
    }
    box()
    for (i in 1:n) {
      polygon(c(base[[i]], rev(base[[i]])),
              c(at[i] - radj*height[[i]], rev(at[i] + ladj*height[[i]])),
              col = col, border = border,
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, lty = lty)
        rect(q1[i], at[i] - radj*boxwidth/2, q3[i], at[i] + ladj*boxwidth/2, 
             col = rectCol)
        lines(y = c(at[i] - radj*med.dens[i],
                    at[i],
                    at[i] + ladj*med.dens[i]),
              x = rep(med[i],3))
      }
    }
  }
  invisible(list(upper = upper, lower = lower, median = med,
                 q1 = q1, q3 = q3))
}
