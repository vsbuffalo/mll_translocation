## plotting-methods.R - Coverage and other plotting methods

plotCoverage <-
function (x, start=1, end=length(x), col="blue", xlab="Position", ylab="Coverage", ...) {
  xWindow <- as.vector(window(x, start, end))
  x <- start:end
  xlim <- c(start, end)
  ylim <- c(0, max(xWindow))
  plot(x=start, y=0, xlim=xlim, ylim=ylim, xlab=xlab,
       ylab=ylab, type="n", ...)
  polygon(c(start, x, end), c(0, xWindow, 0), col=col)
}

PlotRangesCoverage <-
function(x, chr) {
  plotCoverage(coverage(x[[chr]]), min(start(x[[chr]])), max(end(x[[chr]])), main=chr)
}

