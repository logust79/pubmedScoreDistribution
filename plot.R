setwd('~/git/pubmedScoreDistribution')
retinal = read.table('report_retinal',sep='\t',header=T)
nonRetinal = read.table('report_nonretinal',sep='\t',header=T)
dr=density(retinal$pubmed_score,from=0)
dn=density(nonRetinal$pubmed_score,from=0)
# look for points of intersection
ps=list(r=data.frame(x=dr$x,y=dr$y),n=data.frame(x=dn$x,y=dn$y))
## Approximate the functions and find intersection
fs <- sapply(ps, function(x) approxfun(x$x, x$y, yleft=0, yright=0))
f <- function(x) fs[[1]](x) - fs[[2]](x)   # function to minimize (difference b/w curves)
meet <- uniroot(f, interval=c(0, 100))$root  # intersection of the two curves
meet
# plot
plot(dr,col=2,main='Density plots of pubmedScores')
lines(dn,col=4)
# draw a line
abline(v=meet, col="orange", lty=2)

legend('topright', legend = c('Retnet genes','Non-Retnet genes'), col = c(2,4),
       lty = 1, merge = TRUE)   #, trace = TRUE)
axis(side=1, at = meet, labels='')
axis(side=1, at = 30, labels=sprintf("%.3f",meet), cex.axis=0.8, tck=0)
