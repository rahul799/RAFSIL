## ----include=TRUE, eval=TRUE, eval=TRUE,messages=FALSE-------------------
library(RAFSIL)
library(grid)
library(gridGraphics)
library(gridExtra)
dat.file =  system.file("extdata", "example_data.rds", package = "RAFSIL",mustWork = TRUE)
dat      = readRDS(dat.file) 
ord = order(dat$labels) ; dat$expr=dat$expr[,ord] ; dat$labels = dat$labels[ord] ; rm(ord)

## ------------------------------------------------------------------------
table(dat$labels)
dim(dat$expr)

## ---- include=TRUE, eval=TRUE, eval=TRUE, messages=FALSE-----------------

#- run RAFSIL1 with 50 forests
res.r1 = RAFSIL(t(dat$expr),nrep = 50, method="RAFSIL1")
res.r2 = RAFSIL(t(dat$expr),           method="RAFSIL2")

#- retriev the dissimilarities
dis.r1  = res.r1$D
dis.r2  = res.r2$D
dis.cor = sqrt((1 - cor(dat$expr,method="spearman"))/2)

## ---- include=TRUE, eval=TRUE, eval=TRUE, messages=FALSE, fig.show="hide"----
plotSIM(dis.r1  ,labels=dat$labels); grid.echo(); g1 = grid.grab()
plotSIM(dis.r2  ,labels=dat$labels); grid.echo(); g2 = grid.grab()
plotSIM(dis.cor ,labels=dat$labels) ; grid.echo(); g3 = grid.grab()

## ---- include=TRUE, eval=TRUE, eval=TRUE, messages=FALSE,fig.show="asis",fig.dim=c(5,2)----
grid.arrange(g1,g2,g3,ncol=3,clip="on")

## ---- include=TRUE, eval=TRUE, eval=TRUE, messages=FALSE,fig.show="asis",fig.dim=c(6,5)----
par(mfrow=c(2,3))
par(mai=c(.1,.1,.5,.1))
plotTSNE(dis.r1,labels=dat$labels,is_distance=TRUE,verbose=FALSE,perplexity=20)
mtext("rafsil-1", line=1)
plotTSNE(dis.r2,labels=dat$labels,is_distance=TRUE,verbose=FALSE,perplexity=20)
mtext("rafsil-2", line=1)
plotTSNE(dis.cor,labels=dat$labels,is_distance=TRUE,verbose=FALSE,perplexity=20)
mtext("correlation", line=1)

plotTSNE(dis.r1,labels=dat$labels,is_distance=FALSE,verbose=FALSE,perplexity=20)
mtext("rafsil-1 / embedding", line=1)
plotTSNE(dis.r2,labels=dat$labels,is_distance=FALSE,verbose=FALSE,perplexity=20)
mtext("rafsil-2 / embedding", line=1)
plotTSNE(dis.cor,labels=dat$labels,is_distance=FALSE,verbose=FALSE,perplexity=20)
mtext("correlation / embedding", line=1)

