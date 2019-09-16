######################################################
### Code to Replicate Figures of Sim Data          ###
### in Kennedy-Shaffer et al. 2019                 ###
### Update: 9/15/2019                              ###
### Contact: lee_kennedyshaffer@g.harvard.edu      ###
### See Latest Update at:                          ###
### https://github.com/leekshaffer/SW-CRT-analysis ###
######################################################

figfolder <- "Figures"

#### Schematics Figure (Fig 1) ####
if (!require(lattice)) {
  install.packages("lattice")
  library(lattice)
}

schemafolder <- paste0(figfolder,"/SchemaFigs")

NPWPmat <- matrix(c(rep(0,7),1,rep(0,6),rep(1,2),rep(0,5),rep(2,3),rep(3,4),rep(1,4),rep(0,3),rep(1,5),rep(0,2),rep(1,6),0,rep(1,7)), nrow=7, ncol=8, byrow=FALSE)
rownames(NPWPmat) <- seq(1,7,by=1)
colnames(NPWPmat) <- seq(1,8,by=1)

SCmat <- NPWPmat
SCmat[1:3,4] <- 1
SCmat[4,4] <- 4
SCmat[5:6,4] <- 5

COmat <- SCmat
COmat[4:7,3:4] <- 3

COSCmat <- SCmat
COSCmat[4,3] <- 4
COSCmat[5:6,3] <- 5
COSCmat[7,3] <- 3


setEPS()
postscript(file=paste0(figfolder,"/Fig1_Schema.eps"), width=9, height=9, paper="special")

settings <- trellis.par.get()
settings$par.sub.text$font <- 1
settings$axis.components$left$tck <- 0
settings$axis.components$right$tck <- 0
settings$axis.components$top$tck <- 0
settings$axis.components$bottom$tck <- 0
settings$layout.heights$axis.xlab.padding  <- 0
settings$layout.heights$xlab.key.padding <- 1
trellis.par.set(settings)

graphNPWP <- levelplot(x=t(NPWPmat)[1:8,7:1], pretty=TRUE, 
                       xlab="Period", ylab="Cluster",
                       at=seq(-.5,3.5,by=1),
                       col.regions=c("white","darkgreen","green","gray40"),
                       colorkey=FALSE,
                       scales=list(x=list(at=1:8,labels=1:8),
                                   y=list(at=1:7,labels=7:1)),
                       sub="(a) Non-Parametric Within-Period Analysis Method\nfor SW-CRT. The Estimator for Period 4 is the\nUnweighted Mean of the Outcomes in the Upper\nHighlighted Box Compared to the Unweighted Mean\nof the Outcomes in the Lower Highlighted Box.\n",
                       xlim=c(0.5,8.5), ylim=c(0.5,7.5),
                       panel=function(x, ...) {
                         panel.levelplot(x, ...)
                         panel.grid(h=6,v=7,col=1, lty=1, lwd=1)
                         panel.rect(xleft=3.5,ybottom=4.55,xright=4.5,ytop=7.45,
                                    border="orange3", lty=1, lwd=3)
                         panel.rect(xleft=3.5,ybottom=0.55,xright=4.5,ytop=4.45,
                                    border="purple3", lty=1, lwd=3)
                       })

graphSC <- levelplot(x=t(SCmat)[1:8,7:1], pretty=TRUE, 
                       xlab="Period", ylab="Cluster",
                       at=seq(-.5,5.5,by=1),
                       col.regions=c("white","darkgreen","green","gray35","gray15","gray65"),
                       colorkey=FALSE,
                     scales=list(x=list(at=1:8,labels=1:8),
                                 y=list(at=1:7,labels=7:1)),
                       sub="(b) Synthetic Control Analysis Method for SW-CRT.\nThe Estimator for Cluster 2, Period 4 is the Outcome\nin the Upper Highlighted Box Compared to the Weighted\nMean of the Outcomes in the Lower Highlighted Box.\nThe Shading within the Lower Box Indicates \nDifferential Weights.",
                     xlim=c(0.5,8.5), ylim=c(0.5,7.5),
                     panel=function(x, ...) {
                       panel.levelplot(x, ...)
                       panel.grid(h=6,v=7,col=1, lty=1, lwd=1)
                       panel.rect(xleft=3.5,ybottom=5.5,xright=4.5,ytop=6.5,
                                  border="orange3", lty=1, lwd=3)
                       panel.rect(xleft=3.5,ybottom=0.55,xright=4.5,ytop=4.45,
                                  border="purple3", lty=1, lwd=3)})

graphCO <- levelplot(x=t(COmat)[1:8,7:1], pretty=TRUE, 
                     xlab="Period", ylab="Cluster",
                     at=seq(-.5,5.5,by=1),
                     col.regions=c("white","darkgreen","green","gray35","gray15","gray65"),
                     colorkey=FALSE,
                     scales=list(x=list(at=1:8,labels=1:8),
                                 y=list(at=1:7,labels=7:1)),
                     sub="(c) Crossover Analysis Method for SW-CRT.\nThe Estimator for Period 4 is the Difference between\nthe Contrast of Outcomes in the Topmost Highlighted\nBox and the Unweighted Mean of the Contrasts\nof Outcomes in the Lower Highlighted Boxes.\n",
                     xlim=c(0.5,8.5), ylim=c(0.5,7.5),
                     panel=function(x, ...) {
                       panel.levelplot(x, ...)
                       panel.grid(h=6,v=7,col=1, lty=1, lwd=1)
                       panel.rect(xleft=2.5,ybottom=4.55,xright=4.5,ytop=5.5,
                                  border="orange3", lty=1, lwd=3)
                       panel.rect(xleft=2.5,ybottom=3.55,xright=4.5,ytop=4.45,
                                  border="purple3", lty=1, lwd=3)
                       panel.rect(xleft=2.5,ybottom=2.55,xright=4.5,ytop=3.45,
                                  border="purple3", lty=1, lwd=3)
                       panel.rect(xleft=2.5,ybottom=1.55,xright=4.5,ytop=2.45,
                                  border="purple3", lty=1, lwd=3)
                       panel.rect(xleft=2.5,ybottom=0.55,xright=4.5,ytop=1.45,
                                  border="purple3", lty=1, lwd=3)})

graphCOSC <- levelplot(x=t(COSCmat)[1:8,7:1], pretty=TRUE, 
                     xlab="Period", ylab="Cluster",
                     at=seq(-.5,5.5,by=1),
                     col.regions=c("white","darkgreen","green","gray35","gray15","gray65"),
                     colorkey=FALSE,
                     scales=list(x=list(at=1:8,labels=1:8),
                                 y=list(at=1:7,labels=7:1)),
                     sub="(d) Crossover-Synthetic Control Analysis Method\nfor SW-CRT. The Estimator for Period 4 is the Difference\nbetween the Contrast of Outcomes in the Topmost High-\nlighted Box and the Weighted Mean of the Contrasts of\nOutcomes in the Lower Highlighted Boxes. The Shading\nwithin the Lower Boxes Indicates Differential Weights.",
                     xlim=c(0.5,8.5), ylim=c(0.5,7.5),
                     panel=function(x, ...) {
                       panel.levelplot(x, ...)
                       panel.grid(h=6,v=7,col=1, lty=1, lwd=1)
                       panel.rect(xleft=2.5,ybottom=4.55,xright=4.5,ytop=5.5,
                                  border="orange3", lty=1, lwd=3)
                       panel.rect(xleft=2.5,ybottom=3.55,xright=4.5,ytop=4.45,
                                  border="purple3", lty=1, lwd=3)
                       panel.rect(xleft=2.5,ybottom=2.55,xright=4.5,ytop=3.45,
                                  border="purple3", lty=1, lwd=3)
                       panel.rect(xleft=2.5,ybottom=1.55,xright=4.5,ytop=2.45,
                                  border="purple3", lty=1, lwd=3)
                       panel.rect(xleft=2.5,ybottom=0.55,xright=4.5,ytop=1.45,
                                  border="purple3", lty=1, lwd=3)})

plot(graphNPWP, split=c(1,1,2,2))
plot(graphSC, split=c(2,1,2,2), newpage=FALSE)
plot(graphCO, split=c(1,2,2,2), newpage=FALSE)
plot(graphCOSC, split=c(2,2,2,2), newpage=FALSE)

dev.off()



#### RD Figures (Figs 2--6) ####

scentext =rep(c("(a) Scenario 1. Common Time Effects, No Cluster-Period Effect",
                "(b) Scenario 2. Common Time Effecs, Cluster-Period Effect",
                "(c) Scenario 3. Varying Time Effects, No Cluster-Period Effect",
                "(d) Scenario 4. Varying Time Effects, Cluster-Period Effect"),
              6)
subgraphs <- c("a","b","c","d")
setEPS()
postscript(file=paste0(figfolder,"/Fig2_RD_GenData.eps"), width=9, height=9, paper="special")
par(mfrow=c(2,2))
for (i in 1:4) {
  namei <- paste0("Data_Fig2",subgraphs[i])
  load(paste0("Fig_Data/",namei,".Rda"))
  df <- get(namei)
  plot(x=NA, y=NA, xlim=c(.75,8.25),ylim=c(0,.8), 
       xaxt="n", xlab="", ylab="",
       sub=scentext[i], cex.sub=1)
  axis(side=1, at=seq(1,8,by=1), cex.axis=1)
  title(xlab="Period", line=2)
  title(ylab=expression(Y[ij]), line=2)
  for (j in 1:(dim(df)[1])) {
    lines(x=as.numeric(names(df)[-1]), y=df[j,-1],
          type="l", col=j+1, lwd=2)
    points(x=as.numeric(names(df)[-1]), y=df[j,-1],
           pch=19, col=j+1, cex=1.5)
  }
}
dev.off()

PCHs <- c(21,21,22,22,23,24,24,25,25,25,8,8,9)
Types <- c("MEM","CPI","MEM-a","CPI-a","NPWP","SCSWT1","SCSWT2","CO.Ctrl","CO.CtWt","CO.Both",
           "COSC1","COSC2","ENS")
NamesInf <- c("MEM","CPI","MEM-a","CPI-a","NPWP","SC-1","SC-2","CO-1","CO-2","CO-3",
              "COSC","COSC-2","ENS")
NamesEst <- c("MEM","CPI","MEM","CPI","NPWP","SC-1","SC-2","CO-1","CO-2","CO-3",
              "COSC","COSC-2","ENS")
scens <- c(1,2,3,4,9,10,11,12,17,18,19,20)


Plot.Res.Ests <- function(ResComb,scens,plotTitle,PlotTypes,ypositions,ylabs,inlabs,Truth,legend=TRUE,yaxis=TRUE,colorvec) {
  minx <- floor(min(ResComb$MinEst[ResComb$Scen %in% scens & ResComb$Type %in% PlotTypes],na.rm=TRUE)*20)/20
  maxx <- ceiling(max(ResComb$MaxEst[ResComb$Scen %in% scens & ResComb$Type %in% PlotTypes],na.rm=TRUE)*20)/20
  plot(x=NA,y=NA,
       xlim=c(minx,maxx),ylim=c(0,1), 
       xlab="", ylab="", 
       xaxt="n",yaxt="n",
       sub=plotTitle, cex.sub=1)
  title(xlab="Estimated Risk Difference", line=2, cex.lab=.9)
  
  if (yaxis==TRUE) {
    axis(side=2, at=ypositions[order(ypositions)],
         labels=ylabs,
         cex.axis=1.3,las=2)
  }
  axis(side=1, at=round(seq(minx,maxx,by=.025), digits=3), cex.axis=.75)
  abline(v=Truth, col=8, lty=2, lwd=2)
  xmean = mean(c(minx,maxx))
  text(x=xmean, y=1, labels=inlabs[1], cex=.9)
  text(x=xmean, y=.5, labels=inlabs[2], cex=.9)
  
  
  pos2s <- seq(1/15,-1/15,length.out=length(PlotTypes))
  
  for (i in 1:length(scens)) {
    sceni <- scens[i]
    posi <- ypositions[i]
    if (yaxis==FALSE) {text(x=xmean, y=posi+1/10, labels=ylabs[i], cex=.8)}
    for (j in 1:length(PlotTypes)) {
      typej <- PlotTypes[j]
      pchj <- PCHs[match(typej,Types)]
      rowij <- ResComb[ResComb$Scen==sceni & ResComb$Type==typej,]
      lines(x=c(rowij$MinEst,rowij$MaxEst),y=rep(posi+pos2s[j],2), type="l",
            col=colorvec[j], lwd=1.5)
      points(x=rowij$Mean, y=posi+pos2s[j], col=colorvec[j], bg=colorvec[j], pch=pchj, cex=1)
    }
  }
  
  if (legend==TRUE) {
    legnames <- NamesEst[match(PlotTypes,Types)]
    pchvec <- PCHs[match(PlotTypes,Types)]
    legend(x="topright", legend=legnames,
           col=colorvec[seq(1,length(legnames),by=1)],
           pt.bg=colorvec[seq(1,length(legnames),by=1)],
           lty=rep(1,length(legnames)),
           lwd=rep(1,length(legnames)),
           pch=pchvec,
           horiz=FALSE, cex=.75)
  }
}

Plot.Res.Power <- function(ResComb,scens,PlotTypes,ypositions,ylabs,inlabs,xrange=NULL,xlabel="Power",Nominal=NA,legend=TRUE,yaxis=TRUE, rmTIE=1, rmTIEscens=NA, colorvec, binwidth=.01) {
  if (is.null(xrange)) {
    minx <- min(floor(min(ResComb$Power[ResComb$Scen %in% scens & ResComb$Type %in% PlotTypes],na.rm=TRUE)*20)/20)
    maxx <- max(ceiling(max(ResComb$Power[ResComb$Scen %in% scens & ResComb$Type %in% PlotTypes],na.rm=TRUE)*20)/20)
  } else {
    minx <- min(xrange)
    maxx <- max(xrange)
  }
  
  plot(x=NA, y=NA,
       xlim=c(minx,maxx), ylim=c(0,1), 
       xlab="", ylab="", 
       xaxt="n",yaxt="n")
  if (yaxis==TRUE) {
    axis(side=2, at=ypositions[order(ypositions)],
         labels=ylabs,
         cex.axis=1.3,las=2)
  }
  axis(side=1, at=round(seq(minx,maxx,by=.05), digits=2), labels=paste0(seq(minx*100,maxx*100,by=5),"%"), 
       cex.axis=.75)
  title(xlab=xlabel, line=2, cex.lab=.9)
  
  if (!is.na(Nominal)) {abline(v=Nominal, col=8, lwd=2, lty=2)}
  xmean = mean(c(minx,maxx))
  text(x=xmean,y=1,labels=inlabs[1], cex=.9)
  text(x=xmean,y=.5,labels=inlabs[2], cex=.9)
  
  for (i in 1:length(scens)) {
    if (sum(is.na(rmTIEscens))==0) {
      TIEsceni <- rmTIEscens[i]
      TIEclmni <- ResComb[ResComb$Scen==TIEsceni & ResComb$Type %in% PlotTypes,c("Type","Power")]
      rmTIElist <- TIEclmni$Type[TIEclmni$Power > rmTIE]
      PlotTypesi <- PlotTypes[!(PlotTypes %in% rmTIElist)]
    } else {
      PlotTypesi <- PlotTypes
    }
    sceni <- scens[i]
    posi <- ypositions[i]
    rowi <- NULL
    for (j in 1:length(PlotTypesi)) {
      rowi <- append(rowi, ResComb$Power[ResComb$Scen==sceni & ResComb$Type==PlotTypesi[j]])
    }
    bini <- ceiling((rowi-minx)/binwidth)
    bincti <- rep(1,length(bini))
    binmaxi <- rep(0,length(bini))
    for (j in 2:length(bini)) {
      bincti[j] <- sum(bini[1:j]==bini[j])
    }
    for (j in 1:length(bini)) {
      binmaxi[j] <- max(bincti[bini==bini[j]])
    }
    lines(x=c(minx,maxx),y=rep(posi,2),col=8, lwd=2, lty=3)
    if (yaxis==FALSE) {text(x=xmean, y=posi+1/10, labels=ylabs[i], cex=.8)}
    for (j in 1:length(PlotTypesi)) {
      typej <- PlotTypesi[j]
      rowij <- ResComb[ResComb$Scen==sceni & ResComb$Type==typej,]
      pchj <- PCHs[match(typej,Types)]
      yposij <- posi - .014*(binmaxi[j]-1)/2 + .014*(bincti[j]-1)
      points(x=rowij$Power, y=yposij, 
             col=colorvec[j], bg=colorvec[j], pch=pchj, cex=1.5)
    }
  }
  
  if (legend==TRUE) {
    legnames <- NamesInf[match(PlotTypes,Types)]
    pchvec <- PCHs[match(PlotTypes,Types)]
    legend(x="bottom", legend=legnames,
           col=colorvec[seq(1,length(legnames),by=1)],
           pt.bg=colorvec[seq(1,length(legnames),by=1)],
           pch=pchvec,
           horiz=TRUE, cex=.6)
  }
}


Plot.Res.Coverage <- function(ResComb,scens,plotTitle,PlotTypes,ypositions,ylabs,inlabs,Nominal=.95,legend=TRUE,yaxis=TRUE, fixX=FALSE, colorvec, binwidth=.003) {
  if (fixX==TRUE) {
    minx <- .9
    maxx <- 1
  } else {
    minx <- min(floor(min(ResComb$Coverage[ResComb$Scen %in% scens & ResComb$Type %in% PlotTypes],na.rm=TRUE)*20)/20)
    maxx <- max(ceiling(max(ResComb$Coverage[ResComb$Scen %in% scens & ResComb$Type %in% PlotTypes],na.rm=TRUE)*20)/20)
  }
  
  plot(x=NA, y=NA,
       xlim=c(minx,maxx), ylim=c(0,1), 
       xlab="", ylab="", xaxt="n",yaxt="n",
       sub=plotTitle, cex.sub=1)
  if(yaxis==TRUE) {
    axis(side=2, at=ypositions[order(ypositions)],
         labels=ylabs,
         cex.axis=1.3,las=2)
  }
  if (fixX==TRUE) {
    axis(side=1, at=round(seq(minx,maxx,by=.02), digits=2), labels=paste0(seq(minx*100,maxx*100,by=2),"%"),
         cex.axis=.75)
  } else {
    axis(side=1, at=round(seq(minx,maxx,by=.05), digits=2), labels=paste0(seq(minx*100,maxx*100,by=5),"%"),
         cex.axis=.75)
  }
  title(xlab="95% CI Coverage", line=2, cex.lab=.9)
  
  if (!is.na(Nominal)) {abline(v=Nominal, col=8, lwd=2, lty=2)}
  xmean = mean(c(minx,maxx))
  text(x=xmean,y=1,labels=inlabs[1], cex=.9)
  text(x=xmean,y=.5,labels=inlabs[2], cex=.9)
  
  for (i in 1:length(scens)) {
    sceni <- scens[i]
    posi <- ypositions[i]
    rowi <- NULL
    for (j in 1:length(PlotTypes)) {
      rowi <- append(rowi, ResComb$Coverage[ResComb$Scen==sceni & ResComb$Type==PlotTypes[j]])
    }
    bini <- ceiling((rowi-minx)/binwidth)
    bincti <- rep(1,length(bini))
    binmaxi <- rep(0,length(bini))
    for (j in 2:length(bini)) {
      bincti[j] <- sum(bini[1:j]==bini[j])
    }
    for (j in 1:length(bini)) {
      binmaxi[j] <- max(bincti[bini==bini[j]])
    }
    lines(x=c(minx,maxx),y=rep(posi,2),col=8, lwd=2, lty=3)
    if (yaxis==FALSE) {text(x=xmean, y=posi+1/10, labels=ylabs[i], cex=.8)}
    for (j in 1:length(PlotTypes)) {
      typej <- PlotTypes[j]
      rowij <- ResComb[ResComb$Scen==sceni & ResComb$Type==typej,]
      pchj <- PCHs[match(typej,Types)]
      yposij <- posi - .014*(binmaxi[j]-1)/2 + .014*(bincti[j]-1)
      points(x=rowij$Coverage, y=yposij, 
             col=colorvec[j], bg=colorvec[j], pch=pchj, cex=1.5)
    }
  }
  
  if (legend==TRUE) {
    legnames <- NamesInf[match(PlotTypes,Types)]
    pchvec <- PCHs[match(PlotTypes,Types)]
    locat <- ifelse(fixX==TRUE,"right","topleft")
    legend(x=locat, legend=legnames,
           col=colorvec[seq(1,length(legnames),by=1)],
           pt.bg=colorvec[seq(1,length(legnames),by=1)],
           pch=pchvec,
           horiz=FALSE, cex=.75)
  }
}


typesEst <- c("MEM","CPI","NPWP","SCSWT1","SCSWT2","CO.Ctrl","CO.CtWt","CO.Both","COSC1","ENS")
typesInf <- c(typesEst,"MEM-a","CPI-a")
cols <- c("#a6cee3","#1f78b4","#cab2d6","#b2df8a","#33a02c","#fb9a99","#e31a1c","#ff7f00","#fdbf6f","#a6761d","#a6cee3","#1f78b4")
scens1 <- 1:4
scens9 <- 9:12
scens17 <- 17:20
yposs <- c(12/14,9/14,5/14,2/14)
ylabels <- c(expression(nu*"=0"),expression(nu*"=0.01"),
             expression(nu*"=0"),
             expression(nu*"=0.01"))
inlabels <- c(expression("Common Time Effects"),
              expression("Varying Time Effects"))


load("Fig_Data/Data_Fig3a.Rda")
load("Fig_Data/Data_Fig3b.Rda")
load("Fig_Data/Data_Fig4.Rda")
load("Fig_Data/Data_Fig5a.Rda")
load("Fig_Data/Data_Fig5b.Rda")
load("Fig_Data/Data_Fig6.Rda")

scentext =c(expression('(a) Moderate Treatment Effect ('*beta*' = -0.1)'),
            expression('(b) No Treatment Effect ('*beta*' = 0)'))
setEPS()
postscript(file=paste0(figfolder,"/Fig3_RD_Ests.eps"), width=9, height=5, paper="special")
par(mfrow=c(1,2))
Plot.Res.Ests(Data_Fig3a, scens9, scentext[1], typesEst, yposs, ylabels, inlabels, -.1, legend=TRUE, yaxis=FALSE, colorvec=cols)
Plot.Res.Ests(Data_Fig3b, scens17, scentext[2], typesEst, yposs, ylabels, inlabels, 0, legend=TRUE, yaxis=FALSE, colorvec=cols)
dev.off()

setEPS()
postscript(file=paste0(figfolder,"/Fig4_RD_TIE.eps"), width=7.5, height=5, paper="special")
par(mfrow=c(1,1))
Plot.Res.Power(Data_Fig4, scens17, typesInf, yposs, ylabels, inlabels, xrange=NULL, "Type I Error", .05, legend=TRUE, yaxis=FALSE, colorvec=cols)

dev.off()

setEPS()
postscript(file=paste0(figfolder,"/Fig5_RD_Cvg.eps"), width=9, height=5, paper="special")
par(mfrow=c(1,2))
Plot.Res.Coverage(Data_Fig5a, scens9, scentext[1], typesInf, yposs, ylabels, inlabels, legend=TRUE, yaxis=FALSE, fixX=FALSE, colorvec=cols)
Plot.Res.Coverage(Data_Fig5b, scens17, scentext[2], typesInf, yposs, ylabels, inlabels, legend=TRUE, yaxis=FALSE, fixX=FALSE, colorvec=cols)
dev.off()

setEPS()
postscript(file=paste0(figfolder,"/Fig6_RD_Pwr.eps"), width=7.5, height=5, paper="special")
par(mfrow=c(1,1))
Plot.Res.Power(Data_Fig6, scens9, typesInf, yposs, ylabels, inlabels, xrange=c(.35,1), legend=TRUE, yaxis=FALSE, rmTIE=0.1, rmTIEscens=scens17, colorvec=cols)
dev.off()


#### OR Figures (Figs 7--11) ####

scentext =rep(c("(a) Scenario 1. Common Time Effects, No Cluster-Period Effect",
                "(b) Scenario 2. Common Time Effecs, Cluster-Period Effect",
                "(c) Scenario 3. Varying Time Effects, No Cluster-Period Effect",
                "(d) Scenario 4. Varying Time Effects, Cluster-Period Effect"),
              6)
setEPS()
postscript(file=paste0(figfolder,"/Fig7_OR_GenData.eps"), width=9, height=9, paper="special")
par(mfrow=c(2,2))
for (i in 1:4) {
  namei <- paste0("Data_Fig7",subgraphs[i])
  load(paste0("Fig_Data/",namei,".Rda"))
  df <- get(namei)
  
  plot(x=NA, y=NA, xlim=c(.75,8.25),ylim=c(0,.8), 
       xaxt="n", xlab="", ylab="",
       sub=scentext[i], cex.sub=1)
  axis(side=1, at=seq(1,8,by=1), cex.axis=1)
  title(xlab="Period", line=2)
  title(ylab=expression(Y[ij]), line=2)
  for (j in 1:(dim(df)[1])) {
    lines(x=as.numeric(names(df)[-1]), y=df[j,-1],
          type="l", col=j+1, lwd=2)
    points(x=as.numeric(names(df)[-1]), y=df[j,-1],
           pch=19, col=j+1, cex=1.5)
  }
}
dev.off()

PCHs <- c(21,21,22,22,23,24,24,25,25,25,8,8,9)
Types <- c("MEM","CPI","MEM-a","CPI-a","NPWP","SCSWT1","SCSWT2","CO.Ctrl","CO.CtWt","CO.Both",
           "COSC1","COSC2","ENS")
NamesInf <- c("MEM","CPI","MEM-a","CPI-a","NPWP","SC-1","SC-2","CO-1","CO-2","CO-3",
              "COSC","COSC-2","ENS")
NamesEst <- c("MEM","CPI","MEM","CPI","NPWP","SC-1","SC-2","CO-1","CO-2","CO-3",
              "COSC","COSC-2","ENS")
scens <- c(1,2,3,4,9,10,11,12,17,18,19,20)


Plot.Res.Ests <- function(ResComb,scens,plotTitle,PlotTypes,ypositions,ylabs,inlabs,Truth,legend=TRUE,yaxis=TRUE,colorvec,minx=NULL,maxx=NULL) {
  if (is.null(minx)) {
    minx <- floor(min(ResComb$MinEst[ResComb$Scen %in% scens & ResComb$Type %in% PlotTypes],na.rm=TRUE)*20)/20
  }
  if (is.null(maxx)) {
    maxx <- ceiling(max(ResComb$MaxEst[ResComb$Scen %in% scens & ResComb$Type %in% PlotTypes],na.rm=TRUE)*20)/20
  }
  plot(x=NA,y=NA,
       xlim=c(minx,maxx),ylim=c(0,1), 
       xlab="", ylab="", 
       xaxt="n",yaxt="n",
       sub=plotTitle, cex.sub=1)
  title(xlab="Estimated Log Odds Ratio", line=2, cex.lab=.9)
  
  if (yaxis==TRUE) {
    axis(side=2, at=ypositions[order(ypositions)],
         labels=ylabs,
         cex.axis=1.3,las=2)
  }
  xlabs <- round(seq(minx,maxx,by=.05), digits=2)
  axis(side=1, at=xlabs,
       labels=format(xlabs,digits=2,nsmall=2),
       cex.axis=.75)
  abline(v=Truth, col=8, lty=2, lwd=2)
  xmean = mean(c(minx,maxx))
  text(x=xmean, y=1, labels=inlabs[1], cex=.9)
  text(x=xmean, y=.5, labels=inlabs[2], cex=.9)
  
  
  pos2s <- seq(1/15,-1/15,length.out=length(PlotTypes))
  
  for (i in 1:length(scens)) {
    sceni <- scens[i]
    posi <- ypositions[i]
    if (yaxis==FALSE) {text(x=xmean, y=posi+1/10, labels=ylabs[i], cex=.8)}
    for (j in 1:length(PlotTypes)) {
      typej <- PlotTypes[j]
      pchj <- PCHs[match(typej,Types)]
      rowij <- ResComb[ResComb$Scen==sceni & ResComb$Type==typej,]
      lines(x=c(rowij$MinEst,rowij$MaxEst),y=rep(posi+pos2s[j],2), type="l",
            col=colorvec[j], lwd=1.5)
      points(x=rowij$Mean, y=posi+pos2s[j], col=colorvec[j], bg=colorvec[j], pch=pchj, cex=1)
    }
  }
  
  if (legend==TRUE) {
    legnames <- NamesEst[match(PlotTypes,Types)]
    pchvec <- PCHs[match(PlotTypes,Types)]
    legend(x="topright", legend=legnames,
           col=colorvec[seq(1,length(legnames),by=1)],
           pt.bg=colorvec[seq(1,length(legnames),by=1)],
           lty=rep(1,length(legnames)),
           lwd=rep(1,length(legnames)),
           pch=pchvec,
           horiz=FALSE, cex=.75)
  }
}

Plot.Res.Power <- function(ResComb,scens,PlotTypes,ypositions,ylabs,inlabs,xrange=NULL,xlabel="Power",Nominal=NA,legend=TRUE,yaxis=TRUE, rmTIE=1, rmTIEscens=NA, colorvec, binwidth=.01) {
  if (is.null(xrange)) {
    minx <- min(floor(min(ResComb$Power[ResComb$Scen %in% scens & ResComb$Type %in% PlotTypes],na.rm=TRUE)*20)/20)
    maxx <- max(ceiling(max(ResComb$Power[ResComb$Scen %in% scens & ResComb$Type %in% PlotTypes],na.rm=TRUE)*20)/20)
  } else {
    minx <- min(xrange)
    maxx <- max(xrange)
  }
  
  plot(x=NA, y=NA,
       xlim=c(minx,maxx), ylim=c(0,1), 
       xlab="", ylab="", 
       xaxt="n",yaxt="n")
  if (yaxis==TRUE) {
    axis(side=2, at=ypositions[order(ypositions)],
         labels=ylabs,
         cex.axis=1.3,las=2)
  }
  axis(side=1, at=round(seq(minx,maxx,by=.05), digits=2), labels=paste0(seq(minx*100,maxx*100,by=5),"%"), 
       cex.axis=.75)
  title(xlab=xlabel, line=2, cex.lab=.9)
  
  if (!is.na(Nominal)) {abline(v=Nominal, col=8, lwd=2, lty=2)}
  xmean = mean(c(minx,maxx))
  text(x=xmean,y=1,labels=inlabs[1], cex=.9)
  text(x=xmean,y=.5,labels=inlabs[2], cex=.9)
  
  for (i in 1:length(scens)) {
    if (sum(is.na(rmTIEscens))==0) {
      TIEsceni <- rmTIEscens[i]
      TIEclmni <- ResComb[ResComb$Scen==TIEsceni & ResComb$Type %in% PlotTypes,c("Type","Power")]
      rmTIElist <- TIEclmni$Type[TIEclmni$Power > rmTIE]
      PlotTypesi <- PlotTypes[!(PlotTypes %in% rmTIElist)]
    } else {
      PlotTypesi <- PlotTypes
    }
    sceni <- scens[i]
    posi <- ypositions[i]
    rowi <- NULL
    for (j in 1:length(PlotTypesi)) {
      rowi <- append(rowi, ResComb$Power[ResComb$Scen==sceni & ResComb$Type==PlotTypesi[j]])
    }
    bini <- ceiling((rowi-minx)/binwidth)
    bincti <- rep(1,length(bini))
    binmaxi <- rep(0,length(bini))
    for (j in 2:length(bini)) {
      bincti[j] <- sum(bini[1:j]==bini[j])
    }
    for (j in 1:length(bini)) {
      binmaxi[j] <- max(bincti[bini==bini[j]])
    }
    lines(x=c(minx,maxx),y=rep(posi,2),col=8, lwd=2, lty=3)
    if (yaxis==FALSE) {text(x=xmean, y=posi+1/10, labels=ylabs[i], cex=.8)}
    for (j in 1:length(PlotTypesi)) {
      typej <- PlotTypesi[j]
      rowij <- ResComb[ResComb$Scen==sceni & ResComb$Type==typej,]
      pchj <- PCHs[match(typej,Types)]
      yposij <- posi - .014*(binmaxi[j]-1)/2 + .014*(bincti[j]-1)
      points(x=rowij$Power, y=yposij, 
             col=colorvec[j], bg=colorvec[j], pch=pchj, cex=1.5)
    }
  }
  
  if (legend==TRUE) {
    legnames <- NamesInf[match(PlotTypes,Types)]
    pchvec <- PCHs[match(PlotTypes,Types)]
    legend(x="bottom", legend=legnames,
           col=colorvec[seq(1,length(legnames),by=1)],
           pt.bg=colorvec[seq(1,length(legnames),by=1)],
           pch=pchvec,
           horiz=TRUE, cex=.6)
  }
}


Plot.Res.Coverage <- function(ResComb,scens,plotTitle,PlotTypes,ypositions,ylabs,inlabs,Nominal=.95,legend=TRUE,yaxis=TRUE, fixX=FALSE, colorvec, binwidth=.003) {
  if (fixX==TRUE) {
    minx <- .9
    maxx <- 1
  } else {
    minx <- min(floor(min(ResComb$Coverage[ResComb$Scen %in% scens & ResComb$Type %in% PlotTypes],na.rm=TRUE)*20)/20)
    maxx <- max(ceiling(max(ResComb$Coverage[ResComb$Scen %in% scens & ResComb$Type %in% PlotTypes],na.rm=TRUE)*20)/20)
  }
  
  plot(x=NA, y=NA,
       xlim=c(minx,maxx), ylim=c(0,1), 
       xlab="", ylab="", xaxt="n",yaxt="n",
       sub=plotTitle, cex.sub=1)
  if(yaxis==TRUE) {
    axis(side=2, at=ypositions[order(ypositions)],
         labels=ylabs,
         cex.axis=1.3,las=2)
  }
  if (fixX==TRUE) {
    axis(side=1, at=round(seq(minx,maxx,by=.02), digits=2), labels=paste0(seq(minx*100,maxx*100,by=2),"%"),
         cex.axis=.75)
  } else {
    axis(side=1, at=round(seq(minx,maxx,by=.05), digits=2), labels=paste0(seq(minx*100,maxx*100,by=5),"%"),
         cex.axis=.75)
  }
  title(xlab="95% CI Coverage", line=2, cex.lab=.9)
  
  if (!is.na(Nominal)) {abline(v=Nominal, col=8, lwd=2, lty=2)}
  xmean = mean(c(minx,maxx))
  text(x=xmean,y=1,labels=inlabs[1], cex=.9)
  text(x=xmean,y=.5,labels=inlabs[2], cex=.9)
  
  for (i in 1:length(scens)) {
    sceni <- scens[i]
    posi <- ypositions[i]
    rowi <- NULL
    for (j in 1:length(PlotTypes)) {
      rowi <- append(rowi, ResComb$Coverage[ResComb$Scen==sceni & ResComb$Type==PlotTypes[j]])
    }
    bini <- ceiling((rowi-minx)/binwidth)
    bincti <- rep(1,length(bini))
    binmaxi <- rep(0,length(bini))
    for (j in 2:length(bini)) {
      bincti[j] <- sum(bini[1:j]==bini[j])
    }
    for (j in 1:length(bini)) {
      binmaxi[j] <- max(bincti[bini==bini[j]])
    }
    lines(x=c(minx,maxx),y=rep(posi,2),col=8, lwd=2, lty=3)
    if (yaxis==FALSE) {text(x=xmean, y=posi+1/10, labels=ylabs[i], cex=.8)}
    for (j in 1:length(PlotTypes)) {
      typej <- PlotTypes[j]
      rowij <- ResComb[ResComb$Scen==sceni & ResComb$Type==typej,]
      pchj <- PCHs[match(typej,Types)]
      yposij <- posi - .014*(binmaxi[j]-1)/2 + .014*(bincti[j]-1)
      points(x=rowij$Coverage, y=yposij, 
             col=colorvec[j], bg=colorvec[j], pch=pchj, cex=1.5)
    }
  }
  
  if (legend==TRUE) {
    legnames <- NamesInf[match(PlotTypes,Types)]
    pchvec <- PCHs[match(PlotTypes,Types)]
    locat <- ifelse(fixX==TRUE,"right","topleft")
    legend(x=locat, legend=legnames,
           col=colorvec[seq(1,length(legnames),by=1)],
           pt.bg=colorvec[seq(1,length(legnames),by=1)],
           pch=pchvec,
           horiz=FALSE, cex=.75)
  }
}


typesEst <- c("MEM","CPI","NPWP","SCSWT1","SCSWT2","CO.Ctrl","CO.CtWt","CO.Both","COSC1","ENS")
typesInf <- c(typesEst,"MEM-a","CPI-a")
cols <- c("#a6cee3","#1f78b4","#cab2d6","#b2df8a","#33a02c","#fb9a99","#e31a1c","#ff7f00","#fdbf6f","#a6761d","#a6cee3","#1f78b4")
scens1 <- 1:4
scens9 <- 9:12
scens17 <- 17:20
yposs <- c(12/14,9/14,5/14,2/14)
ylabels <- c(expression(nu*"=0"),expression(nu*"=0.01"),
             expression(nu*"=0"),
             expression(nu*"=0.01"))
inlabels <- c(expression("Common Time Effects"),
              expression("Varying Time Effects"))

load("Fig_Data/Data_Fig8a.Rda")
load("Fig_Data/Data_Fig8b.Rda")
load("Fig_Data/Data_Fig9.Rda")
load("Fig_Data/Data_Fig10a.Rda")
load("Fig_Data/Data_Fig10b.Rda")
load("Fig_Data/Data_Fig11.Rda")


scentext =c(expression('(a) Moderate Treatment Effect ('*beta*' = log(0.66) ' %~~% ' -0.416)'),
            expression('(b) No Treatment Effect ('*beta*' = log(1) = 0)'))
setEPS()
postscript(file=paste0(figfolder,"/Fig8_OR_Ests.eps"), width=9, height=5, paper="special")
par(mfrow=c(1,2))
Plot.Res.Ests(Data_Fig8a, scens9, scentext[1], typesEst, yposs, ylabels, inlabels, log(.66), legend=TRUE, yaxis=FALSE, colorvec=cols, minx=-.6, maxx=-.25)
Plot.Res.Ests(Data_Fig8b, scens17, scentext[2], typesEst, yposs, ylabels, inlabels, 0, legend=TRUE, yaxis=FALSE, colorvec=cols, minx=-.15, maxx=.15)
dev.off()

setEPS()
postscript(file=paste0(figfolder,"/Fig9_OR_TIE.eps"), width=7.5, height=5, paper="special")
par(mfrow=c(1,1))
Plot.Res.Power(Data_Fig9, scens17, typesInf, yposs, ylabels, inlabels, xrange=NULL, "Type I Error", .05, legend=TRUE, yaxis=FALSE, colorvec=cols)
dev.off()

setEPS()
postscript(file=paste0(figfolder,"/Fig10_OR_Cvg.eps"), width=9, height=5, paper="special")
par(mfrow=c(1,2))
Plot.Res.Coverage(Data_Fig10a, scens9, scentext[1], typesInf, yposs, ylabels, inlabels, legend=TRUE, yaxis=FALSE, fixX=FALSE, colorvec=cols)
Plot.Res.Coverage(Data_Fig10b, scens17, scentext[2], typesInf, yposs, ylabels, inlabels, legend=TRUE, yaxis=FALSE, fixX=FALSE, colorvec=cols)
dev.off()

setEPS()
postscript(file=paste0(figfolder,"/Fig11_OR_Pwr.eps"), width=7.5, height=5, paper="special")
par(mfrow=c(1,1))
Plot.Res.Power(Data_Fig11, scens9, typesInf, yposs, ylabels, inlabels, xrange=c(.35,1), legend=TRUE, yaxis=FALSE, rmTIE=0.1, rmTIEscens=scens17, colorvec=cols)
dev.off()



#### Covariance Figures (Figs 12--13) ####

load("Fig_Data/Data_Fig12.Rda")
rownames(Data_Fig12) <- c("MEM","CPI","NPWP","SC-1","SC-2","CO-1","CO-2","CO-3","COSC","ENS")
colnames(Data_Fig12) <- rownames(Data_Fig12)

load("Fig_Data/Data_Fig13.Rda")
rownames(Data_Fig13) <- rownames(Data_Fig12)
colnames(Data_Fig13) <- colnames(Data_Fig12)

setEPS()
postscript(file=paste0(figfolder,"/Fig12_RD_Covariance.eps"), width=9, height=9, paper="special")
levelplot(x=Data_Fig12[1:10,10:1],
          pretty=TRUE, 
          xlab="Method", ylab="Method",
          at=seq(0,.00425,by=.00025),
          colorkey=list(at=seq(0,.00425,by=.00025)),
          col.regions=heat.colors(100)[length(heat.colors(100)):1])
dev.off()

setEPS()
postscript(file=paste0(figfolder,"/Fig13_OR_Covariance.eps"), width=9, height=9, paper="special")
levelplot(x=Data_Fig13[1:10,10:1],
          pretty=TRUE, 
          xlab="Method", ylab="Method", 
          at=seq(0,.0425,by=.0025),
          colorkey=list(at=seq(0,.0425,by=.0025)),
          col.regions=heat.colors(100)[length(heat.colors(100)):1])
dev.off()
