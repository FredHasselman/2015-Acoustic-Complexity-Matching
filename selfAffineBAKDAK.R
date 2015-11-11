require("devtools")
source_url("https://raw.githubusercontent.com/FredHasselman/scicuRe/master/scicuRe_source.R")
inIT(c("wmtsa","plyr","lattice","tuneR","EMD"))

#setwd("~/Dropbox/Hasselman2014-PeerJ-Classifying_Acoustic_Signals/STIMULI/") #ARTICLE_TEX/FIGURES")
pre <- "~/Library/Mobile Documents/com~apple~CloudDocs/Projects/Classifying Complex Signals/Hasselman2014-PeerJ-Classifying_Acoustic_Signals"

setwd(paste0(pre,"/STIMULI"))

normal<- paste0("BAKDAK",1:10,".WAV")
slowed <- paste0("BAKDAKV",1:10,".WAV")
amped <- paste0("BAKDAKA",1:10,".WAV")
both <- paste0("BAKDAKB",1:10,".WAV")

stims <- list(normal,slowed,amped,both)
S1  <- llply(stims,function(s) llply(s,readWave))
S1d <- llply(S1,function(s) llply(s,downsample,samp.rate=2000))

tsy    <- list()
W      <- list()
W.tree <- list()
s.n <- 256
cnt <- 0
for(w in 1:length(S1d)){
  for(s in 1:length(S1d[[w]])){
    cnt <- cnt+1
    cat(cnt)
  tt <- cbind((1:length(S1d[[w]][[s]]@left))*(1/S1d[[w]][[s]]@samp.rate))
  y  <- hilbertspec(xt=cbind(S1d[[w]][[s]]@left),tt=tt)
  tsy[[cnt]] <- ts(data=abs(y$amplitude)/max(abs(y$amplitude)),start=0,frequency=S1d[[w]][[s]]@samp.rate)
  s.rng    <- deltat(tsy[[cnt]]) * c(1, length(tsy[[cnt]]))
  W[[cnt]]      <- wavCWT(tsy[[cnt]], wavelet="gaussian2",scale.range=s.rng,n.scale=s.n)
  W.tree[[cnt]] <- wavCWTTree(W[[cnt]],type='extrema') 
  }
}
names(tsy) <- c(normal,slowed,amped,both)

setwd(paste0(pre,"/DATA"))
save.image(file="vd.RData")

load("vd.RData")

w=1

for(w in 1){
  for(s in 1:length(S1[[w]])){
    cnt <- cnt+1
    cat(cnt)
  tt <- cbind((1:length(S1[[w]][[s]]@left))*(1/S1[[w]][[s]]@samp.rate))
  y  <- hilbertspec(xt=cbind(S1[[w]][[s]]@left),tt=tt)
  tsy[[cnt]] <- ts(data=abs(y$amplitude)/max(abs(y$amplitude)),start=0,frequency=S1[[w]][[s]]@samp.rate)
  s.rng    <- deltat(tsy[[cnt]]) * c(1, length(tsy[[cnt]]))
  W[[cnt]]      <- wavCWT(tsy[[cnt]], wavelet="gaussian2",scale.range=s.rng,n.scale=s.n)
  W.tree[[cnt]] <- wavCWTTree(W[[cnt]],type='extrema') 
  }
}


#W[[w]]      <- wavCWT(tsy[[w]], wavelet="gaussian2",scale.range=s.rng,n.scale=s.n)
W.tree[[w]] <- wavCWTTree(W[[w]],type='maxima') 

holder=holderSpectrum(W.tree[[w]])
print(holder)


setwd(paste0(pre,"/ARTICLE_TEX/FIGURES"))
for(w in c(1)){
s.rng    <- deltat(tsy[[w]]) * c(1, length(tsy[[w]]))
pdf(paste0("MFno",w,".pdf"),paper="a4r",width=0,height=0)
plot(W[[w]],grid.size=2^6, col = palette(gray(seq(0.1,.9,len=s.n))),ylim=round(log2(s.rng)+c(2,0.5)),yaxt="n",bty="n",ylab="Temporal Scale (sec.)",xlab="Time (sec.)") #,main=names(tsy)[w])
axis(2,at=seq(-9,log2(.5),length=3),labels=paste(signif(seq(2^-9,2^-.7369,length=3),digits=2)))
plot(W.tree[[w]],pch=20,add=T,ylim=round(log2(s.rng)+c(-.5,.5)),yaxt="n",bty="n")
lines(tsy[[w]]+(log2(s.rng)[1])+3, lwd=2,type="l",col="gray",yaxt="n",bty="n")
dev.off()
}

for(w in c(11,14,15,16,20,21,24,25,26,30,31,34,35,36,40)){
s.rng    <- deltat(tsy[[w]]) * c(1, length(tsy[[w]]))
pdf(paste0("MFno",w,".pdf"),paper="a4r",width=0,height=0)
plot(W[[w]],grid.size=2^10, col = palette(gray(seq(0.1,.9,len=s.n))),ylim=round(log2(s.rng)+c(2,0.5)),ylab="",xlab="Time (sec.)",yaxt="n",bty="n")
     #ylab="Scale (sec.)",main=names(tsy)[w])
#axis(2,at=seq(-10,log2(.5),length=3),labels=paste(signif(seq(2^-10,2^-.7369,length=3),digits=2)))
plot(W.tree[[w]],pch=20,add=T,yaxt="n",bty="n")
lines(tsy[[w]]+(log2(s.rng)[1])+3, lwd=2,type="l",col="gray",yaxt="n",bty="n")
dev.off()
}


wm<-as.matrix(W[[w]])
image(wm)

require(ggplot2)
require(reshape2)
df <- melt(wm)
qplot(Var1, Var2, data = df, fill = value, geom = "raster")

cwtPlot <- function(W,W.tree,tsy,s.n=s.n,s.rng=s.rng,i=c,j=r,col){
   plot(W,grid.size=2^6, col=col)
   plot(W.tree,pch=20,add=T,ylim=round(log2(s.rng)+c(-.5,.5)),yaxt="n",bty="n")
#   popViewport()
   lines(tsy+(log2(s.rng)[1]+3), lwd=2,type="l",col="gray",yaxt="n",bty="n")
}

#palette(gray(seq(0.1,.9,len=s.n))))
#ylim=round(log2(s.rng)+c(2,0)),yaxt="n",bty="n",ylab="Scale (sec.)",xlab="Time (sec.)")
# cols=10
# numPlots=40
# layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),ncol = cols, nrow = ceiling(numPlots/cols))
# grid.show.layout(grid.layout(nrow(layout), ncol(layout)))

# grid.newpage()
# grid.rect(gp=gpar(fill="white"))
# pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
w=0
    ## set the plot grid 
old.plt <- splitplot(4,10,1,new=T)
par(mar = rep(1, 4))
for(r in 1:4){
  for(c in 1:10){
     w=w+1
    if (w>1) splitplot(4,10,w)
    wm <- as.matrix(W[[w]])
    image(wm, col = palette(gray(seq(0.1,.9,len=50))),axes=FALSE)
    #plot(W[[w]],grid.size=2^6, col = palette(gray(seq(0.1,.9,len=s.n))),ylim=round(log2(s.rng)+c(2,0)),yaxt="n",bty="n",ylab="Scale (sec.)",main=names(tsy)[w])
    #cwtPlot(W[[w]],W.tree[[w]],tsy[[w]],s.n,s.rng,c,r,palette(gray(seq(0.1,.9,len=s.n))))
  }
}
 
par(old.plt)
  
  
# if(w%in%c(1,11,21,31)) rcnt<-rcnt+1
# print((p+t+l), vp = viewport(layout.pos.row = rcnt, layout.pos.col = w))
}



grid.newpage()
grid.rect(gp=gpar(fill="grey"))
pushViewport(viewport(layout=grid.layout(2, 2)))
clip.demo <- function(i, j, clip1, clip2) {
  pushViewport(viewport(layout.pos.col=i,
                        layout.pos.row=j))
  pushViewport(viewport(width=0.6, height=0.6, clip=clip1))
  grid.rect(gp=gpar(fill="white"))
  grid.circle(r=0.55, gp=gpar(col="red", fill="pink"))
  popViewport()
  pushViewport(viewport(width=0.6, height=0.6, clip=clip2))
  grid.polygon(x=c(0.5, 1.1, 0.6, 1.1, 0.5, -0.1, 0.4, -0.1),
               y=c(0.6, 1.1, 0.5, -0.1, 0.4, -0.1, 0.5, 1.1),
               gp=gpar(col="blue", fill="light blue"))
  popViewport(2)
}

clip.demo(1, 1, FALSE, FALSE)
clip.demo(1, 2, TRUE, FALSE)
clip.demo(2, 1, FALSE, TRUE)
clip.demo(2, 2, TRUE, TRUE)
popViewport()

dev.off()
}
plot.new()



# MULTIPLOT FUNCTION ------------------------------------------------------------------------------------------------------------------
#
# [copied from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/ ]
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multi.PLOT <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}




x <- make.signal("linchirp", n=1024)
x.dwt <- wavDWT(tsy[[w]], n.levels = 256)
x.modwt <- wavMODWT(tsy[[w]], n.levels = 256)

## calculate the wavelet details for all crystals 
## of the DWT and MODWT 
wavMRD(x.dwt)
wavMRD(x.modwt)

## plot the wavelet details for levels 1 and 3 of 
## the MODWT 
plot(wavMRD(x.modwt, level = c(1,3)))

## plot the wavelet details for all levels of the 
## MODWT of a linear chirp. 
plot(wavMRD(x.modwt))



#par(mfrow=c(1,1))

