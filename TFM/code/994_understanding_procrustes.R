source("tools/load_libraries.R")
# http://www.flutterbys.com.au/stats/tut/tut15.1.html

# Data generation
set.seed(1)
x <- seq(0,50,l=10)
n <- 10
sp1<-coenocline(x=x,A0=5,m=0,r=2,a=1,g=1,int=T, noise=T)
sp2<-coenocline(x=x,A0=70,m=7,r=30,a=1,g=1,int=T, noise=T)
sp3<-coenocline(x=x,A0=50,m=15,r=30,a=1,g=1,int=T, noise=T)
sp4<-coenocline(x=x,A0=7,m=25,r=20,a=0.4,g=0.1,int=T, noise=T)
sp5<-coenocline(x=x,A0=40,m=30,r=30,a=0.6,g=0.5,int=T, noise=T)
sp6<-coenocline(x=x,A0=15,m=35,r=15,a=0.2,g=0.3,int=T, noise=T)
sp7<-coenocline(x=x,A0=20,m=45,r=25,a=0.5,g=0.9,int=T, noise=T)
sp8<-coenocline(x=x,A0=5,m=45,r=5,a=1,g=1,int=T, noise=T)
sp9<-coenocline(x=x,A0=20,m=45,r=15,a=1,g=1,int=T, noise=T)
sp10<-coenocline(x=x,A0=30,m=50,r=5,a=1,g=1,int=T, noise=T)
X <- cbind(sp1, sp10,sp9,sp2,sp3,sp8,sp4,sp5,sp7,sp6)
#X<-X[c(1,10,9,2,3,8,4,5,7,6),] 
colnames(X) <- paste("Sp",1:10,sep="")
rownames(X) <- paste("Site", c(1,10,9,2,3,8,4,5,7,6), sep="")
X <- X[c(1,4,5,7,8,10,9,6,3,2),]
data <- data.frame(Sites=factor(rownames(X),levels=rownames(X)), X)

data.dist <- vegdist(data[,-1], "bray")
