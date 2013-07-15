wgs<-read.csv(file.choose())

library(vegan)

dim(wgs)
head(wgs[,1:5])
wgs.comm<-decostand(wgs[,3:305],"total")
wgs.nmds<-metaMDS(wgs.comm,autotransform=FALSE,k=3)

wgs.corr<-envfit(wgs.nmds, wgs.comm, permutations=9999)
str(wgs.corr)
str(wgs.corr$vectors$arrows)
str(wgs.corr$vectors$pvals)

X<-data.frame(wgs.corr$vectors$arrows,wgs.corr$vectors$r,wgs.corr$vectors$pvals)
head(X)
names(X)<-c("NMDS1","NMDS2","R.2","Pval")



X$P.FDR<-p.adjust(X$Pval,"fdr",n=length(X$Pval))

head(X)

X.sig<-subset(X, P.FDR <0.05)

write.csv(X.sig, "env_fit_OTU.csv")


SoilFrac<-wgs$SoilFrac
plot(wgs.nmds, type="n", xlim=c(-.25,.35),ylim=c(-.3,.2))
plot(wgs.corr, p.max=0.0024, labels=NULL)
ordiellipse(wgs.nmds, SoilFrac, show.groups="Micro",lwd=2)
ordiellipse(wgs.nmds, SoilFrac, show.groups="SM", col="purple",lwd=2)
ordiellipse(wgs.nmds, SoilFrac, show.groups="MM", col="orange",lwd=2)
ordiellipse(wgs.nmds, SoilFrac, show.groups="LM",col="red",lwd=2)
ordiellipse(wgs.nmds, SoilFrac, show.groups="WS",col="grey",lty=2,lwd=2)
legend.text<-c("Micro","SM","MM","LM","WS")
colors<-c("black","purple","orange","red","grey")
legend(-.25,-.15,legend.text,col=colors,pch=16)