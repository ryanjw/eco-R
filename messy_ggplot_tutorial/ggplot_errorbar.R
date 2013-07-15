head(wgs[,1:10])

names(wgs)

summary(aov((k__Bacteria.p__Verrucomicrobia)~SoilFrac, data=wgs))
#install.packages("ggplot2")
library(ggplot2)

p<-ggplot(wgs, aes(SoilFrac, k__Bacteria.p__Verrucomicrobia))+geom_bar(stat="identity")


p+theme_bw()+theme(aspect.ratio=1)

#nstall.packages("plyr")
library(plyr)

stats<-ddply(wgs, .(SoilFrac), summarise,
MEAN=mean(k__Bacteria.p__Verrucomicrobia),
SE=sd(k__Bacteria.p__Verrucomicrobia)/sqrt(length(k__Bacteria.p__Verrucomicrobia))
)

attach(stats)
limits<-aes(ymin=MEAN-SE,ymax=MEAN+SE)

ggplot(stats, aes(SoilFrac, MEAN, colour=SoilFrac))+geom_bar(stat="identity", aes(fill=SoilFrac))+geom_errorbar(limits)+theme(aspect.ratio=1)+scale_fill_discrete(values=colors)

