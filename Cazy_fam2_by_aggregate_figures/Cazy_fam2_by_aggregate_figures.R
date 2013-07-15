data<-read.csv(file.choose())

head(data)

library(ggplot2)


#nstall.packages("plyr")
library(plyr)

stats<-ddply(data, .(agg_frac,Cazy_fam2), summarise,
MEAN=mean(Abundance),
SE=sd(Abundance)/sqrt(length(Abundance))
)
stats$agg_frac<-factor(stats$agg_frac, levels=c("micro","SM","MM","LM","WS"))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)

X1<-ggplot(stats, aes(agg_frac,MEAN))+geom_point(stat="identity")+facet_wrap(~Cazy_fam2,scales="free")+geom_errorbar(limits,width=0)+theme_bw()+theme(aspect.ratio=1)

X1<-ggplot(stats, aes(agg_frac,MEAN))+geom_bar(stat="identity")+facet_wrap(~Cazy_fam2,scales="free")+geom_errorbar(limits,width=0)+theme_bw()+theme(aspect.ratio=1)


X1

ggplot(data, aes(log(Abundance)))+geom_histogram()
list(unique(data$Cazy_fam2))
summary(test<-aov(log(Abundance)~agg_frac,data=subset(data, Cazy_fam2=="PL")))
TukeyHSD(test)


stats.2<-ddply(data, .(agg_frac), summarise,
MEAN=mean(Abundance),
SE=sd(Abundance)/sqrt(length(Abundance))
)

X2<-ggplot(stats.2, aes(agg_frac,MEAN))+geom_point(stat="identity")+geom_errorbar(limits,width=0)+theme_bw()+theme(aspect.ratio=1)
stats.2$agg_frac<-factor(stats.2$agg_frac, levels=c("micro","SM","MM","LM","WS"))
X2