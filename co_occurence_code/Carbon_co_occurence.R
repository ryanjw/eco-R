carbon<-read.csv(file.choose())

(dim(carbon)-4)*5
head(carbon[,1:4])
carbon<-subset(carbon, Crop=="PF")
aggregate<-as.vector((unique((carbon$SoilFrac))))
length(aggregate)


aggregate.results<-matrix(nrow=0,ncol=7)
options(warn=-10)

for(a in 1:length(aggregate)){
	
	aggregate.temp<-aggregate[a]
	temp<-subset(carbon, aggregate==aggregate.temp)
	
	for(b in 4:(dim(temp)[2]-1)){
		
		for(c in (b+1):(dim(temp)[2])){
			
			enz1.ab<-sum(temp[,b])
			enz2.ab<-sum(temp[,c])
			
			if(enz1.ab >0 & enz2.ab >0){
				test<-cor.test(temp[,b],temp[,c],method="spearman",na.action=na.rm)
				rho<-test$estimate
				p.value<-test$p.value
			}
			
			if(enz1.ab ==0 | enz2.ab ==0){
				rho<-0
				p.value<-1
			}	
			
			new.row<-c(aggregate[a],names(temp)[b],names(temp)[c],rho,p.value,enz1.ab,enz2.ab)
			aggregate.results<-rbind(aggregate.results,new.row)			
			
		}
		print((a*315+b-315)/1655)
	}
	
	
	
	
}


aggregate.results<-data.frame(aggregate.results)
aggregate.results$P.FDR<-p.adjust(aggregate.results$P.value,"fdr",n=length(aggregate.results$P.value))
head(aggregate.results)
str(aggregate.results)
names(aggregate.results)<-c("SoilFrac","Order1","Order2","Corr","P.value","Enz1.Abundance","Enz2.Abundance")
write.csv(aggregate.results, "Aggregate_Carbon_Data_Co_Occurence_Results.csv")


####nothing was significant based on correlation and low sample size###

aggregate.results.2<-matrix(nrow=0,ncol=7)
options(warn=-10)


for(b in 4:(dim(carbon)[2]-1)){
		
		for(c in (b+1):(dim(carbon)[2])){
			
			enz1.ab<-sum(carbon[,b])
			enz2.ab<-sum(carbon[,c])
			
			if(enz1.ab >0 & enz2.ab >0){
				test<-cor.test(carbon[,b],carbon[,c],method="spearman",na.action=na.rm)
				rho<-test$estimate
				p.value<-test$p.value
			}
			
			if(enz1.ab ==0 | enz2.ab ==0){
				rho<-0
				p.value<-1
			}	
			
			new.row<-c(aggregate[a],names(carbon)[b],names(carbon)[c],rho,p.value,enz1.ab,enz2.ab)
			aggregate.results.2<-rbind(aggregate.results.2,new.row)			
			
		}
		print((a*315+b-315)/315)
	}
	
aggregate.results.2$P.value
aggregate.results.2<-data.frame(aggregate.results.2)
P.FDR<-p.adjust(as.matrix(aggregate.results.2$P.value),"fdr",n=length(aggregate.results.2$P.value))
aggregate.results.2$P.FDR<-P.FDR
head(aggregate.results)
str(aggregate.results)
names(aggregate.results.2)<-c("SoilFrac","Order1","Order2","Corr","P.value","Enz1.Abundance","Enz2.Abundance","FDR.P")
write.csv(aggregate.results.2, "Aggregate_Carbon_Data_Co_Occurence_Results_no_agg_specific.csv")

library(igraph)
head(aggregate.results.2)
ag.results<-read.csv(file.choose())
ag.sig<-subset(aggregate.results.2, FDR.P < 0.05 )
str(ag.sig)
head(ag.sig)
carbon.graph<-simplify(graph.edgelist(as.matrix(ag.sig[,c(2,3)]),directed=FALSE))
V(carbon.graph)$label<-NA

carbon.network<-data.frame(degree(carbon.graph),betweenness(carbon.graph))
head(carbon.network)
write.csv(carbon.network, "carbon_network_stats.csv")
	

