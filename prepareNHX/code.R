require(phytools)
library(ggtree)
library(treeio)
require(ggplot2)


taxa<-read.table("taxa.corresp.txt", sep=" ")
tr<-read.tree("tree.tre")
tr$tip.label<-taxa$V2[as.numeric(tr$tip.label)]


##continuous character
chias<-read.table("Chiasma_Loc_Data", row.names=1)
chias.small<-chias[tr$tip.label,]


chimean<-as.numeric(chias.small[,1])
names(chimean)<-rownames(chias.small)
fitchimean<-fastAnc(tr,chimean,vars=FALSE,CI=FALSE)
CHIMEAN<-c(chimean, fitchimean)

chimed<-as.numeric(chias.small[,2])
names(chimed)<-rownames(chias.small)
fitchimed<-fastAnc(tr,chimed,vars=FALSE,CI=FALSE)
CHIMED<-c(chimed, fitchimed)

chivar<-as.numeric(chias.small[,3])
names(chivar)<-rownames(chias.small)
fitchivar<-fastAnc(tr,chivar,vars=FALSE,CI=FALSE)
CHIVAR<-c(chivar, fitchivar)

##other characters
karyo<-read.table("Karyo_Data")
karyo2<-as.data.frame(do.call(rbind, strsplit(karyo$V2, "")))
rownames(karyo2)<-karyo$V1

karyo.small<-karyo2[tr$tip.label,]

chrom<-as.integer(karyo.small$V1)
names(chrom)<-rownames(karyo.small)
anc.karyo <- ace(chrom, tr, type="discrete", model="ER")
colnames(anc.karyo$lik.anc)<-c("X0","XX0","XXY","XXXY")
anc.karyo.final<-apply(anc.karyo$lik.anc, 1, function(x) names(which(x==max(x))))
KARYO<-c(colnames(anc.karyo$lik.anc)[chrom+1], anc.karyo.final)

BROWNIANM<-fastBM(tr, internal=TRUE)


DF<-data.frame(CHIMEAN, CHIMED, CHIVAR, KARYO, BROWNIANM)
DF$nodeid<-paste("xx",1:nrow(DF),"xx", sep="")
DF$brlen<-tr$edge.length[match(1:(Ntip(tr)+Nnode(tr)), tr$edge[,2])]
DF$brlen[Ntip(tr)+1]<-tr$root.edge
DF$cleannames<-gsub("\\)","-_",gsub("\\(","_-",gsub(" ","",rownames(DF))))
DF$PARENTALKARYO<-DF$KARYO[tr$edge[match(1:nrow(DF), tr$edge[,2]),1]]
DF$NeoY<-grepl("Y",DF$KARYO)-grepl("Y",DF$PARENTALKARYO) #-1 -> loss, +1 -> gain 0->nodiff

# CREATE NAMES FOR NHX format
DF$NAME<-paste(DF$cleannames,":",DF$brlen,"[&&NHX:","Karyo=",DF$KARYO,":MedChi=",DF$CHIMED,":MeanChi=",DF$CHIMEAN,":VarChi=",CHIVAR,":Brownian=",DF$BROWNIANM,":NeoY=",DF$NeoY,"]", sep="")


#TRANSFROM THE TREE 
TR<-tr
TR$tip.label<-DF$nodeid[1:Ntip(TR)]
TR$node.label<-DF$nodeid[(Ntip(TR)+1):(Ntip(TR)+Nnode(TR))]
TR$edge.length<-NULL
TR$root.edge<-NULL
newicktr<-write.tree(TR)

for (i in 1:nrow(DF)) newicktr<-gsub(DF$nodeid[i],DF$NAME[i], newicktr)

#WRITE TREE TO OUTPUT FILE
cat(newicktr, file="SpiderFinal.nhx")
