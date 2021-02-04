library(DESeq2)
library(edgeR)
library(VennDiagram)
library(genefilter)
library(e1071)

#Venn diagram of DEGs in WT and siIkBa
grid.newpage()
draw.pairwise.venn(371,922,321,category = c("", ""), fill = c("tomato", "limegreen"), cex=1, alpha = rep(0.6,2),cat.pos = c(50,-100),  cat.dist = rep(0.05, 2), scaled = TRUE, rotation.degree = 180)

#TMM normalization of expression using edgeR 
count <- read.table("../Data/counts_all.txt",sep="\t",header= T,row.names = 1,skip = 1)
#count <- read.table("../Data/counts_KO_all.txt",sep="\t",header= T,row.names = 1,skip = 1)
count1=count[,5:31]
count1 <- as.matrix(count1)
count <- apply(count1, 2, as.numeric)
for(k in 2:27){
  for(i in 1:dim(count)[1]){
    count[i,k] <- (count[i,k]/count[i,1])*1000
  }
}
group <- colnames(count[,2:27]) 
d <- DGEList(counts = count[,2:27], group = group)
d <- calcNormFactors(d,method="RLE") 
dn <- d$samples
dn$ln <- dn$lib.size*dn$norm.factors
count <- read.table("../Data/counts_all.txt",sep="\t",header= T,row.names = 1,skip = 1)
#count <- read.table("../Data/counts_KO_all.txt",sep="\t",header= T,row.names = 1,skip = 1)
count1=count[,5:31]
count1 <- as.matrix(count1)
count <- apply(count1, 2, as.numeric)
for(k in 2:27){
  for(i in 1:dim(count)[1]){
    count[i,k] <- (count[i,k]/count[i,1])*1000
  }
}
for(k in 1:26){
  count[,k+1] <- (count[,k+1]/dn$ln[k])*1000000
}
countname <- read.table("../Data/counts_all.txt",sep="\t",header= T,row.names = 1,skip = 1)
#countname <- read.table("../Data/counts_KO_all.txt",sep="\t",header= T,row.names = 1,skip = 1)
rownames(count) <- rownames(countname)
tpm <- matrix(0,nrow = dim(count)[1],ncol = dim(count)[2]/2)
count <- count[,2:27]
for(i in 1:13){
  tpm[,i] <- (count[,2*i-1]+count[,2*i])/2 
}
tpm <- as.data.frame(tpm)
rownames(tpm) <- rownames(count) 
colnames(tpm) <- seq(0,180,15)
wt <- tpm
#kd <- tpm
wtdeg <- read.table("../Data/WT_DEGs.txt",skip=1,row.names=1,sep="\t",stringsAsFactors = FALSE)
kddeg <- read.table("../Data/siIkBa_DEGs.txt",skip=1,row.names=1,sep="\t",stringsAsFactors = FALSE)
wtdeg <- wtdeg[,1]
kddeg <- kddeg[,1]
wtdeg <- as.character(wtdeg)
kddeg <- as.character(kddeg)
wt2 <- wt[which(rownames(wt) %in% wtdeg),]
kd2 <- kd[which(rownames(kd) %in% kddeg),]

#Cluster DEGs into 5 clusters
d <- wt2
d <- genescale(d,axis=1,method="Z")
set.seed(3)
f1 <- cmeans(d,5,2,verbose=FALSE,m=1.3,dist="euclidean",method="cmeans")
q1 <- as.data.frame(f1[3])
a1 <- rownames(q1)[which(f1$cluster=="1")]
a2 <- rownames(q1)[which(f1$cluster=="2")]
a3 <- rownames(q1)[which(f1$cluster=="3")]
a4 <- rownames(q1)[which(f1$cluster=="4")]
a5 <- rownames(q1)[which(f1$cluster=="5")]
length(a1)
length(a2)
length(a3)
length(a4)
length(a5)

erg <- wt2[which(rownames(wt2) %in% a1),]
irg <- wt2[which(rownames(wt2) %in% a2),]
drg <- wt2[which(rownames(wt2) %in% a4),]

d <- kd2[-which(rownames(kd2) %in% wtdeg),]
d <- genescale(d,axis=1,method="Z")
set.seed(11)
f1 <- cmeans(d,3,2,verbose=FALSE,m=1.3,dist="euclidean",method="cmeans")
q1 <- as.data.frame(f1[3])
b1 <- rownames(q1)[which(f1$cluster=="1")]
b2 <- rownames(q1)[which(f1$cluster=="2")]
b3 <- rownames(q1)[which(f1$cluster=="3")]
length(b1)
length(b2)
length(b3)

#Subclustering of ERGs, IRGs and DRGs
setwt <- wtstandardfc[which(rownames(wtstandardfc) %in% rownames(erg)),]
setkd <- kdstandardfc[which(rownames(kdstandardfc) %in% rownames(erg)),]
deg <- cbind(setwt[,1:6],setkd[,1:6]) #ERGs
#deg <- cbind(setwt[,1:10],setkd[,1:10]) #IRGs
#deg <- cbind(setwt[,1:12],setkd[,1:12]) #DRGs
deg <- as.data.frame(deg)
zscore <- genescale(deg,axis=1,method="Z")
set.seed(2) #ERGs
#set.seed(3) #IRGs
#set.seed(1) #DRGs
f1 <- cmeans(zscore,2,2,verbose=FALSE,m=1.3,dist="euclidean",method="cmeans")
q1 <- as.data.frame(f1[3])
w1 <- rownames(q1)[which(f1$cluster=="1")]
w2 <- rownames(q1)[which(f1$cluster=="2")]

write.table(w1,"Subcluster1_ERGs.txt",col.names=T,sep="\t")
write.table(w2,"Subcluster2_ERGs.txt",col.names=T,sep="\t")
w1 <- read.table("Subcluster1_ERGs.txt",sep="\t",stringsAsFactors = FALSE)
w2 <- read.table("Subcluster2_ERGs.txt",sep="\t",stringsAsFactors = FALSE)
w1 <- w1[,1]
w2 <- w2[,1]

#Draw a heatmap of the fold change in expression
data <- wt
for(k in 1:dim(wt)[2]){
  for(i in 1:dim(wt)[1]){
    data[i,k] <- data[i,k]/wt[i,1]
  }  
}
wtstandardfc <- data

data <- kd
for(k in 1:dim(kd)[2]){
  for(i in 1:dim(kd)[1]){
    data[i,k] <- data[i,k]/kd[i,1]
  }  
}
kdstandardfc <- data

cwt <- wtstandardfc[which(rownames(wtstandardfc) %in% b2),]
ckd <- kdstandardfc[which(rownames(kdstandardfc) %in% b2),]
df1 <- cbind(cwt,ckd)
df1 <- as.data.frame(df1)
df1 <- genescale(df1,axis=1,method="Z")

df <- as.data.frame(df1)
colnames(df) <- c(paste("wt_",seq(0,180,15),sep=""),paste("kd_",seq(0,180,15),sep=""))
df <- t(df)
df <- as.data.frame(df)
df$t <- c(paste("wt_",seq(0,180,15),sep=""),paste("kd_",seq(0,180,15),sep=""))
p <- colnames(df)
t <- c(paste("wt_",seq(0,180,15),sep=""),paste("kd_",seq(0,180,15),sep=""))
temp <- melt(df,
             id="t",
             measure=p 
)
temp <- temp[1:(dim(temp)[1]-26),]
temp$value <- as.numeric(temp$value)
q <- rep("",(dim(df1)[2]))
p <- colnames(df1)

g<-ggplot(temp,aes(as.factor(t),as.factor(variable)))+
  geom_tile(aes(fill=value))+
  scale_x_discrete(limits=c(paste("wt_",seq(0,180,15),sep=""),paste("kd_",seq(0,180,15),sep="")),expand=c(0.03,0))+
  scale_fill_gradient(low="blue",high="yellow");g

#Draw a line graph of the mean fold change in expression
wtstandardfc <- wtstandardfc[-which(rownames(wtstandardfc)=="PCSK5"),]
kdstandardfc <- kdstandardfc[-which(rownames(kdstandardfc)=="PCSK5"),]
wtstandardfcextract <- wtstandardfc[which(rownames(wtstandardfc) %in% w2),]
kdstandardfcextract<<- kdstandardfc[which(rownames(kdstandardfc) %in% w2),]
wtstandardfcmean <- apply(wtstandardfcextract,2,mean)
kdstandardfcmean <- apply(kdstandardfcextract,2,mean)
wtstandardfcmean <- as.data.frame(wtstandardfcmean)
kdstandardfcmean <- as.data.frame(kdstandardfcmean)
wtstandardfcmean <- t(wtstandardfcmean)
kdstandardfcmean <- t(kdstandardfcmean)
wtstandardfcmean <- as.data.frame(wtstandardfcmean)
kdstandardfcmean <- as.data.frame(kdstandardfcmean)

colnames(wtstandardfcmean) <- c("0","15","30","45","60","75","90","105","120","135","150","165","180")
wtstandardfcmean$t <- seq(1,dim(wtstandardfcmean)[1],1)
p <- colnames(wtstandardfcmean)
temp <- reshape2::melt(wtstandardfcmean,
             id="t",
             measure=p 
)
temp <- temp[1:(dim(temp)[1]-1),]
temp$group <- rep("siCtrl",dim(temp)[1])
wtstandardfcmeanextract <- temp

colnames(kdstandardfcmean) <- c("0","15","30","45","60","75","90","105","120","135","150","165","180")
kdstandardfcmean$t <- seq(1,dim(kdstandardfcmean)[1],1)
p <- colnames(kdstandardfcmean)
temp <- reshape2::melt(kdstandardfcmean,
             id="t",
             measure=p 
)
temp <- temp[1:(dim(temp)[1]-1),]
temp$group <- rep("siIkBa",dim(temp)[1])
kdstandardfcmeanextract <- temp
wtstandardfcmeanextract$value <- as.numeric(wtstandardfcmeanextract$value)
kdstandardfcmeanextract$value <- as.numeric(kdstandardfcmeanextract$value)
df <- rbind(wtstandardfcmeanextract,kdstandardfcmeanextract)
df$value <- as.numeric(df$value)
df$variable <- as.character(df$variable)

g <- ggplot(df,aes(x=variable,y=value,group=group))+
  theme_bw(base_size = 16)+
  theme(panel.grid = element_blank(),legend.position = "true")+
  labs(x="",y="",title="",colour="")+
  scale_x_discrete(breaks=c("0","15","30","45","60","75","90","105","120","135","150","165","180"),labels=c("0","","","","","","90","","","","","","180"),limits=c("0","15","30","45","60","75","90","105","120","135","150","165","180"))+
  scale_y_continuous(breaks=seq(0.5,3.5,by=0.5),limits=c(0.5,3.5))+
  geom_line(lwd=3.0,aes(color=group))+
  scale_colour_manual(values = c(siCtrl="#DC143C",siIkBa="seagreen"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=20,face="italic"),panel.background=element_rect(color="black", size=2, fill="white" ),panel.grid.minor.y = element_blank(),panel.grid.major.x = element_line(colour="grey", size =1.0),panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),legend.position = "none",axis.title.x = element_text(size=15,margin=margin(t = 25, r = 0, b = 0, l = 0)),axis.title.y = element_text(size=16,margin=margin(t = 0, r = 25, b = 0, l = 0)),axis.text.x = element_text(size=15),axis.text.y = element_text(size=15));g
