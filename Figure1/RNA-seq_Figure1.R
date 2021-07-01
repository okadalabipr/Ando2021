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

#Draw line graphs of mean fold change in expression of 5 DEG clusters in Ctrl
ERG1 <- read.table("../Data/Subcluster1_ERGs.txt",sep="\t",header= F)
ERG2 <- read.table("../Data/Subcluster2_ERGs.txt",sep="\t",header= F)
ERG <- c(as.character(ERG1[,1]),as.character(ERG2[,1]))
ERG <- sort(ERG)

IRG1 <- read.table("../Data/Subcluster1_IRGs.txt",sep="\t",header= F)
IRG2 <- read.table("../Data/Subcluster2_IRGs.txt",sep="\t",header= F)
IRG <- c(as.character(IRG1[,1]),as.character(IRG2[,1]))
IRG <- sort(IRG)
IRG[54] <- c("NKX3-1")

DRG1 <- read.table("../Data/Subcluster1_DRGs.txt",sep="\t",header= F)
DRG2 <- read.table("../Data/Subcluster2_DRGs.txt",sep="\t",header= F)
DRG <- c(as.character(DRG1[,1]),as.character(DRG2[,1]))
DRG <- sort(DRG)

Downregulated <- read.table("../Data/Downregulatedcluster_12022019.txt",sep="\t",header= T)
Downregulated <- Downregulated[,2]
Downregulated <- as.character(Downregulated)
Downregulated[93] <- c("SEPT5")

Others <- read.table("../Data/Otherscluster_12022019.txt",sep="\t",header= T)
Others <- Others[,2]
Others <- as.character(Others)

wt2 <- wt[which(rownames(wt) %in% DRG),]
data <- wt2
for(k in 1:dim(wt2)[2]){
  for(i in 1:dim(wt2)[1]){
    data[i,k] <- data[i,k]/wt2[i,1]
  }  
}
wtstandardfc <- data

kd2 <- kd[which(rownames(kd) %in% DRG),]
data <- kd2
for(k in 1:dim(kd2)[2]){
  for(i in 1:dim(kd2)[1]){
    data[i,k] <- data[i,k]/kd2[i,1]
  }  
}
kdstandardfc <- data

wilcox.exact(x=wtstandardfc[,1],y=kdstandardfc[,1],paired=F,alternative="g")
wilcox.exact(x=wtstandardfc[,2],y=kdstandardfc[,2],paired=F,alternative="g")
wilcox.exact(x=wtstandardfc[,3],y=kdstandardfc[,3],paired=F,alternative="g")
wilcox.exact(x=wtstandardfc[,4],y=kdstandardfc[,4],paired=F,alternative="g")
wilcox.exact(x=wtstandardfc[,5],y=kdstandardfc[,5],paired=F,alternative="g")
wilcox.exact(x=wtstandardfc[,6],y=kdstandardfc[,6],paired=F,alternative="g")
wilcox.exact(x=wtstandardfc[,7],y=kdstandardfc[,7],paired=F,alternative="g")
wilcox.exact(x=wtstandardfc[,8],y=kdstandardfc[,8],paired=F,alternative="g")
wilcox.exact(x=wtstandardfc[,9],y=kdstandardfc[,9],paired=F,alternative="g")
wilcox.exact(x=wtstandardfc[,10],y=kdstandardfc[,10],paired=F,alternative="g")
wilcox.exact(x=wtstandardfc[,11],y=kdstandardfc[,11],paired=F,alternative="g")
wilcox.exact(x=wtstandardfc[,12],y=kdstandardfc[,12],paired=F,alternative="g")
wilcox.exact(x=wtstandardfc[,13],y=kdstandardfc[,13],paired=F,alternative="g")

wtstandardfcmean <- apply(wtstandardfc,2,mean)
kdstandardfcmean <- apply(kdstandardfc,2,mean)
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
  #scale_y_continuous(breaks=seq(0.5,2.5,by=0.5),limits=c(0.5,2.5))+
  geom_line(lwd=3.0,aes(color=group))+
  scale_colour_manual(values = c(siCtrl="magenta",siIkBa="seagreen"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=20,face="italic"),panel.background=element_rect(color="black", size=2, fill="white" ),panel.grid.minor.y = element_blank(),panel.grid.major.x = element_line(colour="grey", size =1.0),panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),legend.position = "none",axis.title.x = element_text(size=15,margin=margin(t = 25, r = 0, b = 0, l = 0)),axis.title.y = element_text(size=16,margin=margin(t = 0, r = 25, b = 0, l = 0)),axis.text.x = element_text(size=15),axis.text.y = element_text(size=15));g

#Draw line graphs of mean fold change in expression of 3 TNF-induced clusters in Ctrl
#w1 <- read.table("Subcluster1_ERGs.txt",sep="\t",stringsAsFactors = FALSE)
#w2 <- read.table("Subcluster2_ERGs.txt",sep="\t",stringsAsFactors = FALSE)
#w1 <- read.table("Subcluster1_IRGs.txt",sep="\t",stringsAsFactors = FALSE)
#w2 <- read.table("Subcluster2_IRGs.txt",sep="\t",stringsAsFactors = FALSE)
#w1 <- read.table("Subcluster1_DRGs.txt",sep="\t",stringsAsFactors = FALSE)
#w2 <- read.table("Subcluster2_DRGs.txt",sep="\t",stringsAsFactors = FALSE)
w1 <- w1[,1]
w2 <- w2[,1]

#Draw a heatmap of the fold change in expression
wt2 <- wt[which(rownames(wt) %in% wtdeg),]
kd2 <- kd[which(rownames(kd) %in% wtdeg),]
wt2 <- wt2[which(rownames(wt2) %in% w1),]
data <- wt2
for(k in 1:dim(wt2)[2]){
  for(i in 1:dim(wt2)[1]){
    data[i,k] <- data[i,k]/wt2[i,1]
  }  
}
wtstandardfc <- data

kd2 <- kd2[which(rownames(kd2) %in% w1),]
data <- kd2
for(k in 1:dim(kd2)[2]){
  for(i in 1:dim(kd2)[1]){
    data[i,k] <- data[i,k]/kd2[i,1]
  }  
}
kdstandardfc <- data

cwt <- wtstandardfc
ckd <- kdstandardfc
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

wilcox.exact(x=wtstandardfcextract[,1],y=kdstandardfcextract[,1],paired=F,alternative="g")
wilcox.exact(x=wtstandardfcextract[,2],y=kdstandardfcextract[,2],paired=F,alternative="g")
wilcox.exact(x=wtstandardfcextract[,3],y=kdstandardfcextract[,3],paired=F,alternative="g")
wilcox.exact(x=wtstandardfcextract[,4],y=kdstandardfcextract[,4],paired=F,alternative="g")
wilcox.exact(x=wtstandardfcextract[,5],y=kdstandardfcextract[,5],paired=F,alternative="g")
wilcox.exact(x=wtstandardfcextract[,6],y=kdstandardfcextract[,6],paired=F,alternative="g")
wilcox.exact(x=wtstandardfcextract[,7],y=kdstandardfcextract[,7],paired=F,alternative="g")
wilcox.exact(x=wtstandardfcextract[,8],y=kdstandardfcextract[,8],paired=F,alternative="g")
wilcox.exact(x=wtstandardfcextract[,9],y=kdstandardfcextract[,9],paired=F,alternative="g")
wilcox.exact(x=wtstandardfcextract[,10],y=kdstandardfcextract[,10],paired=F,alternative="g")
wilcox.exact(x=wtstandardfcextract[,11],y=kdstandardfcextract[,11],paired=F,alternative="g")
wilcox.exact(x=wtstandardfcextract[,12],y=kdstandardfcextract[,12],paired=F,alternative="g")
wilcox.exact(x=wtstandardfcextract[,13],y=kdstandardfcextract[,13],paired=F,alternative="g")

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
  scale_colour_manual(values = c(siCtrl="magenta",siIkBa="seagreen"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=20,face="italic"),panel.background=element_rect(color="black", size=2, fill="white" ),panel.grid.minor.y = element_blank(),panel.grid.major.x = element_line(colour="grey", size =1.0),panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),legend.position = "none",axis.title.x = element_text(size=15,margin=margin(t = 25, r = 0, b = 0, l = 0)),axis.title.y = element_text(size=16,margin=margin(t = 0, r = 25, b = 0, l = 0)),axis.text.x = element_text(size=15),axis.text.y = element_text(size=15));g
