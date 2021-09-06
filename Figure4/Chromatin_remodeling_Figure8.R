library(genefilter)
library(reshape2)
library(ggplot2)
library(rtracklayer)
library(gridExtra)
library(edgeR)
library(scales)
library(DEoptim)
library(scales)

#Draw a heatmap of data
names <- c("GADD45A","IRS2","JUNB","KLF10","NUAK2","PPP1R15A","PTGER4","RND1","SPSB1","TNFAIP3","ZC3H12A")
refall <- c()
maxs <- c()
for(r in 1:length(names)){
  d <- names[r]
  
  ref <- wt2[which(rownames(wt2)==d),]
  
  data <- ref
  zero <- data[1,1]
  for(i in 1:13){
    data[1,i] <- data[1,i]/zero
  }
  
  refwtfc <- data.frame(x=seq(0,180,15),y=as.numeric(data[1,]))
  
  
  ref <- kd2[which(rownames(kd2)==d),]
  
  data <- ref
  zero <- data[1,1]
  for(i in 1:13){
    data[1,i] <- data[1,i]/zero
  }
  refkdfc <- data.frame(x=seq(0,180,15),y=as.numeric(data[1,]))
  
  df <- rbind(refwtfc,refkdfc)
  df <- as.data.frame(df)
  
  refwtfc <- rescale(refwtfc[,2], to=c(0,1))
  refkdfc <- rescale(refkdfc[,2], to=c(0,1))
  refwtfc <- data.frame(x=seq(0,180,15),y=refwtfc[1:13])
  refkdfc <- data.frame(x=seq(0,180,15),y=refkdfc[1:13])
  
  refwt <- refwtfc[,2]
  refkd <- refkdfc[,2]
  
  df <- c(refwt,refkd)
  
  refall <- c(refall, c(df))
}
refall <- data.frame(matrix(unlist(refall),ncol=26,byrow=T))  
rownames(refall) <- names
colnames(refall) <- rep(seq(0,180,15),2)
refall <- t(refall)
refall <- as.data.frame(refall)
colnames(refall) <- names
refall <- refall[1:13,]   #siCtrl
#refall <- refall[14:26,]   #siIkBa

p <- colnames(refall)
t <- colnames(refall)
refall$t <- seq(0,180,15)
temp <- reshape2::melt(refall,
                       id="t",
                       measure=p 
)
temp$variable <- as.character(temp$variable)
temp$value <- as.numeric(temp$value)

g <- ggplot(temp,aes(as.factor(t),as.factor(variable)))+
  geom_tile(aes(fill=value),color="grey40", stat="identity")+
  scale_fill_gradient(low="white",high="magenta",limits=c(0,1.0))+
  #scale_fill_gradient(low="white",high="mediumseagreen",limits=c(0,1.0))+
  labs(x = "",y = "") + 
  theme(panel.grid = element_blank())+
  scale_x_discrete(labels=c("0","15","30","45","60","75","90","105","120","135","150","165","180"),breaks=c("0","15","30","45","60","75","90","105","120","135","150","165","180"),expand = c(0, 0))+ 
  scale_y_discrete(breaks=sort(names,decreasing = T),labels=waiver(),limits=sort(names,decreasing = T),expand = c(0, 0))+
  theme(legend.position = "bottom",axis.ticks = element_blank(), axis.text.x = element_text(size = 7.0, angle = 0, hjust = 0.5, colour = "black"),axis.text.y=element_text(size = 5.6, angle = 0, hjust = 1, colour = "black",face="italic"));g

#Decreased time course chromatin accessibility at least at 30 min or 75 min
data <- read.table("Rawcounts_HOMER_siCtrl_120_DRG.txt",sep="\t",stringsAsFactors = FALSE)
data <- data[order(data[,2]),]
rownames(data) <- seq(1,dim(data)[1],1)

count1 <- data
d <- DGEList(counts = data[,4:dim(data)[2]])
d <- calcNormFactors(d,method="TMM") 
ls <- d$samples
librarysize <- ls$lib.size
normfactors <- ls$norm.factors
librarysize <- as.data.frame(librarysize)
normfactors <- as.data.frame(normfactors)
ln <- cbind(librarysize,normfactors)
ln$multiple <- ln$librarysize*ln$normfactors

data <- read.table("Rawcounts_HOMER_siCtrl_120_DRG.txt",sep="\t",stringsAsFactors = FALSE)
data <- data[order(data[,2]),]
rownames(data) <- seq(1,dim(data)[1],1)

data <- na.omit(data)
for(k in 1:dim(data)[1]){
  for(i in 1:(dim(data)[2]-3)){
    data[k,i+3] <- (data[k,i+3]/ln$multiple[i])*1000000
  }
}

names <- read.table("HOMER_siCtrl_120_DRG.txt",sep="\t",stringsAsFactors = FALSE,header=TRUE,quote="", comment.char="",)
names <- names[order(names$Start),]
rownames(names) <- seq(1,dim(names)[1],1)
names <- names[,2:4]
names <- makeGRangesFromDataFrame(names)
promoter <- read.table("Promoter_up500down500TSS(’uŠ·Ï).txt",sep="\t")
cluster <- read.table("Subcluster2_DRGs.txt",sep="\t")
cluster[,1] <- as.character(cluster[,1])
#cluster[27,1] <- c("NKX3-1")
promoter <- promoter[which(promoter$name %in% cluster[,1]),]
rownames(promoter) <- seq(1,dim(promoter)[1],1)
promoter2 <- promoter[,1:3]
colnames(promoter2) <- c("chr","start","end")
promoter2 <- makeGRangesFromDataFrame(promoter2)
OV <- findOverlaps(promoter2,names,type="any")
OV <- as.data.frame(OV)
promoter$name <- as.character(promoter$name)
target <- c()
for(i in 1:dim(OV)[1]){
  target[i] <- promoter[which(rownames(promoter) %in% OV[i,1]),5]
}


data <- data[which(rownames(data) %in% OV[,2]),]
data <- as.data.frame(data)
rownames(data) <- seq(1,dim(data)[1],1)

data$name <- target
data <- data[order(data$name),]
#data <- data[-which(data$name=="NFKBIA"),]
target <- data$name
rownames(data) <- seq(1,dim(data)[1],1)
data <- as.data.frame(data)
data <- data[,4:11]
data <- as.data.frame(data)

#data <- data[,1:4]  #change
data <- data[,5:8]  #change
data <- genescale(data,axis=1,method="Z")
data <- as.data.frame(data)
colnames(data) <-  c("0","30","75","120")
data <- as.data.frame(data) #Here
data$name <- target
rownames(data) <- seq(1,dim(data)[1],1)

thirty <- c()
seventy <- c()
for(i in 1:dim(data)[1]){
  if(as.character(data[i,2]>data[i,3])=="TRUE"){
    thirty[i] <- rownames(data)[i]
  }else{
    if(as.character(data[i,3]>data[i,4])=="TRUE"){
      seventy[i] <- rownames(data)[i]
    }else{
      next
    }
  }
}

thirty <- na.omit(thirty)
thirty <- as.numeric(thirty)
seventy <- na.omit(seventy)
seventy <- as.numeric(seventy)
thirty_data <- data[which(rownames(data) %in% thirty),]
seventy_data <- data[which(rownames(data) %in% seventy),]

#thirty_wt <- thirty_data  #change
thirty_kd <- thirty_data   #change

#seventy_wt <- seventy_data  #change
seventy_kd <- seventy_data   #change

#peaks <- c(rownames(thirty_wt), rownames(seventy_wt))  #change
peaks <- c(rownames(thirty_kd), rownames(seventy_kd))   #change

peaks <- as.numeric(peaks)
peaks_names <- data[which(rownames(data) %in% peaks),]

#peaks_names_wt <- peaks_names  #change
peaks_names_kd <- peaks_names   #change

peaks_names_ov_wt <- peaks_names_wt[which(peaks_names_wt$name %in% peaks_names_kd$name),]  #chang
#peaks_names_ov_kd <- peaks_names_wt[which(peaks_names_wt$name %in% peaks_names_kd$name),]  #change#

thirty_wt <- thirty_wt[which(thirty_wt$name %in% peaks_names_ov_wt$name),]  #change
thirty_kd <- thirty_kd[which(thirty_kd$name %in% peaks_names_ov_wt$name),]  #change
seventy_wt <- seventy_wt[which(seventy_wt$name %in% peaks_names_ov_wt$name),]  #change
seventy_kd <- seventy_kd[which(seventy_kd$name %in% peaks_names_ov_wt$name),]  #change

#thirty_wt <- thirty_wt[which(thirty_wt$name %in% peaks_names_ov_kd$name),]  #change
#thirty_kd <- thirty_kd[which(thirty_kd$name %in% peaks_names_ov_kd$name),]  #change
#seventy_wt <- seventy_wt[which(seventy_wt$name %in% peaks_names_ov_kd$name),]  #change
#seventy_kd <- seventy_kd[which(seventy_kd$name %in% peaks_names_ov_kd$name),]  #change

we_thirty_wt <- data.frame(high=thirty_wt[,2],low=thirty_wt[,3])
we_seventy_wt <- data.frame(high=seventy_wt[,3],low=seventy_wt[,4])
we_wt <- rbind(we_thirty_wt,we_seventy_wt)
we_wt <- as.data.frame(we_wt)

we_thirty_kd <- data.frame(high=thirty_kd[,2],low=thirty_kd[,3])
we_seventy_kd <- data.frame(high=seventy_kd[,3],low=seventy_kd[,4])
we_kd <- rbind(we_thirty_kd,we_seventy_kd)
we_kd <- as.data.frame(we_kd)

wilcox.exact(x=we_wt[,1],y=we_wt[,2],paired=F,alternative="g")[3]
wilcox.exact(x=we_kd[,1],y=we_kd[,2],paired=F,alternative="g")[3]

thirty_data <- thirty_wt  #change
#thirty_data <- thirty_kd  #change

seventy_data <- seventy_wt  #change
#seventy_data <- seventy_kd  #change

thirty_target <- thirty_data$name
thirty_data <- thirty_data[order(thirty_data$name),]
thirty_data <- thirty_data[,1:4]
p <- colnames(thirty_data)
t <- colnames(thirty_data)
thirty_data$t <- seq(1,dim(thirty_data)[1],1)
temp <- reshape2::melt(thirty_data,
                       id="t",
                       measure=p 
)
temp$variable <- as.character(temp$variable)
temp$value <- as.numeric(temp$value)
temp$t <- as.character(temp$t)
temp$t <- rep(paste(as.character(seq(1,length(thirty_target),1)),"_first",sep=""),4)
thirty_temp <- temp

seventy_target <- seventy_data$name
seventy_data <- seventy_data[order(seventy_data$name),]
seventy_data <- seventy_data[,1:4]
p <- colnames(seventy_data)
t <- colnames(seventy_data)
seventy_data$t <- seq(1,dim(seventy_data)[1],1)
temp <- reshape2::melt(seventy_data,
                       id="t",
                       measure=p 
)
temp$variable <- as.character(temp$variable)
temp$value <- as.numeric(temp$value)
temp$t <- as.character(temp$t)
temp$t <- rep(paste(as.character(seq(1,length(seventy_target),1)),"_second",sep=""),4)
seventy_temp <- temp

temp <- rbind(seventy_temp,thirty_temp)
temp <- as.data.frame(temp)

numbers <- seq(1,length(thirty_target),1)
numbers <- sort(numbers,decreasing = T)
thirty_numbers <- paste(numbers,"_first",sep="")

numbers <- seq(1,length(seventy_target),1)
numbers <- sort(numbers,decreasing = T)
seventy_numbers <- paste(numbers,"_second",sep="")

g <- ggplot(temp,aes(as.factor(variable),as.factor(t)))+
  geom_tile(aes(fill=value))+
  scale_fill_gradient2(name=guide_legend(title="z-score"),low="dodgerblue",mid="mediumpurple1",high="red",limits=c(-1.5,1.5))+
  labs(x = "Time (min)",y = "")+ 
  scale_x_discrete(breaks=c("0","30","75","120"),labels=waiver(),limits=c("0","30","75","120"),expand = c(0, 0))+ 
  scale_y_discrete(breaks=c(seventy_numbers,thirty_numbers),labels=c(sort(seventy_target,decreasing = T),sort(thirty_target,decreasing = T)),limits=c(seventy_numbers,thirty_numbers),expand = c(0, 0))+
  theme(legend.position = "right",axis.ticks = element_blank(), axis.title.x=element_text(size=9.0), axis.text.x = element_text(size = 7.0, angle = 0, hjust = 0.5, colour = "black"),axis.text.y = element_text(size = 6.0, angle = 0, hjust = 1, colour = "black",face="italic"));g


