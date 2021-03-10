library(VennDiagram)
library(genefilter)
library(reshape2)
library(caret)
library(ggplot2)

#Number of overlapping TNF-induced regions
siCtrl <- read.table("../Data/TMM_RPM_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks.txt",sep="\t",stringsAsFactors = FALSE)
siCtrl <- siCtrl[order(siCtrl[,2]),]
rownames(siCtrl) <- seq(1,dim(siCtrl)[1],1)
data <- siCtrl[,4:7]
rowname <- c()
for(i in 1:dim(data)[1]){
  if(as.character(max(data[i,1]) < max(data[i,2:4]))=="TRUE"){
    rowname[i] <- rownames(data)[i] 
  }else{
    next
  }
}
rowname <- na.omit(rowname)
all <- read.table("../Data/DifferentialPeaks_siCtrl_before_and_after_mergedpeaks.txt",sep="\t",stringsAsFactors = FALSE)
all <- all[order(all[,3]),]
rownames(all) <- seq(1,dim(all)[1],1)
all <- all[which(rownames(all) %in% rowname),]
thirty <- c()
seventyfive <- c()
hundred <- c()
for(i in 1:dim(all)[1]){
  if(as.character(nchar(all[i,9])==0)=="TRUE"){
    thirty[i] <- c(0)
  }else{
    thirty[i] <- c(1)
  }
  if(as.character(nchar(all[i,10])==0)=="TRUE"){
    seventyfive[i] <- c(0)
  }else{
    seventyfive[i] <- c(1)
  }
  if(as.character(nchar(all[i,11])==0)=="TRUE"){
    hundred[i] <- c(0)
  }else{
    hundred[i] <- c(1)
  }
}
thirty <- as.data.frame(thirty)
seventyfive <- as.data.frame(seventyfive)
hundred <- as.data.frame(hundred)
all <- cbind(thirty,seventyfive,hundred)
all <- as.data.frame(all)
all$total <- paste(all[,1],all[,2],all[,3],sep="")
one <- dim(all[which(all$total=="100"),])
two <- dim(all[which(all$total=="010"),])
three <- dim(all[which(all$total=="001"),])
four <- dim(all[which(all$total=="101"),])
five <- dim(all[which(all$total=="011"),])
six <- dim(all[which(all$total=="110"),])
seven <- dim(all[which(all$total=="111"),])

thirty<- one+four+six+seven
seventyfive <- two+five+six+seven
hundred <- three+four+five+seven

grid.newpage()
overrideTriple=T
draw.triple.venn(833,704,1889,98,115,48,17,category = c("", "",""), fill = c("darkorchid1", "darkorchid1","darkorchid1"), cex=0, alpha = c(0.5,0.5,0.5),cat.pos = c(50,50,50),scaled = TRUE,rotation.degree = 0)

#Cluster TNF-induced regions based on their time course chromatin accessibility
siCtrl <- read.table("../Data/TMM_RPM_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks.txt",sep="\t",stringsAsFactors = FALSE)
siCtrl <- siCtrl[order(siCtrl[,2]),]
rownames(siCtrl) <- seq(1,dim(siCtrl)[1],1)
data <- siCtrl[,4:7]
data <- na.omit(data)
data <- data[which(rownames(data) %in% rowname),]
data <- genescale(data,axis=1,method="Z")
cluster1 <- c()
cluster2 <- c()
cluster3 <- c()
cluster4 <- c()
cluster5 <- c()
cluster6 <- c()
cluster7 <- c()
cluster8 <- c()
cluster9 <- c()
cluster10 <- c()
cluster11 <- c()
cluster12 <- c()
cluster13 <- c()
cluster14 <- c()
cluster15 <- c()
cluster16 <- c()
for(i in 1:dim(data)[1]){
  if(as.character(data[i,1]>=0)=="TRUE"){
    if(as.character(data[i,2]>=0)=="TRUE"){
      if(as.character(data[i,3]>=0)=="TRUE"){
        if(as.character(data[i,4]>=0)=="TRUE"){
          cluster2[i] <- rownames(data)[i]
        }else{
          cluster3[i] <- rownames(data)[i]
        }
      }else{
        if(as.character(data[i,4]>=0)=="TRUE"){
          cluster4[i] <- rownames(data)[i]    
        }else{
          cluster7[i] <- rownames(data)[i]
        }
      }
    }else{
      if(as.character(data[i,3]>=0)=="TRUE"){  
        if(as.character(data[i,4]>=0)=="TRUE"){
          cluster5[i] <- rownames(data)[i]
        }else{
          cluster10[i] <- rownames(data)[i] 
        }
      }else{
        if(as.character(data[i,4]>=0)=="TRUE"){
          cluster8[i] <- rownames(data)[i]
        }else{
          cluster13[i] <- rownames(data)[i]
        }
      }
    }  
  }else{
    if(as.character(data[i,2]>=0)=="TRUE"){   
      if(as.character(data[i,3]>=0)=="TRUE"){
        if(as.character(data[i,4]>=0)=="TRUE"){
          cluster6[i] <- rownames(data)[i]
        }else{
          cluster12[i] <- rownames(data)[i]   
        }   
      }else{
        if(as.character(data[i,4]>=0)=="TRUE"){
          cluster11[i] <- rownames(data)[i]
        }else{
          cluster16[i] <- rownames(data)[i]
        }
      }
    }else{
      if(as.character(data[i,3]>=0)=="TRUE"){
        if(as.character(data[i,4]>=0)=="TRUE"){
          cluster9[i] <- rownames(data)[i]
        }else{
          cluster15[i] <- rownames(data)[i] 
        }
      }else{
        if(as.character(data[i,4]>=0)=="TRUE"){
          cluster14[i] <- rownames(data)[i]
        }else{
          cluster1[i] <- rownames(data)[i]
        }
      }
    }
  }
}
cluster1 <- na.omit(cluster1)
cluster2 <- na.omit(cluster2)
cluster3 <- na.omit(cluster3)
cluster4 <- na.omit(cluster4)
cluster5 <- na.omit(cluster5)
cluster6 <- na.omit(cluster6)
cluster7 <- na.omit(cluster7)
cluster8 <- na.omit(cluster8)
cluster9 <- na.omit(cluster9)
cluster10 <- na.omit(cluster10)
cluster11 <- na.omit(cluster11)
cluster12 <- na.omit(cluster12)
cluster13 <- na.omit(cluster13)
cluster14 <- na.omit(cluster14)
cluster15 <- na.omit(cluster15)
cluster16 <- na.omit(cluster16)
length(cluster1)
length(cluster2)
length(cluster3)
length(cluster4)
length(cluster5)
length(cluster6)
length(cluster7)
length(cluster8)
length(cluster9)
length(cluster10)
length(cluster11)
length(cluster12)
length(cluster13)
length(cluster14)
length(cluster15)
length(cluster16)

datacluster1 <- data[which(rownames(data) %in% cluster1),]
datacluster2 <- data[which(rownames(data) %in% cluster2),]
datacluster3 <- data[which(rownames(data) %in% cluster3),]
datacluster4 <- data[which(rownames(data) %in% cluster4),]
datacluster5 <- data[which(rownames(data) %in% cluster5),]
datacluster6 <- data[which(rownames(data) %in% cluster6),]
datacluster7 <- data[which(rownames(data) %in% cluster7),]
datacluster8 <- data[which(rownames(data) %in% cluster8),]
datacluster9 <- data[which(rownames(data) %in% cluster9),]
datacluster10 <- data[which(rownames(data) %in% cluster10),]
datacluster11 <- data[which(rownames(data) %in% cluster11),]
datacluster12 <- data[which(rownames(data) %in% cluster12),]
datacluster13 <- data[which(rownames(data) %in% cluster13),]
datacluster14 <- data[which(rownames(data) %in% cluster14),]
datacluster15 <- data[which(rownames(data) %in% cluster15),]
datacluster16 <- data[which(rownames(data) %in% cluster16),]

data <- rbind(datacluster1,datacluster2,datacluster3,datacluster4,datacluster5,datacluster6,datacluster7,datacluster8,datacluster9,datacluster10,datacluster11,datacluster12,datacluster13,datacluster14,datacluster15,datacluster16)
colnames(data) <- c("0min","30min","75min","120min")
data <- as.data.frame(data)

p <- colnames(data)
t <- colnames(data)
data2 <- data
data2$t <- seq(1,dim(data2)[1],1)
temp <- melt(data2,
             id="t",
             measure=p 
)
temp$variable <- as.character(temp$variable)
temp$value <- as.numeric(temp$value)

g <- ggplot(temp,aes(as.factor(variable),as.factor(t)))+
  geom_tile(aes(fill=value),stat="identity")+
  scale_fill_gradient(low="dodgerblue",high="red")+
  labs(x = "",y = "")+ 
  theme(panel.grid = element_blank())+
  scale_x_discrete(labels=c("0min","30min","75min","120min"),breaks=c("0min","30min","75min","120min"),limits=c("0min","30min","75min","120min"),expand = c(0, 0))+ 
  theme(legend.position = "right",axis.ticks = element_blank(), axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, colour = "black"),axis.text.y = element_blank());g

#Draw a heatmap of aggregated chromatin accessibility at each kB site enriched TNF-induced region
one <- read.table("../Data/204regions_regions.txt", sep = "\t")
two <- read.table("../Data/453regions_regions.txt", sep = "\t")
three <- read.table("../Data/348regions_regions.txt", sep = "\t")
four <- read.table("../Data/478regions_regions.txt", sep = "\t")
five <- read.table("../Data/119regions_regions.txt", sep = "\t")

siCtrl <- read.table("../Data/TMM_RPM_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks.txt",sep="\t",stringsAsFactors = FALSE)
siCtrl <- siCtrl[order(siCtrl[,2]),]
rownames(siCtrl) <- seq(1,dim(siCtrl)[1],1)
siCtrl[,2] <- siCtrl[,2]+1
siCtrl_one <- siCtrl[which(siCtrl[,2] %in% one[,3]),4:7]
siCtrl_two <- siCtrl[which(siCtrl[,2] %in% two[,3]),4:7]
siCtrl_three <- siCtrl[which(siCtrl[,2] %in% three[,3]),4:7]
siCtrl_four <- siCtrl[which(siCtrl[,2] %in% four[,3]),4:7]
siCtrl_five <- siCtrl[which(siCtrl[,2] %in% five[,3]),4:7]

dfsiCtrl <- siCtrl_one
dfsiCtrl <- siCtrl_two
dfsiCtrl <- siCtrl_three
dfsiCtrl <- siCtrl_four
dfsiCtrl <- siCtrl_five

dfsiCtrl <- as.data.frame(dfsiCtrl)
rownames(dfsiCtrl) <- seq(1,dim(dfsiCtrl)[1],1)
colnames(dfsiCtrl) <-  c("siCtrl_0","siCtrl_30","siCtrl_75","siCtrl_120")
dfsiCtrl <- genescale(dfsiCtrl,axis=1,method="Z")
dfsiCtrl <- as.data.frame(dfsiCtrl)
dfsiCtrl <- dfsiCtrl[,1:4]
dfsiCtrl <- as.data.frame(dfsiCtrl)

p <- colnames(dfsiCtrl)
t <- colnames(dfsiCtrl)
dfsiCtrl$t <- seq(1,dim(dfsiCtrl)[1],1)
temp <- reshape2::melt(dfsiCtrl,
                       id="t",
                       measure=p 
)
temp$variable <- as.character(temp$variable)
temp$value <- as.numeric(temp$value)
temp$t <- as.character(temp$t)
g <- ggplot(temp,aes(as.factor(variable),as.factor(t)))+
  geom_tile(aes(fill=value),)+
  scale_fill_gradient(low="dodgerblue",high="red")+
  labs(x = "",y = "")+ 
  scale_x_discrete(breaks=c("siCtrl_0","siCtrl_30","siCtrl_75","siCtrl_120"),labels=waiver(),limits=c("siCtrl_0","siCtrl_30","siCtrl_75","siCtrl_120"),expand = c(0, 0))+ 
  #scale_x_discrete(breaks=c("siIkBa_0","siIkBa_30","siIkBa_75","siIkBa_120"),labels=waiver(),limits=c("siIkBa_0","siIkBa_30","siIkBa_75","siIkBa_120"),expand = c(0, 0))+ 
  scale_y_discrete(breaks=order(seq(1,temp$t[dim(temp)[1]],by=1),decreasing = FALSE),labels=order(seq(1,temp$t[dim(temp)[1]],by=1),decreasing = FALSE),limits=order(seq(1,temp$t[dim(temp)[1]],by=1),decreasing = FALSE),expand = c(0, 0))+
  theme(legend.position = "right",axis.ticks = element_blank(), axis.text.x = element_text(size = 5, angle = 0, hjust = 0, colour = "black"),axis.text.y = element_text(size = 3, angle = 0, hjust = 0, colour = "black"));g

#Draw a heatmap of aggregated chromatin accessibility at all kB site enriched TNF-induced regions
one <- read.table("../Data/204regions_regions.txt", sep = "\t")
two <- read.table("../Data/453regions_regions.txt", sep = "\t")
three <- read.table("../Data/348regions_regions.txt", sep = "\t")
four <- read.table("../Data/478regions_regions.txt", sep = "\t")
five <- read.table("../Data/119regions_regions.txt", sep = "\t")

one <- read.table("../Data/177regions_regions.txt", sep = "\t")
two <- read.table("../Data/446regions_regions.txt", sep = "\t")
three <- read.table("../Data/368regions_regions.txt", sep = "\t")
four <- read.table("../Data/397regions_regions.txt", sep = "\t")
five <- read.table("../Data/361regions_regions.txt", sep = "\t")

#siCtrl <- read.table("../Data/TMM_RPM_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks.txt",sep="\t",stringsAsFactors = FALSE)
siCtrl <- read.table("../Data/TMM_RPM_in_siCtrl_at_DifferentialPeaks_mergedpeaks_in_siIkBa.txt",sep="\t",stringsAsFactors = FALSE)
siCtrl <- siCtrl[order(siCtrl[,2]),]
rownames(siCtrl) <- seq(1,dim(siCtrl)[1],1)
siCtrl[,2] <- siCtrl[,2]+1
siCtrl_one <- siCtrl[which(siCtrl[,2] %in% one[,3]),4:7]
siCtrl_two <- siCtrl[which(siCtrl[,2] %in% two[,3]),4:7]
siCtrl_three <- siCtrl[which(siCtrl[,2] %in% three[,3]),4:7]
siCtrl_four <- siCtrl[which(siCtrl[,2] %in% four[,3]),4:7]
siCtrl_five <- siCtrl[which(siCtrl[,2] %in% five[,3]),4:7]

#siIkBa <- read.table("../Data/TMM_RPM_in_siIkBa_at_DifferentialPeaks_mergedpeaks_in_siCtrl.txt",sep="\t",stringsAsFactors = FALSE)
siIkBa <- read.table("../Data/TMM_RPM_DifferentialPeaks_siIkBa_before_and_after_mergedpeaks.txt",sep="\t",stringsAsFactors = FALSE)
siIkBa <- siIkBa[order(siIkBa[,2]),]
rownames(siIkBa) <- seq(1,dim(siIkBa)[1],1)
siIkBa[,2] <- siIkBa[,2]+1
siIkBa_one <- siIkBa[which(siIkBa[,2] %in% one[,3]),4:7]
siIkBa_two <- siIkBa[which(siIkBa[,2] %in% two[,3]),4:7]
siIkBa_three <- siIkBa[which(siIkBa[,2] %in% three[,3]),4:7]
siIkBa_four <- siIkBa[which(siIkBa[,2] %in% four[,3]),4:7]
siIkBa_five <- siIkBa[which(siIkBa[,2] %in% five[,3]),4:7]

dfsiCtrl <- rbind(siCtrl_one,siCtrl_two,siCtrl_three,siCtrl_four,siCtrl_five)
dfsiCtrl <- as.data.frame(dfsiCtrl)
rownames(dfsiCtrl) <- seq(1,dim(dfsiCtrl)[1],1)
colnames(dfsiCtrl) <-  c("siCtrl_0","siCtrl_30","siCtrl_75","siCtrl_120")
dfsiIkBa <- rbind(siIkBa_one,siIkBa_two,siIkBa_three,siIkBa_four,siIkBa_five)
dfsiIkBa <- as.data.frame(dfsiIkBa)
rownames(dfsiIkBa) <- seq(1,dim(dfsiIkBa)[1],1)
colnames(dfsiIkBa) <-  c("siIkBa_0","siIkBa_30","siIkBa_75","siIkBa_120")
df <- cbind(dfsiCtrl,dfsiIkBa)
df <- as.data.frame(df)
df <- genescale(df,axis=1,method="Z")
df <- as.data.frame(df)

df <- df[,1:4]
df <- df[,5:8]

df <- as.data.frame(df)

p <- colnames(df)
t <- colnames(df)
df$t <- seq(1,dim(df)[1],1)
temp <- reshape2::melt(df,
             id="t",
             measure=p 
)
temp$variable <- as.character(temp$variable)
temp$value <- as.numeric(temp$value)
temp$t <- as.character(temp$t)
g <- ggplot(temp,aes(as.factor(variable),as.factor(t)))+
  geom_tile(aes(fill=value),)+
  scale_fill_gradient(low="dodgerblue",high="red")+
  labs(x = "",y = "")+ 
  #scale_x_discrete(breaks=c("siCtrl_0","siCtrl_30","siCtrl_75","siCtrl_120"),labels=waiver(),limits=c("siCtrl_0","siCtrl_30","siCtrl_75","siCtrl_120"),expand = c(0, 0))+ 
  scale_x_discrete(breaks=c("siIkBa_0","siIkBa_30","siIkBa_75","siIkBa_120"),labels=waiver(),limits=c("siIkBa_0","siIkBa_30","siIkBa_75","siIkBa_120"),expand = c(0, 0))+ 
  scale_y_discrete(breaks=order(seq(1,temp$t[dim(temp)[1]],by=1),decreasing = TRUE),labels=order(seq(1,temp$t[dim(temp)[1]],by=1),decreasing = TRUE),limits=order(seq(1,temp$t[dim(temp)[1]],by=1),decreasing = TRUE),expand = c(0, 0))+
  theme(legend.position = "right",axis.ticks = element_blank(), axis.text.x = element_text(size = 5, angle = 0, hjust = 0, colour = "black"),axis.text.y = element_text(size = 3, angle = 0, hjust = 0, colour = "black"));g

#Draw boxplots of aggregated chromatin accessibility between WT and siIkBa
diff2 <- read.table("../Data/478regions_regions.txt", sep = "\t")

siCtrl <- read.table("../Data/TMM_RPM_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks.txt",sep="\t",stringsAsFactors = FALSE)
#siCtrl <- read.table("../Data/TMM_RPM_in_siCtrl_at_DifferentialPeaks_mergedpeaks_in_siIkBa.txt",sep="\t",stringsAsFactors = FALSE)
siCtrl <- siCtrl[order(siCtrl[,2]),]
rownames(siCtrl) <- seq(1,dim(siCtrl)[1],1)
siCtrl[,2] <- siCtrl[,2]+1
siCtrl <- siCtrl[which(siCtrl[,2] %in% diff2[,3]),]

siIkBa <- read.table("../Data/TMM_RPM_in_siIkBa_at_DifferentialPeaks_mergedpeaks_in_siCtrl.txt",sep="\t",stringsAsFactors = FALSE)
#siIkBa <- read.table("../Data/TMM_RPM_DifferentialPeaks_siIkBa_before_and_after_mergedpeaks.txt",sep="\t",stringsAsFactors = FALSE)
siIkBa <- siIkBa[order(siIkBa[,2]),]
rownames(siIkBa) <- seq(1,dim(siIkBa)[1],1)
siIkBa[,2] <- siIkBa[,2]+1
siIkBa <- siIkBa[which(siIkBa[,2] %in% diff2[,3]),]

siCtrl <- siCtrl[,4:7]
siCtrl <- as.data.frame(siCtrl)
siIkBa <- siIkBa[,4:7]
siIkBa <- as.data.frame(siIkBa)
df0 <- cbind(siCtrl,siIkBa)
df0 <- as.data.frame(df0)
df0 <- genescale(df0,axis=1,method="Z")
df0 <- as.data.frame(df0)
siCtrl <- df0[,1:4]
siIkBa <- df0[,5:8]

siCtrl <- t(siCtrl)
siCtrl <- as.data.frame(siCtrl)
siCtrl$t <- c("0","30","75","120")
p <- colnames(siCtrl)
temp <- melt(siCtrl,
             id="t",
             measure=p 
)
temp <- temp[1:(dim(temp)[1]-4),]
temp$group <- rep("siCtrl",dim(temp)[1])
siCtrl <- temp

siIkBa <- t(siIkBa)
siIkBa <- as.data.frame(siIkBa)
siIkBa$t <- c("0","30","75","120")
p <- colnames(siIkBa)
temp <- melt(siIkBa,
             id="t",
             measure=p 
)
temp <- temp[1:(dim(temp)[1]-4),]
temp$group <- rep("siIkBa",dim(temp)[1])
siIkBa <- temp

siCtrl$value <- as.numeric(siCtrl$value)
siIkBa$value <- as.numeric(siIkBa$value)
df <- rbind(siCtrl,siIkBa)
df$value <- as.numeric(df$value)

g <- ggplot(df,
            aes(x=t,
                y=value,
                fill=group))+ 
  theme_bw(base_size = 16)+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values = c(siCtrl="magenta",siIkBa="mediumseagreen"))+
  scale_x_discrete(breaks=c("0","30","75","120"),labels=c("0","30","75","120"),limits=c("0","30","75","120"),expand=c(0.2,0.2))+
  scale_y_continuous(breaks=seq(-3,6,by=1),limits=c(-3,6),labels=waiver())+
  labs(x="Time (min)",y="ATAC-seq signal",title="",colour="")+
  geom_boxplot(position=position_dodge(0.9))+
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),legend.position = "true",panel.border=element_rect(fill=NA,color="black", size=0.5),axis.title.x = element_text(size=1,margin=margin(t = 35, r = 0, b = 0, l = 0)),axis.title.y = element_text(size=1,margin=margin(t = 0, r = 25, b = 0, l = 0)),axis.text.x = element_text(size=1),axis.text.y = element_text(size=1));g

#Paired t-test of aggregated chromatin accessibility between WT and siIkBa
diff2 <- read.table("../Data/453regions_regions.txt", sep = "\t")

siCtrl <- read.table("../Data/TMM_RPM_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks.txt",sep="\t",stringsAsFactors = FALSE)
#siCtrl <- read.table("../Data/TMM_RPM_in_siCtrl_at_DifferentialPeaks_mergedpeaks_in_siIkBa.txt",sep="\t",stringsAsFactors = FALSE)
siCtrl <- siCtrl[order(siCtrl[,2]),]
rownames(siCtrl) <- seq(1,dim(siCtrl)[1],1)
siCtrl[,2] <- siCtrl[,2]+1
siCtrl <- siCtrl[which(siCtrl[,2] %in% diff2[,3]),]

siIkBa <- read.table("../Data/TMM_RPM_in_siIkBa_at_DifferentialPeaks_mergedpeaks_in_siCtrl.txt",sep="\t",stringsAsFactors = FALSE)
#siIkBa <- read.table("../Data/TMM_RPM_DifferentialPeaks_siIkBa_before_and_after_mergedpeaks.txt",sep="\t",stringsAsFactors = FALSE)
siIkBa <- siIkBa[order(siIkBa[,2]),]
rownames(siIkBa) <- seq(1,dim(siIkBa)[1],1)
siIkBa[,2] <- siIkBa[,2]+1
siIkBa <- siIkBa[which(siIkBa[,2] %in% diff2[,3]),]

siCtrl <- siCtrl[,4:7]
siCtrl <- as.data.frame(siCtrl)
siIkBa <- siIkBa[,4:7]
siIkBa <- as.data.frame(siIkBa)

#Paired t-test of chromatin accessibility in aggregated cells
t.test(x=siCtrl[,1],y=siIkBa[,1],paired=T,alternative="greater")
t.test(x=siCtrl[,2],y=siIkBa[,2],paired=T,alternative="greater")
t.test(x=siCtrl[,3],y=siIkBa[,3],paired=T,alternative="greater")
t.test(x=siCtrl[,4],y=siIkBa[,4],paired=T,alternative="greater")

t.test(x=siCtrl[,1],y=siIkBa[,1],paired=T,alternative="less")
t.test(x=siCtrl[,2],y=siIkBa[,2],paired=T,alternative="less")
t.test(x=siCtrl[,3],y=siIkBa[,3],paired=T,alternative="less")
t.test(x=siCtrl[,4],y=siIkBa[,4],paired=T,alternative="less")

#Accuracy of classifying single-cell chromatin accessibility between WT and siIkBa using 10-fold cross validation of Bayesian generalized linear model
data <- read.table("../Data/DifferentialPeaks_siCtrl_before_and_after_mergedpeaks.txt", sep = "\t",header=FALSE)
data <- data[order(data[,3]),]
rownames(data) <- seq(1,dim(data)[1],1)
diff2 <- read.table("../Data/453regions_regions.txt", sep = "\t")
data <- data[which(data[,3] %in% diff2[,3]),]
hit <- data

zero_1 <- read.table("../Data/TMM_RPM_single-cell_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks_siCtrl_0_rep1.txt",sep="\t",stringsAsFactors = FALSE)
zero_2 <- read.table("../Data/TMM_RPM_single-cell_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks_siCtrl_0_rep2.txt",sep="\t",stringsAsFactors = FALSE)
zero_3 <- read.table("../Data/TMM_RPM_single-cell_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks_siCtrl_0_rep3.txt",sep="\t",stringsAsFactors = FALSE)
three_1 <- read.table("../Data/TMM_RPM_single-cell_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks_siCtrl_30_rep1.txt",sep="\t",stringsAsFactors = FALSE)
three_2 <- read.table("../Data/TMM_RPM_single-cell_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks_siCtrl_30_rep2.txt",sep="\t",stringsAsFactors = FALSE)
three_3 <- read.table("../Data/TMM_RPM_single-cell_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks_siCtrl_30_rep3.txt",sep="\t",stringsAsFactors = FALSE)
seven_1 <- read.table("../Data/TMM_RPM_single-cell_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks_siCtrl_75_rep1.txt",sep="\t",stringsAsFactors = FALSE)
seven_2 <- read.table("../Data/TMM_RPM_single-cell_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks_siCtrl_75_rep2.txt",sep="\t",stringsAsFactors = FALSE)
seven_3 <- read.table("../Data/TMM_RPM_single-cell_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks_siCtrl_75_rep3.txt",sep="\t",stringsAsFactors = FALSE)
hundred_1 <- read.table("../Data/TMM_RPM_single-cell_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks_siCtrl_120_rep1.txt",sep="\t",stringsAsFactors = FALSE)
hundred_2 <- read.table("../Data/TMM_RPM_single-cell_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks_siCtrl_120_rep2.txt",sep="\t",stringsAsFactors = FALSE)
hundred_3 <- read.table("../Data/TMM_RPM_single-cell_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks_siCtrl_120_rep3.txt",sep="\t",stringsAsFactors = FALSE)

data <- cbind(zero_1[,4:dim(zero_1)[2]],zero_2[,4:dim(zero_2)[2]],zero_3[,4:dim(zero_3)[2]])
data <- cbind(three_1[,4:dim(three_1)[2]],three_2[,4:dim(three_2)[2]],three_3[,4:dim(three_3)[2]])
data <- cbind(seven_1[,4:dim(seven_1)[2]],seven_2[,4:dim(seven_2)[2]],seven_3[,4:dim(seven_3)[2]])
data <- cbind(hundred_1[,4:dim(hundred_1)[2]],hundred_2[,4:dim(hundred_2)[2]],hundred_3[,4:dim(hundred_3)[2]])

data <- data[which(rownames(data) %in% rownames(hit)),]
siCtrl <- data

zero_1 <- read.table("../Data/TMM_RPM_single-cell_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks_siIkBa_0_rep1.txt",sep="\t",stringsAsFactors = FALSE)
zero_2 <- read.table("../Data/TMM_RPM_single-cell_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks_siIkBa_0_rep2.txt",sep="\t",stringsAsFactors = FALSE)
zero_3 <- read.table("../Data/TMM_RPM_single-cell_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks_siIkBa_0_rep3.txt",sep="\t",stringsAsFactors = FALSE)
three_1 <- read.table("../Data/TMM_RPM_single-cell_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks_siIkBa_30_rep1.txt",sep="\t",stringsAsFactors = FALSE)
three_2 <- read.table("../Data/TMM_RPM_single-cell_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks_siIkBa_30_rep2.txt",sep="\t",stringsAsFactors = FALSE)
three_3 <- read.table("../Data/TMM_RPM_single-cell_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks_siIkBa_30_rep3.txt",sep="\t",stringsAsFactors = FALSE)
seven_1 <- read.table("../Data/TMM_RPM_single-cell_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks_siIkBa_75_rep1.txt",sep="\t",stringsAsFactors = FALSE)
seven_2 <- read.table("../Data/TMM_RPM_single-cell_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks_siIkBa_75_rep2.txt",sep="\t",stringsAsFactors = FALSE)
seven_3 <- read.table("../Data/TMM_RPM_single-cell_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks_siIkBa_75_rep3.txt",sep="\t",stringsAsFactors = FALSE)
hundred_1 <- read.table("../Data/TMM_RPM_single-cell_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks_siIkBa_120_rep1.txt",sep="\t",stringsAsFactors = FALSE)
hundred_2 <- read.table("../Data/TMM_RPM_single-cell_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks_siIkBa_120_rep2.txt",sep="\t",stringsAsFactors = FALSE)
hundred_3 <- read.table("../Data/TMM_RPM_single-cell_DifferentialPeaks_siCtrl_before_and_after_mergedpeaks_siIkBa_120_rep3.txt",sep="\t",stringsAsFactors = FALSE)

data <- cbind(zero_1[,4:dim(zero_1)[2]],zero_2[,4:dim(zero_2)[2]],zero_3[,4:dim(zero_3)[2]])
data <- cbind(three_1[,4:dim(three_1)[2]],three_2[,4:dim(three_2)[2]],three_3[,4:dim(three_3)[2]])
data <- cbind(seven_1[,4:dim(seven_1)[2]],seven_2[,4:dim(seven_2)[2]],seven_3[,4:dim(seven_3)[2]])
data <- cbind(hundred_1[,4:dim(hundred_1)[2]],hundred_2[,4:dim(hundred_2)[2]],hundred_3[,4:dim(hundred_3)[2]])

data <- data[which(rownames(data) %in% rownames(hit)),]
siIkBa <- data

data <- cbind(siCtrl,siIkBa)
data <- as.data.frame(data)
df <- data
df <- t(df)
df <- as.data.frame(df)
df$group <- c(rep("siCtrl",dim(siCtrl)[2]),rep("siIkBa",dim(siIkBa)[2]))

set.seed(1)
model <- train(group ~ ., df,method = "bayesglm",
               trControl = trainControl(
                 method = "cv", 
                 number = 10,
                 verboseIter = TRUE
               )
)


