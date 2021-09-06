library(ggplot2)
library(scales)
library(pracma)
library(openxlsx)

count <- read.csv("160712_ikbakd.csv")
count1=count[,2:12]
count1 <- data.frame(count1)
mean <- count1[,2]
condition <- count1[,1]
time <- count1[,11]
df <- data.frame(Time=time,Mean=mean,Con=condition)
df1 <- df[1:13,]
x <- df1[,1]                                         
y <- df1[,2] 
signal <- data.frame(time = x,nfkb = y)
wtsignal <- signal
wtsignal$group <- rep("siCtrl",dim(wtsignal)[1])
wtsignal$time <- as.numeric(wtsignal$time)
wtrep1signal <- wtsignal

count <- read.csv("Ando_160617_ikbakd.csv")
count1=count[,2:12]
count1 <- data.frame(count1)
mean <- count1[,2]
condition <- count1[,1]
time <- count1[,11]
df <- data.frame(Time=time,Mean=mean,Con=condition)
df1 <- df[1:13,]
x <- df1[,1]                                         
y <- df1[,2] 
signal <- data.frame(time = x,nfkb = y)
wtsignal <- signal
wtsignal$group <- rep("siCtrl",dim(wtsignal)[1])
wtsignal$time <- as.numeric(wtsignal$time)
wtrep2signal <- wtsignal

wtsignal <- wtrep1signal
two <- cbind(wtrep1signal$nfkb,wtrep2signal$nfkb)
two <- as.data.frame(two)
two$mean <- apply(two,1,mean)
wtsignal$nfkb <- two$mean

fcwtsignal <- c()
first <- wtsignal$nfkb[1]
for(i in 1:13){
  fcwtsignal[i] <- wtsignal$nfkb[i]/first
}
wtsignal$nfkb <- fcwtsignal

count <- read.csv("160712_ikbakd.csv")
count1=count[,2:12]
count1 <- data.frame(count1)
mean <- count1[,2]
condition <- count1[,1]
time <- count1[,11]
df <- data.frame(Time=time,Mean=mean,Con=condition)
df2 <- df[26:38,]
x <- df2[,1]
y <- df2[,2] 
signal <- data.frame(time = x,nfkb = y)
kdsignal <- signal
kdsignal$group <- rep("siIkBa",dim(kdsignal)[1])
kdsignal$time <- as.numeric(kdsignal$time)
kdrep1signal <- kdsignal

count <- read.csv("Ando_160617_ikbakd.csv")
count1=count[,2:12]
count1 <- data.frame(count1)
mean <- count1[,2]
condition <- count1[,1]
time <- count1[,11]
df <- data.frame(Time=time,Mean=mean,Con=condition)
df2 <- df[26:38,]
x <- df2[,1]
y <- df2[,2] 
signal <- data.frame(time = x,nfkb = y)
kdsignal <- signal
kdsignal$group <- rep("siIkBa",dim(kdsignal)[1])
kdsignal$time <- as.numeric(kdsignal$time)
kdrep2signal <- kdsignal

kdsignal <- kdrep1signal
two <- cbind(kdrep1signal$nfkb,kdrep2signal$nfkb)
two <- as.data.frame(two)
two$mean <- apply(two,1,mean)
kdsignal$nfkb <- two$mean

fckdsignal <- c()
first <- kdsignal$nfkb[1]
for(i in 1:13){
  fckdsignal[i] <- kdsignal$nfkb[i]/first
}
kdsignal$nfkb <- fckdsignal

df <- data.frame(time=seq(0,10800,900),mrna=wtsignal$nfkb)
df <- df[1:13,]                                         
fp <- pchipfun(df[,1], df[,2]) 
f <- fp(seq(0, 10800, by = 1))
df <- data.frame(time=seq(0,10800,1),mrna=f)
signal <- data.frame(time = df[,1],nfkb = df[,2])
nfkb_wt <- signal

df <- data.frame(time=seq(0,10800,900),mrna=kdsignal$nfkb)
df <- df[1:13,]                                         
fp <- pchipfun(df[,1], df[,2]) 
f <- fp(seq(0, 10800, by = 1))
df <- data.frame(time=seq(0,10800,1),mrna=f)
signal <- data.frame(time = df[,1],nfkb = df[,2])
nfkb_kd <- signal

df <- rbind(nfkb_wt,nfkb_kd)
df <- as.data.frame(df)
df <- rescale(df[,2], to=c(0,1))
df <- data.frame(time=c(nfkb_wt[,1],nfkb_kd[,1]),nfkb=df)
nfkb_wt <- df[1:10801,]
nfkb_kd <- df[10802:21602,]

dots <- 60*seq(0,180,15)
nfkb_wt$group <- rep("none",dim(nfkb_wt)[1])
nfkb_wt[which(nfkb_wt$time %in% dots),3]="dots"
nfkb_kd$group <- rep("none",dim(nfkb_kd)[1])
nfkb_kd[which(nfkb_kd$time %in% dots),3]="dots"

g <- ggplot(nfkb_wt,aes(x=time,y=nfkb,size=group))+
  theme_bw(base_size = 16)+
  theme(panel.grid = element_blank(),legend.position = "true")+
  labs(x="",y="",title="",colour="")+
  scale_x_continuous(breaks=seq(0,10800,900),labels=waiver(),limits=c(0,10800))+
  scale_y_continuous(breaks=seq(0.0,1.0,0.25),labels=waiver(),limits=c(0.0,1.0),expand=c(0,0.1))+
  geom_point(color="magenta",shape=19,aes(size=group))+
  #geom_point(color="seagreen",shape=19,aes(size=group))+
  scale_size_manual(values=c(none=0.0,dots=4.0))+
  geom_line(lwd=1.2,color="magenta")+
  #geom_line(lwd=1.2,color="seagreen")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=20,face="italic"),panel.background=element_rect(color="black", size=1.9, fill="white" ),panel.grid.minor.y = element_blank(),panel.grid.major.x = element_line(colour="grey86", size =1.0),panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),legend.position = "none",axis.title.x = element_text(size=15,margin=margin(t = 25, r = 0, b = 0, l = 0)),axis.title.y = element_text(size=16,margin=margin(t = 0, r = 25, b = 0, l = 0)),axis.text.x = element_text(size=15),axis.text.y = element_text(size=15));g