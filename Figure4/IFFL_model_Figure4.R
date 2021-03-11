library(edgeR)
library(pracma)
library(deSolve)
library(DEoptim)
library(scales)
library(bio3d)
library(car)
library(openxlsx)
library(ggplot2)

#Draw a heatmap of data
names <- read.table("../Data/ERGsubcluster2_Promoter_up500down500TSS_12012019.txt", sep = "\t",header=FALSE)
names <- names[-which(names[,1]=="NFKBIA"),]
names <- read.table("../Data/IRGsubcluster2_Promoter_up500down500TSS_12012019.txt", sep = "\t",header=FALSE)
names[,1] <- as.character(names[,1])    
names[which(names[,1]=="NKX3_1"),1] <- c("NKX3-1")
names <- read.table("../Data/DRGsubcluster2_Promoter_up500down500TSS_12012019.txt", sep = "\t",header=FALSE)
names <- names[,1]
names <- sort(names)

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
  
  df <- rescale(df[,2], to=c(0,1))
  refwtfc <- data.frame(x=seq(0,180,15),y=df[1:13])
  refkdfc <- data.frame(x=seq(0,180,15),y=df[14:26])
  
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
#refall <- refall[1:13,]   #siCtrl
refall <- refall[14:26,]   #siIkBa     

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
  #scale_fill_gradient(low="white",high="magenta",limits=c(0,1.0))+
  scale_fill_gradient(low="white",high="mediumseagreen",limits=c(0,1.0))+
  labs(x = "",y = "") + 
  theme(panel.grid = element_blank())+
  scale_x_discrete(labels=c("0","15","30","45","60","75","90","105","120","135","150","165","180"),breaks=c("0","15","30","45","60","75","90","105","120","135","150","165","180"),expand = c(0, 0))+ 
  scale_y_discrete(breaks=sort(names,decreasing = T),labels=waiver(),limits=sort(names,decreasing = T),expand = c(0, 0))+
  theme(legend.position = "bottom",axis.ticks = element_blank(), axis.text.x = element_text(size = 7.0, angle = 0, hjust = 0.5, colour = "black"),axis.text.y=element_text(size = 5.6, angle = 0, hjust = 1, colour = "black",face="italic"));g

#Draw a heatmap of simulation from IFFL model
kdegTF <- log(2)/(24*60*60)
KDTF <- 100

D2wt <- function(kdegTF,KDTF) {
  parameters2 <- c(kdegTF=kdegTF,KDTF=KDTF)   
  Lorenz2 <- function(times, state, parameters2){
    with(as.list(c(state, parameters2)), {
      
      dTF <- (((KDTF*inputwt(times))/(1+KDTF*inputwt(times))))-kdegTF*TF
      
      return(list(c(dTF)))
    })
  }
  
  TF0  <- (((KDTF*inputwt(0))/(1+KDTF*inputwt(0))))/kdegTF
  
  state <- c(TF=TF0)
  
  out0 <- ode(y = state, times = seq(-3600,0,by=1),
              func = Lorenz2, parms = parameters2)
  out0 <- as.data.frame(out0)
  return(out0)
}

D4wt <- function(kdegTF,KDTF) {
  parameters4 <- c(kdegTF=kdegTF,KDTF=KDTF)   
  
  TFfirst <- D2wt(kdegTF,KDTF)[3601,2]
  
  Lorenz4 <- function(times, state, parameters4){
    with(as.list(c(state, parameters4)), {
      
      dTF <- (((KDTF*inputwt(times))/(1+KDTF*inputwt(times))))-kdegTF*TF
      
      return(list(c(dTF)))
      
    })
  }
  
  state <- c(TF=TFfirst)
  
  out0 <- ode(y = state, times = seq(0,10800,by=1),
              func = Lorenz4, parms = parameters4)
  out0 <- as.data.frame(out0)
  return(out0)
}



D3wt <- function(kdeg,KD1,KD2) {
  parameters3 <- c(kdeg=kdeg,KD1=KD1,KD2=KD2)  
  
  signal2 <- data.frame(time = D2wt(kdegTF,KDTF)[,1],tf = D2wt(kdegTF,KDTF)[,2])
  sigimp2 <- approxfun(signal2$time, signal2$tf, rule = 2)
  input2wt <- sigimp2
  
  Lorenz3 <- function(times, state, parameters3){
    with(as.list(c(state, parameters3)), {
      
      dmrna = ((KD1*inputwt(times))^2)/(((KD1*inputwt(times))^2)+((KD2*input2wt(times))^2)+1)-kdeg*mrna
      
      return(list(c(dmrna)))
    })
  }
  
  mrna0 = ((KD1*inputwt(0))^2)/(((KD1*inputwt(0))^2)+((KD2*input2wt(0))^2)+1)/kdeg
  
  state <- c(mrna=mrna0)
  
  out0 <- ode(y = state, times = seq(-3600,0,by=1),
              func = Lorenz3, parms = parameters3)
  out0 <- as.data.frame(out0)
  return(out0)
}

D5wt <- function(kdeg,KD1,KD2,tau) {
  parameters5 <- c(kdeg=kdeg,KD1=KD1,KD2=KD2,tau=tau) 
  
  signal3 <- data.frame(time = D4wt(kdegTF,KDTF)[,1],tf = D4wt(kdegTF,KDTF)[,2])
  sigimp3 <- approxfun(signal3$time, signal3$tf, rule = 2)
  input3wt <- sigimp3
  
  mrnafirst <- D3wt(kdeg,KD1,KD2)[3601,2]
  
  Lorenz5 <- function(times, state, parameters5){
    with(as.list(c(state, parameters5)), {
      
      dmrna = ((KD1*inputwt(times-tau))^2)/(((KD1*inputwt(times-tau))^2)+((KD2*input3wt(times-7200))^2)+1)-kdeg*mrna
      
      return(list(c(dmrna)))
    })
  }
  
  state <- c(mrna=mrnafirst)
  
  out0 <- ode(y = state, times = seq(0,10800,by=1),
              func = Lorenz5, parms = parameters5)
  out0 <- as.data.frame(out0)
  return(out0)
}

D1wt <- function(kdeg,KD1,KD2,tau) {
  parameters <- c(kdeg=kdeg,KD1=KD1,KD2=KD2,tau=tau)
  
  D3wt <- function(kdeg,KD1,KD2) {
    parameters3 <- c(kdeg=kdeg,KD1=KD1,KD2=KD2)  
    
    signal2 <- data.frame(time = D2wt(kdegTF,KDTF)[,1],tf = D2wt(kdegTF,KDTF)[,2])
    sigimp2 <- approxfun(signal2$time, signal2$tf, rule = 2)
    input2wt <- sigimp2
    
    Lorenz3 <- function(times, state, parameters3){
      with(as.list(c(state, parameters3)), {
        
        dmrna = ((KD1*inputwt(times))^2)/(((KD1*inputwt(times))^2)+((KD2*input2wt(times))^2)+1)-kdeg*mrna
        
        return(list(c(dmrna)))
      })
    }
    
    mrna0 = ((KD1*inputwt(0))^2)/(((KD1*inputwt(0))^2)+((KD2*input2wt(0))^2)+1)/kdeg
    
    state <- c(mrna=mrna0)
    
    out0 <- ode(y = state, times = seq(-3600,0,by=1),
                func = Lorenz3, parms = parameters3)
    out0 <- as.data.frame(out0)
    return(out0)
  }
  
  D5wt <- function(kdeg,KD1,KD2,tau) {
    parameters5 <- c(kdeg=kdeg,KD1=KD1,KD2=KD2,tau=tau) 
    
    signal3 <- data.frame(time = D4wt(kdegTF,KDTF)[,1],tf = D4wt(kdegTF,KDTF)[,2])
    sigimp3 <- approxfun(signal3$time, signal3$tf, rule = 2)
    input3wt <- sigimp3
    
    mrnafirst <- D3wt(kdeg,KD1,KD2)[3601,2]
    
    Lorenz5 <- function(times, state, parameters5){
      with(as.list(c(state, parameters5)), {
        
        dmrna = ((KD1*inputwt(times-tau))^2)/(((KD1*inputwt(times-tau))^2)+((KD2*input3wt(times-7200))^2)+1)-kdeg*mrna
        
        return(list(c(dmrna)))
      })
    }
    
    state <- c(mrna=mrnafirst)
    
    out0 <- ode(y = state, times = seq(0,10800,by=1),
                func = Lorenz5, parms = parameters5)
    out0 <- as.data.frame(out0)
    return(out0)
  }
  
  dfTF <- rbind(D2wt(kdegTF,KDTF),D4wt(kdegTF,KDTF)[2:dim(D4wt(kdegTF,KDTF))[1],])
  dftarget <- rbind(D3wt(kdeg,KD1,KD2),D5wt(kdeg,KD1,KD2,tau)[2:dim(D5wt(kdeg,KD1,KD2,tau))[1],])
  
  x <- dftarget[which(dftarget[,1] %in% seq(0,10800,by=1)),1]
  y <- dftarget[which(dftarget[,1] %in% seq(0,10800,by=1)),2]
  df <- data.frame(time=x,mrna=y)
  
  return(df)
}

D2kd <- function(kdegTF,KDTF) {
  parameters2 <- c(kdegTF=kdegTF,KDTF=KDTF)   
  Lorenz2 <- function(times, state, parameters2){
    with(as.list(c(state, parameters2)), {
      
      dTF <- (((KDTF*inputkd(times))/(1+KDTF*inputkd(times))))-kdegTF*TF
      
      return(list(c(dTF)))
    })
  }
  
  TF0  <- (((KDTF*inputkd(0))/(1+KDTF*inputkd(0))))/kdegTF
  
  state <- c(TF=TF0)
  
  out0 <- ode(y = state, times = seq(-3600,0,by=1),
              func = Lorenz2, parms = parameters2)
  out0 <- as.data.frame(out0)
  return(out0)
}

D4kd <- function(kdegTF,KDTF) {
  parameters4 <- c(kdegTF=kdegTF,KDTF=KDTF)   
  
  TFfirst <- D2kd(kdegTF,KDTF)[3601,2]
  
  Lorenz4 <- function(times, state, parameters4){
    with(as.list(c(state, parameters4)), {
      
      dTF <- (((KDTF*inputkd(times))/(1+KDTF*inputkd(times))))-kdegTF*TF
      
      return(list(c(dTF)))
      
    })
  }
  
  state <- c(TF=TFfirst)
  
  out0 <- ode(y = state, times = seq(0,10800,by=1),
              func = Lorenz4, parms = parameters4)
  out0 <- as.data.frame(out0)
  return(out0)
}



D3kd <- function(kdeg,KD1,KD2) {
  parameters3 <- c(kdeg=kdeg,KD1=KD1,KD2=KD2)  
  
  signal2 <- data.frame(time = D2kd(kdegTF,KDTF)[,1],tf = D2kd(kdegTF,KDTF)[,2])
  sigimp2 <- approxfun(signal2$time, signal2$tf, rule = 2)
  input2kd <- sigimp2
  
  Lorenz3 <- function(times, state, parameters3){
    with(as.list(c(state, parameters3)), {
      
      dmrna = ((KD1*inputkd(times))^2)/(((KD1*inputkd(times))^2)+((KD2*input2kd(times))^2)+1)-kdeg*mrna
      
      return(list(c(dmrna)))
    })
  }
  
  mrna0 = ((KD1*inputkd(0))^2)/(((KD1*inputkd(0))^2)+((KD2*input2kd(0))^2)+1)/kdeg
  
  state <- c(mrna=mrna0)
  
  out0 <- ode(y = state, times = seq(-3600,0,by=1),
              func = Lorenz3, parms = parameters3)
  out0 <- as.data.frame(out0)
  return(out0)
}

D5kd <- function(kdeg,KD1,KD2,tau) {
  parameters5 <- c(kdeg=kdeg,KD1=KD1,KD2=KD2,tau=tau) 
  
  signal3 <- data.frame(time = D4kd(kdegTF,KDTF)[,1],tf = D4kd(kdegTF,KDTF)[,2])
  sigimp3 <- approxfun(signal3$time, signal3$tf, rule = 2)
  input3kd <- sigimp3
  
  mrnafirst <- D3kd(kdeg,KD1,KD2)[3601,2]
  
  Lorenz5 <- function(times, state, parameters5){
    with(as.list(c(state, parameters5)), {
      
      dmrna = ((KD1*inputkd(times-tau))^2)/(((KD1*inputkd(times-tau))^2)+((KD2*input3kd(times-7200))^2)+1)-kdeg*mrna
      
      return(list(c(dmrna)))
    })
  }
  
  state <- c(mrna=mrnafirst)
  
  out0 <- ode(y = state, times = seq(0,10800,by=1),
              func = Lorenz5, parms = parameters5)
  out0 <- as.data.frame(out0)
  return(out0)
}

D1kd <- function(kdeg,KD1,KD2,tau) {
  parameters <- c(kdeg=kdeg,KD1=KD1,KD2=KD2,tau=tau)
  
  D3kd <- function(kdeg,KD1,KD2) {
    parameters3 <- c(kdeg=kdeg,KD1=KD1,KD2=KD2)  
    
    signal2 <- data.frame(time = D2kd(kdegTF,KDTF)[,1],tf = D2kd(kdegTF,KDTF)[,2])
    sigimp2 <- approxfun(signal2$time, signal2$tf, rule = 2)
    input2kd <- sigimp2
    
    Lorenz3 <- function(times, state, parameters3){
      with(as.list(c(state, parameters3)), {
        
        dmrna = ((KD1*inputkd(times))^2)/(((KD1*inputkd(times))^2)+((KD2*input2kd(times))^2)+1)-kdeg*mrna
        
        return(list(c(dmrna)))
      })
    }
    
    mrna0 = ((KD1*inputkd(0))^2)/(((KD1*inputkd(0))^2)+((KD2*input2kd(0))^2)+1)/kdeg
    
    state <- c(mrna=mrna0)
    
    out0 <- ode(y = state, times = seq(-3600,0,by=1),
                func = Lorenz3, parms = parameters3)
    out0 <- as.data.frame(out0)
    return(out0)
  }
  
  D5kd <- function(kdeg,KD1,KD2,tau) {
    parameters5 <- c(kdeg=kdeg,KD1=KD1,KD2=KD2,tau=tau) 
    
    signal3 <- data.frame(time = D4kd(kdegTF,KDTF)[,1],tf = D4kd(kdegTF,KDTF)[,2])
    sigimp3 <- approxfun(signal3$time, signal3$tf, rule = 2)
    input3kd <- sigimp3
    
    mrnafirst <- D3kd(kdeg,KD1,KD2)[3601,2]
    
    Lorenz5 <- function(times, state, parameters5){
      with(as.list(c(state, parameters5)), {
        
        dmrna = ((KD1*inputkd(times-tau))^2)/(((KD1*inputkd(times-tau))^2)+((KD2*input3kd(times-7200))^2)+1)-kdeg*mrna
        
        return(list(c(dmrna)))
      })
    }
    
    state <- c(mrna=mrnafirst)
    
    out0 <- ode(y = state, times = seq(0,10800,by=1),
                func = Lorenz5, parms = parameters5)
    out0 <- as.data.frame(out0)
    return(out0)
  }
  
  dfTF <- rbind(D2kd(kdegTF,KDTF),D4kd(kdegTF,KDTF)[2:dim(D4kd(kdegTF,KDTF))[1],])
  dftarget <- rbind(D3kd(kdeg,KD1,KD2),D5kd(kdeg,KD1,KD2,tau)[2:dim(D5kd(kdeg,KD1,KD2,tau))[1],])
  
  x <- dftarget[which(dftarget[,1] %in% seq(0,10800,by=1)),1]
  y <- dftarget[which(dftarget[,1] %in% seq(0,10800,by=1)),2]
  df <- data.frame(time=x,mrna=y)
  
  return(df)
}
smallest <- read.table("../Data/IFFL_ERGsmallest.txt", sep = "\t",skip=1)
smallest <- read.table("../Data/IFFL_IRGsmallest.txt", sep = "\t",skip=1)
smallest <- read.table("../Data/IFFL_DRGsmallest.txt", sep = "\t",skip=1)
names <- smallest[,1]
smallest <- smallest[,2:3]

refall <- c()
for(r in 1:length(names)[1]){
  d <- names[r]
  
  #stackkari <- c()
  #stack <- c()
  #dfde <- c()
  #file <- paste("../Optimized_parameters/IFFL_model/ERGsubcluster2/rep1/2_IFFL_",d,"_subplex_rep1",".txt",sep="")
  #datafile <- read.table(file,sep="\t",skip=1) 
  #datafile <- as.data.frame(datafile)
  #dfde <- datafile[,2:6]
  #dfde <- as.data.frame(dfde)
  #colnames(dfde) <- c("RMSD","kdeg","kd1","kd2","tau")
  #
  #kdeg <- dfde[smallest[r,1],2]
  #KD1 <- dfde[smallest[r,1],3]
  #KD2 <- dfde[smallest[r,1],4]
  #tau <- dfde[smallest[r,1],5]
  #rmsd <- dfde[smallest[r,1],1]
#
  #count <- read.csv("../Data/160712_ikbakd.csv")
  #count1=count[,2:12]
  #count1 <- data.frame(count1)
  #mean <- count1[,2]
  #condition <- count1[,1]
  #time <- count1[,11]
  #df <- data.frame(Time=time,Mean=mean,Con=condition)
  #df1 <- df[1:13,]
  ##df2 <- df[26:38,]
  #x <- df1[,1]                                         
  #y <- df1[,2] 
  #signal <- data.frame(time = x,nfkb = y)
  #wtsignal <- signal
  #wtsignal$group <- rep("wt",dim(wtsignal)[1])
  #
  #df2 <- df[26:38,]
  #x <- df2[,1]
  #y <- df2[,2] 
  #signal <- data.frame(time = x,nfkb = y)
  #kdsignal <- signal
  #kdsignal$group <- rep("kd",dim(wtsignal)[1]) 
  #
  #wtb <- wtsignal[1,2]
  #for(i in 1:13){
  #  wtsignal[i,2] <- wtsignal[i,2]/wtb
  #}
  #
  #kdb <- kdsignal[1,2]
  #for(i in 1:13){
  #  kdsignal[i,2] <- kdsignal[i,2]/kdb
  #}
  #
  #bothsignal <- rbind(wtsignal,kdsignal)
  #bothsignal <- as.data.frame(bothsignal)
  #
  #df <- rescale(bothsignal[,2], to=c(2,100))
  #df <- data.frame(time=rep(seq(0,10800,900),2),mrna=df)
  #df <- df[1:13,]                                         
  #
  #fp <- pchipfun(df[,1], df[,2]) 
  #f <- fp(seq(0, 10800, by = 1))
  #inputwt <- f
  #df <- data.frame(time=seq(0,10800,1),mrna=f)
  #signal <- data.frame(time = df[,1],nfkb = df[,2])
  #sigimp <- approxfun(signal$t, signal$nfkb, rule = 2)
  #inputwt <- sigimp
  #
  #bothsignal <- rbind(wtsignal,kdsignal)
  #bothsignal <- as.data.frame(bothsignal)
  #
  #df <- rescale(bothsignal[,2], to=c(2,100))
  #df <- data.frame(time=rep(seq(0,10800,900),2),mrna=df)
  #df <- df[14:26,]                                         
  #
  #fp <- pchipfun(df[,1], df[,2]) 
  #f <- fp(seq(0, 10800, by = 1))
  #inputkd <- f
  #df <- data.frame(time=seq(0,10800,1),mrna=f)
  #signal <- data.frame(time = df[,1],nfkb = df[,2])
  #sigimp <- approxfun(signal$t, signal$nfkb, rule = 2)
  #inputkd <- sigimp
  #
  #normwt <- data.frame(x=seq(0,10800,1),y=D1wt(kdeg,KD1,KD2,tau)[,2])
  #minnorm <- normwt[1,2]
  #for(i in 1:dim(normwt)[1]){
  #  normwt[i,2] <- normwt[i,2]/minnorm
  #}
  #
  #normkd <- data.frame(x=seq(0,10800,1),y=D1kd(kdeg,KD1,KD2,tau)[,2])
  #minnorm <- normkd[1,2]
  #for(i in 1:dim(normkd)[1]){
  #  normkd[i,2] <- normkd[i,2]/minnorm
  #}
  #
  #
  #normdf <- rbind(normwt,normkd)
  #normdf <- as.data.frame(normdf)
  #
  #normdf <- rescale(normdf[,2], to=c(0,1))
  #norm <- c()
  #norm <- data.frame(x=rep(seq(0,10800,1),2),y=normdf)
  #
  #alloutswt <- norm[1:10801,]   
  #alloutskd <- norm[10802:21602,]
  
 
  
  stackkari <- c()
  stack <- c()
  dfde <- c()
  file <- paste("../Optimized_parameters/IFFL_model/ERGsubcluster2/rep2/2_IFFL_",d,"_subplex_rep2",".txt",sep="")
  datafile <- read.table(file,sep="\t",skip=1) 
  datafile <- as.data.frame(datafile)
  dfde <- datafile[,2:6]
  dfde <- as.data.frame(dfde)
  colnames(dfde) <- c("RMSD","kdeg","kd1","kd2","tau")
  
  kdeg <- dfde[smallest[r,2],2]
  KD1 <- dfde[smallest[r,2],3]
  KD2 <- dfde[smallest[r,2],4]
  tau <- dfde[smallest[r,2],5]
  rmsd <- dfde[smallest[r,2],1]
  
  count <- read.csv("../Data/160617_ikbakd.csv")
  count1=count[,2:12]
  count1 <- data.frame(count1)
  mean <- count1[,2]
  condition <- count1[,1]
  time <- count1[,11]
  df <- data.frame(Time=time,Mean=mean,Con=condition)
  df1 <- df[1:13,]
  #df2 <- df[26:38,]
  x <- df1[,1]                                         
  y <- df1[,2] 
  signal <- data.frame(time = x,nfkb = y)
  wtsignal <- signal
  wtsignal$group <- rep("wt",dim(wtsignal)[1])
  
  df2 <- df[26:38,]
  x <- df2[,1]
  y <- df2[,2] 
  signal <- data.frame(time = x,nfkb = y)
  kdsignal <- signal
  kdsignal$group <- rep("kd",dim(wtsignal)[1]) 
  
  wtb <- wtsignal[1,2]
  for(i in 1:13){
    wtsignal[i,2] <- wtsignal[i,2]/wtb
  }
  
  kdb <- kdsignal[1,2]
  for(i in 1:13){
    kdsignal[i,2] <- kdsignal[i,2]/kdb
  }
  
  df <- rescale(wtsignal[1:4,2], to=c(2,100))
  unit <- (wtsignal[4,2]-wtsignal[1,2])/(df[4]-df[1])
  wts <- c()
  for(i in 5:13){
    kari <- wtsignal[i,2]-wtsignal[1,2]
    final <- kari/unit+2
    wts[i] <- final
  }
  wts <- na.omit(wts)
  wts <- c(df,wts)
  
  kds <- c()
  for(i in 1:dim(kdsignal)[1]){
    kari <- kdsignal[i,2]-wtsignal[1,2]
    final <- kari/unit+2
    kds[i] <- final
  }
  
  df <- data.frame(time=seq(0,10800,900),mrna=wts)
  
  fp <- pchipfun(df[,1], df[,2]) 
  f <- fp(seq(0, 10800, by = 1))
  inputwt <- f
  df <- data.frame(time=seq(0,10800,1),mrna=f)
  signal <- data.frame(time = df[,1],nfkb = df[,2])
  sigimp <- approxfun(signal$t, signal$nfkb, rule = 2)
  inputwt <- sigimp
  
  df <- data.frame(time=seq(0,10800,900),mrna=kds)
  
  fp <- pchipfun(df[,1], df[,2]) 
  f <- fp(seq(0, 10800, by = 1))
  inputkd <- f
  df <- data.frame(time=seq(0,10800,1),mrna=f)
  signal <- data.frame(time = df[,1],nfkb = df[,2])
  sigimp <- approxfun(signal$t, signal$nfkb, rule = 2)
  inputkd <- sigimp
  
  normwt <- data.frame(x=seq(0,10800,1),y=D1wt(kdeg,KD1,KD2,tau)[,2])
  minnorm <- normwt[1,2]
  for(i in 1:dim(normwt)[1]){
    normwt[i,2] <- normwt[i,2]/minnorm
  }
  
  normkd <- data.frame(x=seq(0,10800,1),y=D1kd(kdeg,KD1,KD2,tau)[,2])
  minnorm <- normkd[1,2]
  for(i in 1:dim(normkd)[1]){
    normkd[i,2] <- normkd[i,2]/minnorm
  }
  
  
  normdf <- rbind(normwt,normkd)
  normdf <- as.data.frame(normdf)
  
  normdf <- rescale(normdf[,2], to=c(0,1))
  norm <- c()
  norm <- data.frame(x=rep(seq(0,10800,1),2),y=normdf)
  
  alloutswt <- norm[1:10801,]   
  alloutskd <- norm[10802:21602,]
  
  
  allout <- rbind(alloutswt,alloutskd)
  allout <- as.data.frame(allout)
  
  outfinal2 <- allout[which(allout[,1] %in% seq(0,10800,by=900)),]
  outfinal2 <- outfinal2[,2]
  outfinal2 <- as.data.frame(outfinal2)
  outfinal2 <- t(outfinal2)
  outfinal2 <- as.data.frame(outfinal2)
  
  refall <- c(refall, c(outfinal2))
}

refall2 <- data.frame(matrix(unlist(refall),ncol=26,byrow=T))

refall3 <- refall2[,1:13]    #siCtrl
#refall3 <- refall2[,14:26]   #siIkBa

refall3 <- t(refall3)
refall3 <- as.data.frame(refall3)
colnames(refall3) <- names

p <- colnames(refall3)
t <- colnames(refall3)
refall3$t <- seq(0,180,15)
temp <- reshape2::melt(refall3,
                       id="t",
                       measure=p 
)
temp$variable <- as.character(temp$variable)
temp$value <- as.numeric(temp$value)

g <- ggplot(temp,aes(as.factor(t),as.factor(variable)))+
  geom_tile(aes(fill=value),color="grey50",stat="identity")+
  scale_fill_gradient(low="white",high="magenta",limits=c(0,1.0))+
  #scale_fill_gradient(low="white",high="mediumseagreen",limits=c(0,1.0))+
  labs(x = "",y = "") + 
  theme(panel.grid = element_blank())+#
  scale_x_discrete(labels=c("0","15","30","45","60","75","90","105","120","135","150","165","180"),breaks=c("0","15","30","45","60","75","90","105","120","135","150","165","180"),expand = c(0, 0))+ 
  scale_y_discrete(breaks=sort(names,decreasing = T),labels=waiver(),limits=sort(names,decreasing = T),expand = c(0, 0))+
  theme(legend.position = "bottom",axis.ticks = element_blank(), axis.text.x = element_text(size = 7, angle = 0, hjust = 0.5, colour = "black"),axis.text.y = element_text(size = 6.0, angle = 0, hjust = 1.0, colour = "black"));g

#Draw a heatmap of nRMSD
names <- read.table("../Data/ERGsubcluster2_Promoter_up500down500TSS_12012019.txt", sep = "\t",header=FALSE)
names <- names[-which(names[,1]=="NFKBIA"),]
names <- read.table("../Data/IRGsubcluster2_Promoter_up500down500TSS_12012019.txt", sep = "\t",header=FALSE)
names[,1] <- as.character(names[,1])    
names[which(names[,1]=="NKX3_1"),1] <- c("NKX3-1")
names <- read.table("../Data/DRGsubcluster2_Promoter_up500down500TSS_12012019.txt", sep = "\t",header=FALSE)
names <- names[,1]
names <- sort(names)


df <- data.frame(name=names,group=rep("black",length(names)))
df$name <- as.character(df$name)
df$group <- as.character(df$group)
df[which(df$name=="EFNA1"),2]="yellow"
df[which(df$name=="IRF1"),2]="yellow"
df[which(df$name=="PLEKHF2"),2]="yellow"
df[which(df$name=="SOX9"),2]="yellow"
df[which(df$name=="SPRY4"),2]="yellow"
#df[which(df$name=="ANKRD33B"),2]="yellow"
#df[which(df$name=="BID"),2]="yellow"
#df[which(df$name=="CLIC4"),2]="yellow"
#df[which(df$name=="NFKBIE"),2]="yellow"
#df[which(df$name=="SDC4"),2]="yellow"
#simple_genes <- c("ABTB2","ARID5B","CREBZF","CRIM1","EFEMP1","FAM117A","HIVEP2","IFNGR2","IL6ST","KLHL5","KYNU","LAMC2","MGAT4A","NFATC2","NFKB1","PHLDB2","PRKAR2B","RAB3IP","RBM25","SEMA3C","SGPP2","STAT5A","TRIP10")
df[which(df$name %in% simple_genes),2]="yellow"
df$x <- rep("all",dim(df)[1])

g <- ggplot(df,aes(as.factor(x),as.factor(name)))+
  geom_tile(aes(fill=group),color="grey50",stat="identity",na.rm = TRUE)+
  scale_fill_manual(values=c(black="grey85", yellow="yellow"))+
  labs(x = "",y = "") + 
  theme(panel.grid = element_blank())+
  scale_x_discrete(labels=waiver(),breaks=group,expand = c(0, 0))+ 
  scale_y_discrete(breaks=sort(df$name,decreasing=T),labels=waiver(),limits=sort(df$name,decreasing=T),expand = c(0, 0))+
  theme(legend.position = "right",axis.ticks = element_blank(), axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, colour = "black"),axis.text.y = element_text(size = 12, angle = 0, hjust = 0.5, colour = "black"));g

#Draw lineplots of example gene
smallest <- read.table("../Data/IFFL_ERGsmallest.txt", sep = "\t",skip=1)
smallest <- read.table("../Data/IFFL_IRGsmallest.txt", sep = "\t",skip=1)
smallest <- read.table("../Data/IFFL_DRGsmallest.txt", sep = "\t",skip=1)
names <- smallest[,1]
smallest <- smallest[,2:3]

d <- names[r]

#stackkari <- c()
#stack <- c()
#dfde <- c()
#file <- paste("../Optimized_parameters/IFFL_model/ERGsubcluster2/rep1/2_IFFL_",d,"_subplex_rep1",".txt",sep="")
#datafile <- read.table(file,sep="\t",skip=1) 
#datafile <- as.data.frame(datafile)
#dfde <- datafile[,2:6]
#dfde <- as.data.frame(dfde)
#colnames(dfde) <- c("RMSD","kdeg","kd1","kd2","tau")
#
#kdeg <- dfde[smallest[r,1],2]
#KD1 <- dfde[smallest[r,1],3]
#KD2 <- dfde[smallest[r,1],4]
#tau <- dfde[smallest[r,1],5]
#rmsd <- dfde[smallest[r,1],1]

stackkari <- c()
stack <- c()
dfde <- c()
file <- paste("../Optimized_parameters/IFFL_model/ERGsubcluster2/rep2/2_IFFL_",d,"_subplex_rep2",".txt",sep="")
datafile <- read.table(file,sep="\t",skip=1) 
datafile <- as.data.frame(datafile)
dfde <- datafile[,2:6]
dfde <- as.data.frame(dfde)
colnames(dfde) <- c("RMSD","kdeg","kd1","kd2","tau")

kdeg <- dfde[smallest[r,2],2]
KD1 <- dfde[smallest[r,2],3]
KD2 <- dfde[smallest[r,2],4]
tau <- dfde[smallest[r,2],5]
rmsd <- dfde[smallest[r,2],1]

#count <- read.csv("../Data/160712_ikbakd.csv")
#count1=count[,2:12]
#count1 <- data.frame(count1)
#mean <- count1[,2]
#condition <- count1[,1]
#time <- count1[,11]
#df <- data.frame(Time=time,Mean=mean,Con=condition)
#df1 <- df[1:13,]
##df2 <- df[26:38,]
#x <- df1[,1]                                         
#y <- df1[,2] 
#signal <- data.frame(time = x,nfkb = y)
#wtsignal <- signal
#wtsignal$group <- rep("wt",dim(wtsignal)[1])
#
#df2 <- df[26:38,]
#x <- df2[,1]
#y <- df2[,2] 
#signal <- data.frame(time = x,nfkb = y)
#kdsignal <- signal
#kdsignal$group <- rep("kd",dim(wtsignal)[1]) 
#
#
#wtb <- wtsignal[1,2]
#for(i in 1:13){
#  wtsignal[i,2] <- wtsignal[i,2]/wtb
#}
#
#kdb <- kdsignal[1,2]
#for(i in 1:13){
#  kdsignal[i,2] <- kdsignal[i,2]/kdb
#}
#
#
#bothsignal <- rbind(wtsignal,kdsignal)
#bothsignal <- as.data.frame(bothsignal)
#
#df <- rescale(bothsignal[,2], to=c(2,100))
#df <- data.frame(time=rep(seq(0,10800,900),2),mrna=df)
#df <- df[1:13,]                                         
#
#fp <- pchipfun(df[,1], df[,2]) 
#f <- fp(seq(0, 10800, by = 1))
#inputwt <- f
#df <- data.frame(time=seq(0,10800,1),mrna=f)
#signal <- data.frame(time = df[,1],nfkb = df[,2])
#sigimp <- approxfun(signal$t, signal$nfkb, rule = 2)
#inputwt <- sigimp
#nfkb_wt <- signal
#colnames(nfkb_wt) <- c("Time","Value")
#
#times <- seq(-3600,0,1)
#df <- rep(nfkb_wt[1,2],length(times))
#fp <- pchipfun(times, df) 
#sigimp <- approxfun(times, df, rule = 2)
#steadywt <- sigimp
#
#bothsignal <- rbind(wtsignal,kdsignal)
#bothsignal <- as.data.frame(bothsignal)
#
#df <- rescale(bothsignal[,2], to=c(2,100))
#df <- data.frame(time=rep(seq(0,10800,900),2),mrna=df)
#df <- df[14:26,]                                         
#
#fp <- pchipfun(df[,1], df[,2]) 
#f <- fp(seq(0, 10800, by = 1))
#inputkd <- f
#df <- data.frame(time=seq(0,10800,1),mrna=f)
#signal <- data.frame(time = df[,1],nfkb = df[,2])
#sigimp <- approxfun(signal$t, signal$nfkb, rule = 2)
#inputkd <- sigimp
#nfkb_kd <- signal
#colnames(nfkb_kd) <- c("Time","Value")
#
#
#times <- seq(-3600,0,1)
#df <- rep(nfkb_kd[1,2],length(times))
#fp <- pchipfun(times, df) 
#sigimp <- approxfun(times, df, rule = 2)
#steadykd <- sigimp



count <- read.csv("../Data/160617_ikbakd.csv")
count1=count[,2:12]
count1 <- data.frame(count1)
mean <- count1[,2]
condition <- count1[,1]
time <- count1[,11]
df <- data.frame(Time=time,Mean=mean,Con=condition)
df1 <- df[1:13,]
#df2 <- df[26:38,]
x <- df1[,1]                                         
y <- df1[,2] 
signal <- data.frame(time = x,nfkb = y)
wtsignal <- signal
wtsignal$group <- rep("wt",dim(wtsignal)[1])

df2 <- df[26:38,]
x <- df2[,1]
y <- df2[,2] 
signal <- data.frame(time = x,nfkb = y)
kdsignal <- signal
kdsignal$group <- rep("kd",dim(wtsignal)[1]) 


wtb <- wtsignal[1,2]
for(i in 1:13){
  wtsignal[i,2] <- wtsignal[i,2]/wtb
}

kdb <- kdsignal[1,2]
for(i in 1:13){
  kdsignal[i,2] <- kdsignal[i,2]/kdb
}

df <- rescale(wtsignal[1:4,2], to=c(2,100))
unit <- (wtsignal[4,2]-wtsignal[1,2])/(df[4]-df[1])
wts <- c()
for(i in 5:13){
  kari <- wtsignal[i,2]-wtsignal[1,2]
  final <- kari/unit+2
  wts[i] <- final
}
wts <- na.omit(wts)
wts <- c(df,wts)

kds <- c()
for(i in 1:dim(kdsignal)[1]){
  kari <- kdsignal[i,2]-wtsignal[1,2]
  final <- kari/unit+2
  kds[i] <- final
}

df <- data.frame(time=seq(0,10800,900),mrna=wts)

fp <- pchipfun(df[,1], df[,2]) 
f <- fp(seq(0, 10800, by = 1))
inputwt <- f
df <- data.frame(time=seq(0,10800,1),mrna=f)
signal <- data.frame(time = df[,1],nfkb = df[,2])
sigimp <- approxfun(signal$t, signal$nfkb, rule = 2)
inputwt <- sigimp
nfkb_wt <- signal
colnames(nfkb_wt) <- c("Time","Value")

times <- seq(-3600,0,1)
df <- rep(nfkb_wt[1,2],length(times))
fp <- pchipfun(times, df) 
sigimp <- approxfun(times, df, rule = 2)
steadywt <- sigimp

df <- data.frame(time=seq(0,10800,900),mrna=kds)

fp <- pchipfun(df[,1], df[,2]) 
f <- fp(seq(0, 10800, by = 1))
inputkd <- f
df <- data.frame(time=seq(0,10800,1),mrna=f)
signal <- data.frame(time = df[,1],nfkb = df[,2])
sigimp <- approxfun(signal$t, signal$nfkb, rule = 2)
inputkd <- sigimp
nfkb_kd <- signal
colnames(nfkb_kd) <- c("Time","Value")

times <- seq(-3600,0,1)
df <- rep(nfkb_kd[1,2],length(times))
fp <- pchipfun(times, df) 
sigimp <- approxfun(times, df, rule = 2)
steadykd <- sigimp


normwt <- data.frame(x=seq(0,10800,1),y=D1wt(kdeg,KD1,KD2,tau)[,2])
minnorm <- normwt[1,2]
for(i in 1:dim(normwt)[1]){
  normwt[i,2] <- normwt[i,2]/minnorm
}

normkd <- data.frame(x=seq(0,10800,1),y=D1kd(kdeg,KD1,KD2,tau)[,2])
minnorm <- normkd[1,2]
for(i in 1:dim(normkd)[1]){
  normkd[i,2] <- normkd[i,2]/minnorm
}


normdf <- rbind(normwt,normkd)
normdf <- as.data.frame(normdf)

normdf <- rescale(normdf[,2], to=c(0,1))
norm <- c()
norm <- data.frame(x=rep(seq(0,10800,1),2),y=normdf)


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

df <- rescale(df[,2], to=c(0,1))
refwtfc <- data.frame(x=seq(0,180,15),y=df[1:13])
refkdfc <- data.frame(x=seq(0,180,15),y=df[14:26])

refwt <- refwtfc[,2]
refkd <- refkdfc[,2]

sim <- norm[which(norm[,1] %in% seq(0,10800,900)),]
simwt <- sim[1:13,2]
simkd <- sim[14:26,2]

alloutswt <- norm[1:10801,]   
alloutskd <- norm[10802:21602,]

alloutswt$condition <- rep("wt",dim(alloutswt)[1])
alloutskd$condition <- rep("kd",dim(alloutskd)[1])
alloutswt$line <- rep("solid",dim(alloutswt)[1])
alloutskd$line <- rep("solid",dim(alloutskd)[1])
refwt <- data.frame(time=seq(0,10800,900),mrna=refwt)
refkd <- data.frame(time=seq(0,10800,900),mrna=refkd)
refwt$condition <- rep("wt",dim(refwt)[1])
refkd$condition <- rep("kd",dim(refkd)[1])
refwt$line <- rep("dashed",dim(refwt)[1])
refkd$line <- rep("dashed",dim(refkd)[1])

colnames(alloutswt) <- c("time","mrna","condition","line")
colnames(alloutskd) <- c("time","mrna","condition","line")


df <- rbind(alloutswt,alloutskd,refwt,refkd)
df <- as.data.frame(df)

g <- ggplot(df,aes(x=time,y=mrna,color=condition,linetype=line,alpha=line))+
  theme_bw(base_size = 16)+
  theme(panel.grid = element_blank())+
  labs(x="",y="",title=c(""),colour="")+
  scale_colour_manual(values = c(wt="magenta",kd="green4"))+
  scale_x_continuous(breaks=seq(0,10800,1800),labels=seq(0,180,30))+
  scale_y_continuous(breaks=c(0.0,0.2,0.4,0.6,0.8,1.0),labels=waiver(),limits=c(0.0,1.0))+
  scale_linetype_manual(values=c(solid="solid", dashed="dotted"))+
  scale_alpha_manual(values=c(solid=0.7, dashed=1.0))+
  theme(panel.border=element_rect(fill=NA,color="black", size=0.6),plot.title = element_text(hjust = 0.5,size=15.0,margin=margin(t = 0, r = 0, b = 0.1, l = 0)),panel.grid.major.x = element_line(colour="grey90", size =0.001),panel.grid.major.y = element_line(colour="grey90", size =0.001),legend.position = "none",title = element_text(size=7),legend.text = element_text(size=7),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+
  geom_line(lwd=2.0,aes(color=condition));g



