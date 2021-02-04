library(edgeR)
library(pracma)
library(deSolve)
library(scales)
library(bio3d)
library(nloptr)


count <- read.csv("../Data/160617_ikbakd.csv")
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
 
  count <- read.table("../Data/counts_all.txt",sep="\t",header= T,row.names = 1,skip = 1)
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
  
  count <- read.table("../Data/counts_KO_all.txt",sep="\t",header= T,row.names = 1,skip = 1)
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
  count <- read.table("../Data/counts_KO_all.txt",sep="\t",header= T,row.names = 1,skip = 1)
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
  countname <- read.table("../Data/counts_KO_all.txt",sep="\t",header= T,row.names = 1,skip = 1)
  rownames(count) <- rownames(countname)
  tpm <- matrix(0,nrow = dim(count)[1],ncol = dim(count)[2]/2)
  count <- count[,2:27]
  for(i in 1:13){
    tpm[,i] <- (count[,2*i-1]+count[,2*i])/2 
  }
  tpm <- as.data.frame(tpm)
  rownames(tpm) <- rownames(count)
  colnames(tpm) <- seq(0,180,15)
  kd <- tpm
  
  wtdeg <- read.table("../Data/WT_DEGs.txt",skip=1,row.names=1,sep="\t",stringsAsFactors = FALSE)
  wtdeg <- wtdeg[,1]
  wtdeg <- as.character(wtdeg)
  
  wt2 <- wt[which(rownames(wt) %in% wtdeg),]
  kd2 <- kd[which(rownames(kd) %in% wtdeg),]
  
  ref <- wt2[which(rownames(wt2)=="gene.name"),]
  
  data <- ref
  zero <- data[1,1]
  for(i in 1:13){
    data[1,i] <- data[1,i]/zero
  }
  
  refwtfc <- data.frame(x=seq(0,180,15),y=as.numeric(data[1,]))
  
  
  ref <- kd2[which(rownames(kd2)=="gene.name"),]
  
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
  
  
  
  ewt <- function(kdeg,KD1,KD2,tau){  
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
    norm <- data.frame(x=seq(0,10800,1),y=normdf[1:10801])
    
    return(norm)
  }
  
  ekd <- function(kdeg,KD1,KD2,tau){  
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
    norm <- data.frame(x=seq(0,10800,1),y=normdf[10802:21602])
    
    
    return(norm)
  }
  
  g <- function(kdeg,KD1,KD2,tau){
    
    rmsdwt <- rmsd(refwtfc[,2],ewt(kdeg,KD1,KD2,tau)[which(ewt(kdeg,KD1,KD2,tau)[,1] %in% seq(0,10800,900)),2],fit=FALSE)
    mul <- rmsdwt
    
    return(mul)
  }
  
  h <- function(p){
    g(p[1], p[2], p[3], p[4])
  }
    
  cm <- read.table("../Data/IFFL_Subplex_random_initial.txt",sep="\t",skip=1)
  cm <- cm[1:50,2:5]
  cm <- as.data.frame(cm)
  
  ps <- c()
  ps0 <- c()
  for(i in 1:dim(cm)[1]){
    pa1 <- cm[i,1]
    pa2 <- cm[i,2]
    pa3 <- cm[i,3]
    pa4 <- cm[i,4]
    
    result <- c()
    result <- sbplx(c(pa1,pa2,pa3,pa4),h,control=list(maxeval=200),lower=c(kdeg=0.00002,KD1=0.001,KD2=0.001,tau=0),upper=c(kdeg=0.002,KD1=1000,KD2=1000,tau=7200))
   
    kdegpar <- c()
    KD1par <- c()
    KD2par <- c()
    taupar <- c()
    rmsd <- c()
    kdegpar <- result$par[1]
    KD1par <- result$par[2]
    KD2par <- result$par[3]
    taupar <- result$par[4]
    rmsd <- result$value
    
    ps0 <- c(rmsd,kdegpar,KD1par,KD2par,taupar)
    
    ps <- c(ps,c(ps0))
  }
  ps <- data.frame(matrix(unlist(ps),ncol=5,byrow=T))

  write.table(ps, "IFFL_gene.name_rep2_1.txt", sep = "\t")
  