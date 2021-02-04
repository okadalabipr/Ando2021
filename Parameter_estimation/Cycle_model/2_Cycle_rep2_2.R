library(edgeR)
library(pracma)
library(deSolve)
library(scales)
library(bio3d)
library(nloptr)
library(ggplot2)

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

wtBase0 <- function(k1,k2,KD1,KD2) { 
  
  pars <- c(k1=k1,k2=k2,KD1=KD1,KD2=KD2)
  
  wtBaseModel0 <- function (times, State, pars) {
    
    pars <- c(k1=k1,k2=k2,KD1=KD1,KD2=KD2)
    
    
    with(as.list(c(State, pars)),{
      
      A = 1-C-O
      
      k_1 <- 0.01*k1
      k_2 <- 0.01*k2
      k_3 <- 1
      
      f1 <- k1*((KD1*steadywt(times))/((KD1*steadywt(times))+1))*C
      f_1 <- k_1 * O
      f_3 <- k_3 *((0.5*A)/((0.5*A)+1))
      
      f2 <- k2*((KD2*steadywt(times))/((KD2*steadywt(times))+1))*O
      f_2 <- k_2 * A
      
      
      dC <- -f1+f_1+f_3
      dO <- f1 - f2 - f_1 + f_2
      
      # return
      return(list(c(dC, dO)))
    })
  }
  
  wtBaseModel <- function (times, State, pars) {
    
    pars <- c(k1=k1,k2=k2,KD1=KD1,KD2=KD2)
    
    
    with(as.list(c(State, pars)),{
      
      A = 1-C-O
      k_1 <- 0.01*k1
      k_2 <- 0.01*k2
      k_3 <- 1
      
      f1 <- k1*((KD1*inputwt(times))/((KD1*inputwt(times))+1))*C
      f_1 <- k_1 * O
      f_3 <- k_3 *((0.5*A)/((0.5*A)+1))
      
      f2 <- k2*((KD2*inputwt(times))/((KD2*inputwt(times))+1))*O
      f_2 <- k_2 * A
      
      
      
      dC <- -f1+f_1+f_3
      dO <- f1 - f2 - f_1 + f_2
      
      return(list(c(dC, dO)))
    })
  }
  
  ss.init <- ode(y=c(C=1,O=0), times=seq(-3600,0,1), wtBaseModel0, parms=pars)
  ss.init <- as.data.frame(ss.init)
  
  
  
  out <- ode(y=c(C=ss.init[3601,2],O=ss.init[3601,3]), times=seq(0,10800,1), wtBaseModel, parms=pars)
  
  out <- as.data.frame(out)
  out$A <- 1-(out$C+out$O) 
  
  out
  
  return(out)
}

wtModel0 <- function(kdeg) { 
  
  
  signal2 <- data.frame(time = wtBase0(k1,k2,KD1,KD2)[,1],C = wtBase0(k1,k2,KD1,KD2)[,2],O = wtBase0(k1,k2,KD1,KD2)[,3])
  signal2$A <- 1-(signal2$C+signal2$O)
  sigimp2 <- approxfun(signal2$time, signal2$A, rule = 2)
  active <- sigimp2
  
  mrnafirst <- (1*active(0))/kdeg
  
  odeModel <- function (times, state, pars) {
    with(as.list(c(state, pars)),{
      
      dmrna <- (1*active(times))-(kdeg*mrna)
      
      return(list(c(dmrna)))
    })
  }
  
  state <- c(mrna=mrnafirst)
  
  out0 <- ode(y = state, times = seq(-3600,0,by=1),
              func = odeModel, parms = pars)
  out0 <- as.data.frame(out0)
  return(out0)
}

wtModel <- function(kdeg,tau) { 
  
  
  signal2 <- data.frame(time = wtBase0(k1,k2,KD1,KD2)[,1],C = wtBase0(k1,k2,KD1,KD2)[,2],O = wtBase0(k1,k2,KD1,KD2)[,3])
  signal2$A <- 1-(signal2$C+signal2$O)
  sigimp2 <- approxfun(signal2$time, signal2$A, rule = 2)
  active <- sigimp2
  
  mrnafirst <- wtModel0(kdeg)[3601,2]
  
  odeModel <- function (times, state, pars) {
    with(as.list(c(state, pars)),{
      
      dmrna <- (1*active(times-tau))-(kdeg*mrna)
      
      return(list(c(dmrna)))
    })
  }
  
  state <- c(mrna=mrnafirst)
  
  out0 <- ode(y = state, times = seq(0,10800,by=1),
              func = odeModel, parms = pars)
  out0 <- as.data.frame(out0)
  return(out0)
}

wtModels <- function(k1,k2,KD1,KD2,kdeg,tau) { 
  
  pars <- c(k1=k1,k2=k2,KD1=KD1,KD2=KD2,kdeg=kdeg,tau=tau)
  
  wtBase0 <- function(k1,k2,KD1,KD2) { 
    
    pars <- c(k1=k1,k2=k2,KD1=KD1,KD2=KD2)
    
    wtBaseModel0 <- function (times, State, pars) {
      
      pars <- c(k1=k1,k2=k2,KD1=KD1,KD2=KD2)
      
      
      with(as.list(c(State, pars)),{
        
        A = 1-C-O
        
        k_1 <- 0.01*k1
        k_2 <- 0.01*k2
        k_3 <- 1
        
        f1 <- k1*((KD1*steadywt(times))/((KD1*steadywt(times))+1))*C
        f_1 <- k_1 * O
        f_3 <- k_3 *((0.5*A)/((0.5*A)+1))
        
        f2 <- k2*((KD2*steadywt(times))/((KD2*steadywt(times))+1))*O
        f_2 <- k_2 * A
        
        
        
        dC <- -f1+f_1+f_3
        
        dO <- f1 - f2 - f_1 + f_2
        
        # return
        return(list(c(dC, dO)))
      })
    }
    
    wtBaseModel <- function (times, State, pars) {
      
      pars <- c(k1=k1,k2=k2,KD1=KD1,KD2=KD2)
      
      
      with(as.list(c(State, pars)),{
        
        A = 1-C-O
        k_1 <- 0.01*k1
        k_2 <- 0.01*k2
        k_3 <- 1  
        
        f1 <- k1*((KD1*inputwt(times))/((KD1*inputwt(times))+1))*C
        f_1 <- k_1 * O
        f_3 <- k_3 *((0.5*A)/((0.5*A)+1))
        
        f2 <- k2*((KD2*inputwt(times))/((KD2*inputwt(times))+1))*O
        f_2 <- k_2 * A
        
        
        
        
        
        dC <- -f1+f_1+f_3
        dO <- f1 - f2 - f_1 + f_2
        
        # return
        return(list(c(dC, dO)))
      })
    }
    
    
    ss.init <- ode(y=c(C=1,O=0), times=seq(-3600,0,1), wtBaseModel0, parms=pars)
    ss.init <- as.data.frame(ss.init)
    
    
    out <- ode(y=c(C=ss.init[3601,2],O=ss.init[3601,3]), times=seq(0,10800,1), wtBaseModel, parms=pars)
    
    out <- as.data.frame(out)
    out$A <- 1-(out$C+out$O) 
    
    out
    
    return(out)
  }
  
  wtModel0 <- function(kdeg) { 
    
    
    signal2 <- data.frame(time = wtBase0(k1,k2,KD1,KD2)[,1],C = wtBase0(k1,k2,KD1,KD2)[,2],O = wtBase0(k1,k2,KD1,KD2)[,3])
    signal2$A <- 1-(signal2$C+signal2$O)
    sigimp2 <- approxfun(signal2$time, signal2$A, rule = 2)
    active <- sigimp2
    
    mrnafirst <- (1*active(0))/kdeg
    
    odeModel <- function (times, state, pars) {
      with(as.list(c(state, pars)),{
        
        dmrna <- (1*active(times))-(kdeg*mrna)
        
        return(list(c(dmrna)))
      })
    }
    
    state <- c(mrna=mrnafirst)
    
    out0 <- ode(y = state, times = seq(-3600,0,by=1),
                func = odeModel, parms = pars)
    out0 <- as.data.frame(out0)
    return(out0)
  }
  
  wtModel <- function(kdeg,tau) { 
    
    
    
    signal2 <- data.frame(time = wtBase0(k1,k2,KD1,KD2)[,1],C = wtBase0(k1,k2,KD1,KD2)[,2],O = wtBase0(k1,k2,KD1,KD2)[,3])
    signal2$A <- 1-(signal2$C+signal2$O)
    sigimp2 <- approxfun(signal2$time, signal2$A, rule = 2)
    active <- sigimp2
    
    mrnafirst <- wtModel0(kdeg)[3601,2]
    
    odeModel <- function (times, state, pars) {
      with(as.list(c(state, pars)),{
        
        dmrna <- (1*active(times-tau))-(kdeg*mrna)
        
        return(list(c(dmrna)))
      })
    }
    
    state <- c(mrna=mrnafirst)
    
    out0 <- ode(y = state, times = seq(0,10800,by=1),
                func = odeModel, parms = pars)
    out0 <- as.data.frame(out0)
    return(out0)
  }
  
  
  
  dftarget <- wtModel(kdeg,tau)
  
  x <- dftarget[which(dftarget[,1] %in% seq(0,10800,by=1)),1]
  y <- dftarget[which(dftarget[,1] %in% seq(0,10800,by=1)),2]
  df <- data.frame(time=x,mrna=y)
  
  return(df)
}


kdBase0 <- function(k1,k2,KD1,KD2) { 
  
  pars <- c(k1=k1,k2=k2,KD1=KD1,KD2=KD2)
  
  kdBaseModel0 <- function (times, State, pars) {
    
    pars <- c(k1=k1,k2=k2,KD1=KD1,KD2=KD2)
    
    
    with(as.list(c(State, pars)),{
      
      A = 1-C-O
      
      k_1 <- 0.01*k1
      k_2 <- 0.01*k2
      k_3 <- 1 
      
      f1 <- k1*((KD1*steadykd(times))/((KD1*steadykd(times))+1))*C
      f_1 <- k_1 * O
      f_3 <- k_3 *((0.5*A)/((0.5*A)+1))
      
      f2 <- k2*((KD2*steadykd(times))/((KD2*steadykd(times))+1))*O
      f_2 <- k_2 * A
      
      
      dC <- -f1+f_1+f_3
      dO <- f1 - f2 - f_1 + f_2
      
      return(list(c(dC, dO)))
    })
  }
  
  kdBaseModel <- function (times, State, pars) {
    
    pars <- c(k1=k1,k2=k2,KD1=KD1,KD2=KD2)
    
    
    with(as.list(c(State, pars)),{
      
      A = 1-C-O
      
      k_1 <- 0.01*k1
      k_2 <- 0.01*k2
      k_3 <- 1
      
      f1 <- k1*((KD1*inputkd(times))/((KD1*inputkd(times))+1))*C
      f_1 <- k_1 * O
      f_3 <- k_3 *((0.5*A)/((0.5*A)+1))
      
      f2 <- k2*((KD2*inputkd(times))/((KD2*inputkd(times))+1))*O
      f_2 <- k_2 * A
      
      
      
      dC <- -f1+f_1+f_3
      dO <- f1 - f2 - f_1 + f_2
      
      return(list(c(dC, dO)))
    })
  }
  
  ss.init <- ode(y=c(C=1,O=0), times=seq(-3600,0,1), kdBaseModel0, parms=pars)
  ss.init <- as.data.frame(ss.init)
  
  out <- ode(y=c(C=ss.init[3601,2],O=ss.init[3601,3]), times=seq(0,10800,1), kdBaseModel, parms=pars)
  
  out <- as.data.frame(out)
  out$A <- 1-(out$C+out$O) 
  
  out
  
  return(out)
}

kdModel0 <- function(kdeg) { 
  
  
  signal2 <- data.frame(time = kdBase0(k1,k2,KD1,KD2)[,1],C = kdBase0(k1,k2,KD1,KD2)[,2],O = kdBase0(k1,k2,KD1,KD2)[,3])
  signal2$A <- 1-(signal2$C+signal2$O)
  sigimp2 <- approxfun(signal2$time, signal2$A, rule = 2)
  active <- sigimp2
  
  mrnafirst <- (1*active(0))/kdeg
  
  odeModel <- function (times, state, pars) {
    with(as.list(c(state, pars)),{
      
      dmrna <- (1*active(times))-(kdeg*mrna)
      
      return(list(c(dmrna)))
    })
  }
  
  state <- c(mrna=mrnafirst)
  
  out0 <- ode(y = state, times = seq(-3600,0,by=1),
              func = odeModel, parms = pars)
  out0 <- as.data.frame(out0)
  return(out0)
}

kdModel <- function(kdeg,tau) { 
  
  
  
  signal2 <- data.frame(time = kdBase0(k1,k2,KD1,KD2)[,1],C = kdBase0(k1,k2,KD1,KD2)[,2],O = kdBase0(k1,k2,KD1,KD2)[,3])
  signal2$A <- 1-(signal2$C+signal2$O)
  sigimp2 <- approxfun(signal2$time, signal2$A, rule = 2)
  active <- sigimp2
  
  mrnafirst <- kdModel0(kdeg)[3601,2]
  
  odeModel <- function (times, state, pars) {
    with(as.list(c(state, pars)),{
      
      dmrna <- (1*active(times-tau))-(kdeg*mrna)
      
      return(list(c(dmrna)))
    })
  }
  
  state <- c(mrna=mrnafirst)
  
  out0 <- ode(y = state, times = seq(0,10800,by=1),
              func = odeModel, parms = pars)
  out0 <- as.data.frame(out0)
  return(out0)
}

kdModels <- function(k1,k2,KD1,KD2,kdeg,tau) { 
  
  pars <- c(k1=k1,k2=k2,KD1=KD1,KD2=KD2,kdeg=kdeg,tau=tau)
  
  kdBase0 <- function(k1,k2,KD1,KD2) { 
    
    pars <- c(k1=k1,k2=k2,KD1=KD1,KD2=KD2)
    
    kdBaseModel0 <- function (times, State, pars) {
      
      pars <- c(k1=k1,k2=k2,KD1=KD1,KD2=KD2)
      
      
      with(as.list(c(State, pars)),{
        
        A = 1-C-O
        
        k_1 <- 0.01*k1
        k_2 <- 0.01*k2
        k_3 <- 1
        
        f1 <- k1*((KD1*steadykd(times))/((KD1*steadykd(times))+1))*C
        f_1 <- k_1 * O
        f_3 <- k_3 *((0.5*A)/((0.5*A)+1))
        
        f2 <- k2*((KD2*steadykd(times))/((KD2*steadykd(times))+1))*O
        f_2 <- k_2 * A
        
        
        dC <- -f1+f_1+f_3
        dO <- f1 - f2 - f_1 + f_2
        
        return(list(c(dC, dO)))
      })
    }
    
    kdBaseModel <- function (times, State, pars) {
      
      pars <- c(k1=k1,k2=k2,KD1=KD1,KD2=KD2)
      
      with(as.list(c(State, pars)),{
        
        A = 1-C-O
        
        k_1 <- 0.01*k1
        k_2 <- 0.01*k2
        k_3 <- 1
        
        f1 <- k1*((KD1*inputkd(times))/((KD1*inputkd(times))+1))*C
        f_1 <- k_1 * O
        f_3 <- k_3 *((0.5*A)/((0.5*A)+1))
        
        f2 <- k2*((KD2*inputkd(times))/((KD2*inputkd(times))+1))*O
        f_2 <- k_2 * A
        
        dC <- -f1+f_1+f_3
        dO <- f1 - f2 - f_1 + f_2
        
        
        return(list(c(dC, dO)))
      })
    }
    
    
    ss.init <- ode(y=c(C=1,O=0), times=seq(-3600,0,1), kdBaseModel0, parms=pars)
    ss.init <- as.data.frame(ss.init)
    
    
    out <- ode(y=c(C=ss.init[3601,2],O=ss.init[3601,3]), times=seq(0,10800,1), kdBaseModel, parms=pars)
    
    
    out <- as.data.frame(out)
    out$A <- 1-(out$C+out$O) 
    
    out
    
    return(out)
  }
  
  kdModel0 <- function(kdeg) { 
    
    
    signal2 <- data.frame(time = kdBase0(k1,k2,KD1,KD2)[,1],C = kdBase0(k1,k2,KD1,KD2)[,2],O = kdBase0(k1,k2,KD1,KD2)[,3])
    signal2$A <- 1-(signal2$C+signal2$O)
    sigimp2 <- approxfun(signal2$time, signal2$A, rule = 2)
    active <- sigimp2
    
    mrnafirst <- (1*active(0))/kdeg
    
    odeModel <- function (times, state, pars) {
      with(as.list(c(state, pars)),{
        
        dmrna <- (1*active(times))-(kdeg*mrna)
        
        return(list(c(dmrna)))
      })
    }
    
    state <- c(mrna=mrnafirst)
    
    out0 <- ode(y = state, times = seq(-3600,0,by=1),
                func = odeModel, parms = pars)
    out0 <- as.data.frame(out0)
    return(out0)
  }
  
  kdModel <- function(kdeg,tau) { 
    
    
    signal2 <- data.frame(time = kdBase0(k1,k2,KD1,KD2)[,1],C = kdBase0(k1,k2,KD1,KD2)[,2],O = kdBase0(k1,k2,KD1,KD2)[,3])
    signal2$A <- 1-(signal2$C+signal2$O)
    sigimp2 <- approxfun(signal2$time, signal2$A, rule = 2)
    active <- sigimp2
    
    mrnafirst <- kdModel0(kdeg)[3601,2]
    
    odeModel <- function (times, state, pars) {
      with(as.list(c(state, pars)),{
        
        dmrna <- (1*active(times-tau))-(kdeg*mrna)
        
        return(list(c(dmrna)))
      })
    }
    
    state <- c(mrna=mrnafirst)
    
    out0 <- ode(y = state, times = seq(0,10800,by=1),
                func = odeModel, parms = pars)
    out0 <- as.data.frame(out0)
    return(out0)
  }
  
  
  dftarget <- kdModel(kdeg,tau)
  
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

e <- function(k1,k2,KD1,KD2,kdeg,tau){  
  normwt <- data.frame(x=seq(0,10800,1),y=wtModels(k1,k2,KD1,KD2,kdeg,tau)[,2])
  minnorm <- normwt[1,2]
  for(i in 1:dim(normwt)[1]){
    normwt[i,2] <- normwt[i,2]/minnorm
  }
  
  normkd <- data.frame(x=seq(0,10800,1),y=kdModels(k1,k2,KD1,KD2,kdeg,tau)[,2])
  minnorm <- normkd[1,2]
  for(i in 1:dim(normkd)[1]){
    normkd[i,2] <- normkd[i,2]/minnorm
  }
  
  
  normdf <- rbind(normwt,normkd)
  normdf <- as.data.frame(normdf)
  
  normdf <- rescale(normdf[,2], to=c(0,1))
  norm <- c()
  norm <- data.frame(x=rep(seq(0,10800,1),2),y=normdf)
  
  return(norm)
}

g <- function(k1,k2,KD1,KD2,kdeg,tau){
  sim <- e(k1,k2,KD1,KD2,kdeg,tau)
  sim <- sim[which(sim[,1] %in% seq(0,10800,900)),]
  simwt <- sim[1:13,2]
  simkd <- sim[14:26,2]
  rmsdwt <- rmsd(refwtfc[,2],simwt,fit=FALSE)
  rmsdkd <- rmsd(refkdfc[,2],simkd,fit=FALSE)
  
  mul <- rmsdkd
  
  return(mul)
}

h <- function(p){
  g(p[1], p[2], p[3], p[4], p[5], p[6])
}

cm <- read.table("../First_step_optimization/Cycle_model/Bounds/Cycle_gene.name_random_initial.txt",sep="\t",skip=1)
cm <- cm[21:40,2:8]
cm <- as.data.frame(cm)


data <- read.table("../First_step_optimization/Cycle_model/Optimized_parameters/Cycle_gene.name_rep1.txt",sep="",skip=1)
data <- as.data.frame(data)
data <- data[which(data[,2]<0.5),]

k1min1 <- min(data[,3])
k1max1 <- max(data[,3])

k2min1 <- min(data[,4])
k2max1 <- max(data[,4])

KD1min1 <- min(data[,5])
KD1max1 <- max(data[,5])

KD2min1 <- min(data[,6])
KD2max1 <- max(data[,6])

kdegmin1 <- min(data[,7])
kdegmax1 <- max(data[,7])

taumin1 <- min(data[,8])
taumax1 <- max(data[,8])

data <- read.table("../First_step_optimization/Cycle_model/Optimized_parameters/Cycle_gene.name_rep2.txt",sep="",skip=1)
data <- as.data.frame(data)
data <- data[which(data[,2]<0.5),]

k1min2 <- min(data[,3])
k1max2 <- max(data[,3])

k2min2 <- min(data[,4])
k2max2 <- max(data[,4])

KD1min2 <- min(data[,5])
KD1max2 <- max(data[,5])

KD2min2 <- min(data[,6])
KD2max2 <- max(data[,6])

kdegmin2 <- min(data[,7])
kdegmax2 <- max(data[,7])

taumin2 <- min(data[,8])
taumax2 <- max(data[,8])

k1min <- min(k1min1,k1min2,k1max1,k1max2)
k1max <- max(k1min1,k1min2,k1max1,k1max2)

k2min <- min(k2min1,k2min2,k2max1,k2max2)
k2max <- max(k2min1,k2min2,k2max1,k2max2)

KD1min <- min(KD1min1,KD1min2,KD1max1,KD1max2)
KD1max <- max(KD1min1,KD1min2,KD1max1,KD1max2)

KD2min <- min(KD2min1,KD2min2,KD2max1,KD2max2)
KD2max <- max(KD2min1,KD2min2,KD2max1,KD2max2)

kdegmin <- min(kdegmin1,kdegmin2,kdegmax1,kdegmax2)
kdegmax <- max(kdegmin1,kdegmin2,kdegmax1,kdegmax2)

taumin <- min(taumin1,taumin2,taumax1,taumax2)
taumax <- max(taumin1,taumin2,taumax1,taumax2)

ps <- c()
ps0 <- c()
for(i in 1:dim(cm)[1]){
  pa1 <- cm[i,1]
  pa2 <- cm[i,2]
  pa3 <- cm[i,3]
  pa4 <- cm[i,4]
  pa5 <- cm[i,5]
  pa6 <- cm[i,6]
  
  pars <- c(k1=pa1,k2=pa2,KD1=pa3,KD2=pa4,kdeg=pa5,tau=pa6)
  
  result <- c()
  result <- sbplx(c(pa1,pa2,pa3,pa4,pa5,pa6),h,control=list(maxeval=200),lower=c(k1=k1min,k2=k2min,KD1=KD1min,KD2=KD2min,kdeg=kdegmin,tau=taumin),upper=c(k1=k1max,k2=k2max,KD1=KD1max,KD2=KD2max,kdeg=kdegmax,tau=taumax))
  
  kdegpar <- c()
  KD1par <- c()
  KD2par <- c()
  taupar <- c()
  rmsd <- c()
  k1par <- result$par[1]
  k2par <- result$par[2]
  KD1par <- result$par[3]
  KD2par <- result$par[4]
  kdegpar <- result$par[5]
  taupar <- result$par[6]
  rmsd <- result$value
  
  ps0 <- c(rmsd,k1par,k2par,KD1par,KD2par,kdegpar,taupar)
  
  ps <- c(ps,c(ps0))
}

ps <- data.frame(matrix(unlist(ps),ncol=7,byrow=T))

write.table(ps, "2_Cycle_gene.name_rep1_2.txt", sep = "\t")

