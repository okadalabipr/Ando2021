library(edgeR)
library(pracma)
library(deSolve)
library(DEoptim)
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

X2wt <- function(kdeg,kd) {
  parameters2 <- c(kdeg=kdeg,kd=kd)
  Lorenz2 <- function(times,state, parameters) {
    with(as.list(c(state, parameters)), {
      f <- ((kd*inputwt(times))/(1+kd*inputwt(times)))
      dmrna <- f-(kdeg*mrna)
      
      return(list(c(dmrna)))
    })
  }
  
  f <- ((kd*inputwt(0))/(1+kd*inputwt(0)))
  mrna0  <- f/kdeg
  
  state <- c(mrna=mrna0)
  
  out0 <- ode(y = state, times = seq(-3600,0,by=1),
              func = Lorenz2, parms = parameters2)
  out0 <- as.data.frame(out0)
  return(out0)
}

X1wt <- function(kdeg,kd,tau) {
  parameters <- c(kdeg=kdeg,kd=kd,tau=tau)
  
  f <- ((kd*inputwt(0))/(1+kd*inputwt(0)))
  mrna0  <- f/kdeg
  
  X2wt <- function(kdeg,kd) {
    parameters2 <- c(kdeg=kdeg,kd=kd)
    Lorenz2 <- function(times,state, parameters) {
      with(as.list(c(state, parameters)), {
        f <- ((kd*inputwt(times))/(1+kd*inputwt(times)))
        dmrna <- f-(kdeg*mrna)
        
        return(list(c(dmrna)))
      })
    }
    
    f <- ((kd*inputwt(0))/(1+kd*inputwt(0)))
    mrna0  <- f/kdeg
    
    state <- c(mrna=mrna0)
    
    out0 <- ode(y = state, times = seq(-3600,0,by=1),
                func = Lorenz2, parms = parameters2)
    out0 <- as.data.frame(out0)
    return(out0)
  }
  
  mrnafirst <- X2wt(kdeg,kd)[3601,2]
  
  Lorenz <- function(times,state, parameters) {
    with(as.list(c(state, parameters)), {
      f <- ((kd*inputwt(times-tau))/(1+kd*inputwt(times-tau)))
      dmrna <- f-kdeg*mrna
      
      return(list(c(dmrna)))
    })
  }
  
  state <- c(mrna=mrnafirst)
  
  out <- ode(y = state, times = seq(0,10800,by=1),
             func = Lorenz, parms = parameters)
  out <- as.data.frame(out)
  
  df <- rbind(X2wt(kdeg,kd),out[2:dim(out)[1],])
  x <- df[which(df[,1] %in% seq(0,10800,by=1)),1]
  y <- df[which(df[,1] %in% seq(0,10800,by=1)),2]
  df <- data.frame(time=x,mrna=y)
  
  return(df)
}

X2kd <- function(kdeg,kd) {
  parameters2 <- c(kdeg=kdeg,kd=kd)
  Lorenz2 <- function(times,state, parameters) {
    with(as.list(c(state, parameters)), {
      f <- ((kd*inputkd(times))/(1+kd*inputkd(times)))
      dmrna <- f-kdeg*mrna
      
      return(list(c(dmrna)))
    })
  }
  
  f <- ((kd*inputkd(0))/(1+kd*inputkd(0)))
  mrna0  <- f/kdeg
  
  state <- c(mrna=mrna0)
  
  out0 <- ode(y = state, times = seq(-3600,0,by=1),
              func = Lorenz2, parms = parameters2)
  out0 <- as.data.frame(out0)
  return(out0)
}

X1kd <- function(kdeg,kd,tau) {
  parameters <- c(kdeg=kdeg,kd=kd,tau=tau)
  
  f <- ((kd*inputkd(0))/(1+kd*inputkd(0)))
  mrna0  <- f/kdeg
  
  X2kd <- function(kdeg,kd) {
    parameters2 <- c(kdeg=kdeg,kd=kd)
    Lorenz2 <- function(times,state, parameters) {
      with(as.list(c(state, parameters)), {
        f <- ((kd*inputkd(times))/(1+kd*inputkd(times)))
        dmrna <- f-kdeg*mrna
        
        return(list(c(dmrna)))
      })
    }
    
    f <- ((kd*inputkd(0))/(1+kd*inputkd(0)))
    mrna0  <- f/kdeg
    
    state <- c(mrna=mrna0)
    
    out0 <- ode(y = state, times = seq(-3600,0,by=1),
                func = Lorenz2, parms = parameters2)
    out0 <- as.data.frame(out0)
    return(out0)
  }
  
  mrnafirst <- X2kd(kdeg,kd)[3601,2]
  
  Lorenz <- function(times,state, parameters) {
    with(as.list(c(state, parameters)), {
      f <- ((kd*inputkd(times-tau))/(1+kd*inputkd(times-tau)))
      dmrna <- f-kdeg*mrna
      
      return(list(c(dmrna)))
    })
  }
  
  state <- c(mrna=mrnafirst)
  
  out <- ode(y = state, times = seq(0,10800,by=1),
             func = Lorenz, parms = parameters)
  out <- as.data.frame(out)
  
  
  df <- rbind(X2kd(kdeg,kd),out[2:dim(out)[1],])
  x <- df[which(df[,1] %in% seq(0,10800,by=1)),1]
  y <- df[which(df[,1] %in% seq(0,10800,by=1)),2]
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




ewtfc <- function(kdeg,kd,tau){ 
  normwt <- data.frame(x=seq(0,10800,1),y=X1wt(kdeg,kd,tau)[,2])
  minnorm <- normwt[1,2]
  for(i in 1:dim(normwt)[1]){
    normwt[i,2] <- normwt[i,2]/minnorm
  }
  
  normkd <- data.frame(x=seq(0,10800,1),y=X1kd(kdeg,kd,tau)[,2])
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

ekdfc <- function(kdeg,kd,tau){  
  normwt <- data.frame(x=seq(0,10800,1),y=X1wt(kdeg,kd,tau)[,2])
  minnorm <- normwt[1,2]
  for(i in 1:dim(normwt)[1]){
    normwt[i,2] <- normwt[i,2]/minnorm
  }
  
  normkd <- data.frame(x=seq(0,10800,1),y=X1kd(kdeg,kd,tau)[,2])
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

g <- function(kdeg,kd,tau){
  
  rmsdwt <- rmsd(refwtfc[,2],ewtfc(kdeg,kd,tau)[which(ewtfc(kdeg,kd,tau)[,1] %in% seq(0,10800,900)),2],fit=FALSE)
  #rmsdkd <- rmsd(refkdfc[,2],ekdfc(kdeg,kd,tau)[which(ekdfc(kdeg,kd,tau)[,1] %in% seq(0,10800,900)),2],fit=FALSE)
  
  mul <- rmsdwt
  
  return(mul)
}



h <- function(p){
  g(p[1], p[2], p[3])
}

cm <- read.table("../Data/Subplex_random_initial.txt",sep="\t",skip=1)
cm <- cm[1:100,2:4]
cm <- as.data.frame(cm)

ps <- c()
ps0 <- c()
for(i in 1:dim(cm)[1]){
  pa1 <- cm[i,1]
  pa2 <- cm[i,2]
  pa3 <- cm[i,3]
  
  
  result <- c()
  result <- sbplx(c(pa1,pa2,pa3),h,control=list(maxeval=200),lower=c(kdeg=0.00002,kd=0.001,tau=0),upper=c(kdeg=0.002,kd=1000,tau=7200))
  
  kdegpar <- c()
  kdpar <- c()
  taupar <- c()
  rmsd <- c()
  kdegpar <- result$par[1]
  kdpar <- result$par[2]
  taupar <- result$par[3]
  rmsd <- result$value
  ps0 <- c(rmsd,kdegpar,kdpar,taupar)
  
  ps <- c(ps,c(ps0))
}
ps <- data.frame(matrix(unlist(ps),ncol=4,byrow=T))
write.table(ps, "gene.name_subplex_rep2.txt", sep = "\t")