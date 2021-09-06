library(edgeR)
library(pracma)
library(deSolve)
library(DEoptim)
library(scales)
library(bio3d)
library(car)
library(openxlsx)
library(ggplot2)


#Draw a heatmap of simulation from Simple model
smallest <- read.table("../Data/Simple_ERGsmallest.txt", sep = "\t",skip=1)
except_simple_ERG <- c("EFNA1","GADD45A","IRF1","IRS2","JUNB","KLF10","NUAK2","PHLDA1","PPP1R15A","PTGER4","RND1","SOX9","SPRY4","SPSB1","TNFAIP3","ZC3H12A")
smallest <- smallest[-which(smallest[,1] %in% except_simple_ERG),]
smallest <- smallest[-which(smallest[,1]=="NFKBIA"),]
Simple_ERG <- smallest[,1]

smallest <- read.table("../Data/Simple_IRGsmallest.txt", sep = "\t",skip=1)
Simple_IRG <- c("ADGRF4","AKR1C2","ANKRD18B","AREG","BAG1","BCL3","DMXL2","ETS2","FAM107B","FNBP4","ICAM1","IL17C","ITGAV","ITPKC","LYPD6B","NFKBIE","NKX3_1","NMD3","PCNA","PLAU","PLPP6","RAB9A","REL","S100A9","SDC4","SDCBP","SPTSSB","TET2","TXNRD1","VMP1")
smallest <- smallest[which(smallest[,1] %in% Simple_IRG),]
Simple_IRG <- smallest[,1]

smallest <- read.table("../Data/Simple_DRGsmallest.txt", sep = "\t",skip=1)
Simple_DRG <- c("APOL2","ARID5B","CHST15","CREBZF","CRIM1","DAXX","DRAM1","FAM117A","FGD6","HMGCR","IER5L","ITGB8","KYNU","LINC00052","LTB","OPTN","PHF23","PLEKHG3","PPP1R18","RELB","RFX5","SEMA4B","SGPP2","SLC6A14","STAT5A","TRAF3")
smallest <- smallest[which(smallest[,1] %in% Simple_DRG),]
Simple_DRG <- smallest[,1]

names <- smallest[,1]
smallest <- smallest[,2:3]

wtrmsdvalue <- c()
kdrmsdvalue <- c()
refall <- c()
for(r in 1:length(names)[1]){
  d <- names[r]
  
  #stackkari <- c()
  #stack <- c()
  #dfde <- c()
  #file <- paste("../Optimized_parameters/Simple_model/ERGsubcluster2/rep1/2_",d,"_subplex_rep1",".txt",sep="")
  #datafile <- read.table(file,sep="\t",skip=1) 
  #dfde <- datafile[,2:5]
  #dfde <- as.data.frame(dfde)
  #colnames(dfde) <- c("RMSD","kdeg","KD","tau")
  #
  #kdeg <- dfde[smallest[r,1],2]
  #kd <- dfde[smallest[r,1],3]
  #tau <- dfde[smallest[r,1],4]
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
  #normwt <- data.frame(x=seq(0,10800,1),y=X1wt(kdeg,kd,tau)[,2])
  #minnorm <- normwt[1,2]
  #for(i in 1:dim(normwt)[1]){
  #  normwt[i,2] <- normwt[i,2]/minnorm
  #}
  #
  #normkd <- data.frame(x=seq(0,10800,1),y=X1kd(kdeg,kd,tau)[,2])
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
  file <- paste("../Optimized_parameters/Simple_model/ERGsubcluster2/rep2/2_",d,"_subplex_rep2",".txt",sep="")
  datafile <- read.table(file,sep="\t",skip=1) 
  dfde <- datafile[,2:5]
  dfde <- as.data.frame(dfde)
  colnames(dfde) <- c("RMSD","kdeg","KD","tau")
  
  kdeg <- dfde[smallest[r,2],2]
  kd <- dfde[smallest[r,2],3]
  tau <- dfde[smallest[r,2],4]
  rmsd <- dfde[smallest[r,2],1]
  
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

Simple_refall2 <- data.frame(matrix(unlist(refall),ncol=26,byrow=T))

Simple_refall3 <- Simple_refall2[,1:13] 
Simple_refall3 <- Simple_refall2[,14:26]

#Draw a heatmap of simulation from IFFL model
smallest <- read.table("../Data/IFFL_ERGsmallest.txt", sep = "\t",skip=1)
IFFL_ERG <- c("EFNA1","IRF1","SPRY4")
smallest <- smallest[which(smallest[,1] %in% IFFL_ERG),]

smallest <- read.table("../Data/IFFL_IRGsmallest.txt", sep = "\t",skip=1)
IFFL_IRG <- c("ANKRD33B")
smallest <- smallest[which(smallest[,1] %in% IFFL_IRG),]

smallest <- read.table("../Data/IFFL_DRGsmallest.txt", sep = "\t",skip=1)
IFFL_DRG <- c("EFEMP1","IL6ST","RAB3IP","RBM25","SEMA3C","TRIP10")
smallest <- smallest[which(smallest[,1] %in% IFFL_DRG),]

names <- smallest[,1]
smallest <- smallest[,2:3]

wtrmsdvalue <- c()
kdrmsdvalue <- c()
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

IFFL_refall2 <- data.frame(matrix(unlist(refall),ncol=26,byrow=T))

IFFL_refall3 <- IFFL_refall2[,1:13] 
IFFL_refall3 <- IFFL_refall2[,14:26]

#Draw a heatmap of simulation from Cycle model
smallest <- read.table("../Data/Cycle_ERGsmallest.txt", sep = "\t",skip=1)
Cycle_ERG <- c("ZC3H12A")
smallest <- smallest[which(smallest[,1] %in% Cycle_ERG),]

smallest <- read.table("../Data/Cycle_IRGsmallest.txt", sep = "\t",skip=1)
Cycle_IRG <- c("ADGRG6","AMIGO2","ATP1B1","BCL2L11","BID","BIRC3","CLIC4","RBM3","RHOV","TNFRSF11B")
smallest <- smallest[which(smallest[,1] %in% Cycle_IRG),]

smallest <- read.table("../Data/Cycle_DRGsmallest.txt", sep = "\t",skip=1)
Cycle_DRG <- c("ABAT","ABTB2","AOX1","FOXO1","IFNGR2","IRF2","LINC02015","MGAT4A","NFATC2","NFKB1","PRKAR2B","RUNX2","TAP1","TBC1D9","TNFAIP2","VDR")
smallest <- smallest[which(smallest[,1] %in% Cycle_DRG),]

names <- smallest[,1]
smallest <- smallest[,2:3]

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

wtrmsdvalue <- c()
kdrmsdvalue <- c()
refall <- c()
for(r in 1:length(names)[1]){
  d <- names[r]
  
  #stackkari <- c()
  #stack <- c()
  #dfde <- c()
  #file <- paste("../Optimized_parameters/Cycle_model/ERGsubcluster2/rep1/2_Cycle_",d,"_rep1",".txt",sep="")
  #datafile <- read.table(file,sep="\t",skip=1) 
  #datafile <- as.data.frame(datafile)
  #dfde <- datafile[,2:8]
  #dfde <- as.data.frame(dfde)
  #colnames(dfde) <- c("RMSD","k1","k2","KD1","KD2","kdeg","tau")
  #rep1dfde <- dfde
  #
  #k1 <- dfde[smallest[r,1],2]
  #k2 <- dfde[smallest[r,1],3]
  #KD1 <- dfde[smallest[r,1],4]
  #KD2 <- dfde[smallest[r,1],5]
  #kdeg <- dfde[smallest[r,1],6]
  #tau <- dfde[smallest[r,1],7]
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
  #normwt <- data.frame(x=seq(0,10800,1),y=wtModels(k1,k2,KD1,KD2,kdeg,tau)[,2])
  #minnorm <- normwt[1,2]
  #for(i in 1:dim(normwt)[1]){
  #  normwt[i,2] <- normwt[i,2]/minnorm
  #}
  #
  #normkd <- data.frame(x=seq(0,10800,1),y=kdModels(k1,k2,KD1,KD2,kdeg,tau)[,2])
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
  file <- paste("../Optimized_parameters/Cycle_model/ERGsubcluster2/rep2/2_Cycle_",d,"_rep2",".txt",sep="")
  datafile <- read.table(file,sep="\t",skip=1) 
  datafile <- as.data.frame(datafile)
  dfde <- datafile[,2:8]
  dfde <- as.data.frame(dfde)
  colnames(dfde) <- c("RMSD","k1","k2","KD1","KD2","kdeg","tau")
  rep2dfde <- dfde
  
  k1 <- dfde[smallest[r,2],2]
  k2 <- dfde[smallest[r,2],3]
  KD1 <- dfde[smallest[r,2],4]
  KD2 <- dfde[smallest[r,2],5]
  kdeg <- dfde[smallest[r,2],6]
  tau <- dfde[smallest[r,2],7]
  rmsd <- dfde[smallest[r,2],1]
  
  count <- read.csv("../Data/160617_ikbakd.csv")
  count1=count[,2:12]#6`11ñÚðo
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

Cycle_refall2 <- data.frame(matrix(unlist(refall),ncol=26,byrow=T))

Cycle_refall3 <- Cycle_refall2[,1:13] 
Cycle_refall3 <- Cycle_refall2[,14:26]

#Draw a heatmap of simulation from Model v4
smallest <- read.table("../Data/Model_v4_ERGsmallest.txt", sep = "\t",skip=1)
Combination_ERG <- c("GADD45A","IRS2","JUNB","KLF10","NUAK2","PHLDA1","PPP1R15A","PTGER4","RND1","SOX9","SPSB1","TNFAIP3")
smallest <- smallest[which(smallest[,1] %in% Combination_ERG),]

smallest <- read.table("../Data/Model_v4_IRGsmallest.txt", sep = "\t",skip=1)
Combination_IRG <- c("ANXA8","CDKN2B","EIF5","MAFF","NCOA7","NEDD9","TICAM1")
smallest <- smallest[which(smallest[,1] %in% Combination_IRG),]

smallest <- read.table("../Data/Model_v4_DRGsmallest.txt", sep = "\t",skip=1)
Combination_DRG <- c("B4GALT5","BAZ1A","CLK4","DDX58","HIVEP2","KLHL5","LAMC2","NFKB2","PHLDB2","RSRC2","SERPINB8","TUBD1")
smallest <- smallest[which(smallest[,1] %in% Combination_DRG),]

names <- smallest[,1]
smallest <- smallest[,2:3]

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

wtBase0 <- function(k1,k2,KD1,KD2,tau){ 
  
  pars <- c(k1=k1,k2=k2,KD1=KD1,KD2=KD2,tau=tau)
  
  
  
  wtBaseModel0 <- function (times, State, pars) {
    
    pars <- c(k1=k1,k2=k2,KD1=KD1,KD2=KD2,tau=tau)
    
    
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
    
    pars <- c(k1=k1,k2=k2,KD1=KD1,KD2=KD2,tau=tau)
    
    with(as.list(c(State, pars)),{
      
      A = 1-C-O
      k_1 <- 0.01*k1
      k_2 <- 0.01*k2
      k_3 <- 1
      
      f1 <- k1*((KD1*inputwt(times-tau))/((KD1*inputwt(times-tau))+1))*C
      f_1 <- k_1 * O
      f_3 <- k_3 *((0.5*A)/((0.5*A)+1))
      
      f2 <- k2*((KD2*inputwt(times-tau))/((KD2*inputwt(times-tau))+1))*O
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

wtModel0 <- function(kdeg,KDTF2) { 
  
  
  signal2 <- data.frame(time = wtBase0(k1,k2,KD1,KD2,tau)[,1],C = wtBase0(k1,k2,KD1,KD2,tau)[,2],O = wtBase0(k1,k2,KD1,KD2,tau)[,3])
  signal2$A <- 1-(signal2$C+signal2$O)
  sigimp2 <- approxfun(signal2$time, signal2$A, rule = 2)
  active <- sigimp2
  
  signal2 <- data.frame(time = D2wt(kdegTF,KDTF)[,1],tf = D2wt(kdegTF,KDTF)[,2])
  sigimp2 <- approxfun(signal2$time, signal2$tf, rule = 2)
  competitor0wt <- sigimp2
  
  
  mrnafirst <- (1*((active(0))^2)/(((active(0))^2)+((KDTF2*competitor0wt(0))^2)+1))/kdeg
  
  odeModel <- function (times, state, pars) {
    with(as.list(c(state, pars)),{
      
      dmrna <- (1*((active(times))^2)/(((active(times))^2)+((KDTF2*competitor0wt(times))^2)+1))-(kdeg*mrna)
      
      return(list(c(dmrna)))
    })
  }
  
  state <- c(mrna=mrnafirst)
  
  out0 <- ode(y = state, times = seq(-3600,0,by=1),
              func = odeModel, parms = pars)
  out0 <- as.data.frame(out0)
  return(out0)
}

wtModel <- function(kdeg,KDTF2,tau) { 
  
  
  signal2 <- data.frame(time = wtBase0(k1,k2,KD1,KD2,tau)[,1],C = wtBase0(k1,k2,KD1,KD2,tau)[,2],O = wtBase0(k1,k2,KD1,KD2,tau)[,3])
  signal2$A <- 1-(signal2$C+signal2$O)
  sigimp2 <- approxfun(signal2$time, signal2$A, rule = 2)
  active <- sigimp2
  
  signal3 <- data.frame(time = D4wt(kdegTF,KDTF)[,1],tf = D4wt(kdegTF,KDTF)[,2])
  sigimp3 <- approxfun(signal3$time, signal3$tf, rule = 2)
  competitor1wt <- sigimp3
  
  mrnafirst <- wtModel0(kdeg,KDTF2)[3601,2]
  
  odeModel <- function (times, state, pars) {
    with(as.list(c(state, pars)),{
      
      dmrna <- (1*((active(times-tau))^2)/(((active(times-tau))^2)+((KDTF2*competitor1wt(times-7200))^2)+1))-(kdeg*mrna)
      
      return(list(c(dmrna)))
    })
  }
  
  state <- c(mrna=mrnafirst)
  
  out0 <- ode(y = state, times = seq(0,10800,by=1),
              func = odeModel, parms = pars)
  out0 <- as.data.frame(out0)
  return(out0)
}

wtModels <- function(k1,k2,KD1,KD2,KDTF2,kdeg,tau) { 
  
  pars <- c(k1=k1,k2=k2,KD1=KD1,KD2=KD2,KDTF2=KDTF2,kdeg=kdeg,tau=tau)
  
  wtBase0 <- function(k1,k2,KD1,KD2,tau) { 
    
    pars <- c(k1=k1,k2=k2,KD1=KD1,KD2=KD2,tau=tau)
    
    
    wtBaseModel0 <- function (times, State, pars) {
      
      pars <- c(k1=k1,k2=k2,KD1=KD1,KD2=KD2,tau=tau)
      
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
      
      pars <- c(k1=k1,k2=k2,KD1=KD1,KD2=KD2,tau=tau)
      
      with(as.list(c(State, pars)),{
        
        A = 1-C-O
        k_1 <- 0.01*k1
        k_2 <- 0.01*k2
        k_3 <- 1
        
        f1 <- k1*((KD1*inputwt(times-tau))/((KD1*inputwt(times-tau))+1))*C
        f_1 <- k_1 * O
        f_3 <- k_3 *((0.5*A)/((0.5*A)+1))
        
        f2 <- k2*((KD2*inputwt(times-tau))/((KD2*inputwt(times-tau))+1))*O
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
  
  wtModel0 <- function(kdeg,KDTF2) { 
    
    
    signal2 <- data.frame(time = wtBase0(k1,k2,KD1,KD2,tau)[,1],C = wtBase0(k1,k2,KD1,KD2,tau)[,2],O = wtBase0(k1,k2,KD1,KD2,tau)[,3])
    signal2$A <- 1-(signal2$C+signal2$O)
    sigimp2 <- approxfun(signal2$time, signal2$A, rule = 2)
    active <- sigimp2
    
    signal2 <- data.frame(time = D2wt(kdegTF,KDTF)[,1],tf = D2wt(kdegTF,KDTF)[,2])
    sigimp2 <- approxfun(signal2$time, signal2$tf, rule = 2)
    competitor0wt <- sigimp2
    
    mrnafirst <- (1*((active(0))^2)/(((active(0))^2)+((KDTF2*competitor0wt(0))^2)+1))/kdeg
    
    odeModel <- function (times, state, pars) {
      with(as.list(c(state, pars)),{
        
        dmrna <- (1*((active(times))^2)/(((active(times))^2)+((KDTF2*competitor0wt(times))^2)+1))-(kdeg*mrna)
        
        return(list(c(dmrna)))
      })
    }
    
    state <- c(mrna=mrnafirst)
    
    out0 <- ode(y = state, times = seq(-3600,0,by=1),
                func = odeModel, parms = pars)
    out0 <- as.data.frame(out0)
    return(out0)
  }
  
  wtModel <- function(kdeg,KDTF2,tau) { 
    
    
    
    signal2 <- data.frame(time = wtBase0(k1,k2,KD1,KD2,tau)[,1],C = wtBase0(k1,k2,KD1,KD2,tau)[,2],O = wtBase0(k1,k2,KD1,KD2,tau)[,3])
    signal2$A <- 1-(signal2$C+signal2$O)
    sigimp2 <- approxfun(signal2$time, signal2$A, rule = 2)
    active <- sigimp2
    
    signal3 <- data.frame(time = D4wt(kdegTF,KDTF)[,1],tf = D4wt(kdegTF,KDTF)[,2])
    sigimp3 <- approxfun(signal3$time, signal3$tf, rule = 2)
    competitor1wt <- sigimp3
    
    mrnafirst <- wtModel0(kdeg,KDTF2)[3601,2]
    
    odeModel <- function (times, state, pars) {
      with(as.list(c(state, pars)),{
        
        dmrna <- (1*((active(times-tau))^2)/(((active(times-tau))^2)+((KDTF2*competitor1wt(times-7200))^2)+1))-(kdeg*mrna)
        
        return(list(c(dmrna)))
      })
    }
    
    state <- c(mrna=mrnafirst)
    
    out0 <- ode(y = state, times = seq(0,10800,by=1),
                func = odeModel, parms = pars)
    out0 <- as.data.frame(out0)
    return(out0)
  }
  
  
  
  dftarget <- wtModel(kdeg,KDTF2,tau)
  
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

kdBase0 <- function(k1,k2,KD1,KD2,tau) { 
  
  pars <- c(k1=k1,k2=k2,KD1=KD1,KD2=KD2,tau=tau)
  
  kdBaseModel0 <- function (times, State, pars) {
    
    pars <- c(k1=k1,k2=k2,KD1=KD1,KD2=KD2,tau=tau)
    
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
    
    pars <- c(k1=k1,k2=k2,KD1=KD1,KD2=KD2,tau=tau)
    
    with(as.list(c(State, pars)),{
      
      A = 1-C-O
      
      k_1 <- 0.01*k1
      k_2 <- 0.01*k2
      k_3 <- 1
      
      f1 <- k1*((KD1*inputkd(times-tau))/((KD1*inputkd(times-tau))+1))*C
      f_1 <- k_1 * O
      f_3 <- k_3 *((0.5*A)/((0.5*A)+1))
      
      f2 <- k2*((KD2*inputkd(times-tau))/((KD2*inputkd(times-tau))+1))*O
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

kdModel0 <- function(kdeg,KDTF2) { 
  
  
  signal2 <- data.frame(time = kdBase0(k1,k2,KD1,KD2,tau)[,1],C = kdBase0(k1,k2,KD1,KD2,tau)[,2],O = kdBase0(k1,k2,KD1,KD2,tau)[,3])
  signal2$A <- 1-(signal2$C+signal2$O)
  sigimp2 <- approxfun(signal2$time, signal2$A, rule = 2)
  active <- sigimp2
  
  signal2 <- data.frame(time = D2kd(kdegTF,KDTF)[,1],tf = D2kd(kdegTF,KDTF)[,2])
  sigimp2 <- approxfun(signal2$time, signal2$tf, rule = 2)
  competitor0kd <- sigimp2
  
  mrnafirst <- (1*((active(0))^2)/(((active(0))^2)+((KDTF2*competitor0kd(0))^2)+1))/kdeg
  
  odeModel <- function (times, state, pars) {
    with(as.list(c(state, pars)),{
      
      dmrna <- (1*((active(times))^2)/(((active(times))^2)+((KDTF2*competitor0kd(times))^2)+1))-(kdeg*mrna)
      
      return(list(c(dmrna)))
    })
  }
  
  state <- c(mrna=mrnafirst)
  
  out0 <- ode(y = state, times = seq(-3600,0,by=1),
              func = odeModel, parms = pars)
  out0 <- as.data.frame(out0)
  return(out0)
}

kdModel <- function(kdeg,KDTF2,tau) { 
  
  signal2 <- data.frame(time = kdBase0(k1,k2,KD1,KD2,tau)[,1],C = kdBase0(k1,k2,KD1,KD2,tau)[,2],O = kdBase0(k1,k2,KD1,KD2,tau)[,3])
  signal2$A <- 1-(signal2$C+signal2$O)
  sigimp2 <- approxfun(signal2$time, signal2$A, rule = 2)
  active <- sigimp2
  
  signal2 <- data.frame(time = D4kd(kdegTF,KDTF)[,1],tf = D4kd(kdegTF,KDTF)[,2])
  sigimp2 <- approxfun(signal2$time, signal2$tf, rule = 2)
  competitor1kd <- sigimp2
  
  mrnafirst <- kdModel0(kdeg,KDTF2)[3601,2]
  
  odeModel <- function (times, state, pars) {
    with(as.list(c(state, pars)),{
      
      dmrna <- (1*((active(times-tau))^2)/(((active(times-tau))^2)+((KDTF2*competitor1kd(times-7200))^2)+1))-(kdeg*mrna)
      
      return(list(c(dmrna)))
    })
  }
  
  state <- c(mrna=mrnafirst)
  
  out0 <- ode(y = state, times = seq(0,10800,by=1),
              func = odeModel, parms = pars)
  out0 <- as.data.frame(out0)
  return(out0)
}

kdModels <- function(k1,k2,KD1,KD2,KDTF2,kdeg,tau) { 
  
  pars <- c(k1=k1,k2=k2,KD1=KD1,KD2=KD2,KDTF2=KDTF2,kdeg=kdeg,tau=tau)
  
  kdBase0 <- function(k1,k2,KD1,KD2,tau) { 
    
    pars <- c(k1=k1,k2=k2,KD1=KD1,KD2=KD2,tau=tau)
    
    kdBaseModel0 <- function (times, State, pars) {
      
      pars <- c(k1=k1,k2=k2,KD1=KD1,KD2=KD2,tau=tau)
      
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
      
      pars <- c(k1=k1,k2=k2,KD1=KD1,KD2=KD2,tau=tau)
      
      with(as.list(c(State, pars)),{
        
        A = 1-C-O
        
        k_1 <- 0.01*k1
        k_2 <- 0.01*k2
        k_3 <- 1
        
        f1 <- k1*((KD1*inputkd(times-tau))/((KD1*inputkd(times-tau))+1))*C
        f_1 <- k_1 * O
        f_3 <- k_3 *((0.5*A)/((0.5*A)+1))
        
        f2 <- k2*((KD2*inputkd(times-tau))/((KD2*inputkd(times-tau))+1))*O
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
  
  kdModel0 <- function(kdeg,KDTF2) { 
    
    
    signal2 <- data.frame(time = kdBase0(k1,k2,KD1,KD2,tau)[,1],C = kdBase0(k1,k2,KD1,KD2,tau)[,2],O = kdBase0(k1,k2,KD1,KD2,tau)[,3])
    signal2$A <- 1-(signal2$C+signal2$O)
    sigimp2 <- approxfun(signal2$time, signal2$A, rule = 2)
    active <- sigimp2
    
    signal2 <- data.frame(time = D2kd(kdegTF,KDTF)[,1],tf = D2kd(kdegTF,KDTF)[,2])
    sigimp2 <- approxfun(signal2$time, signal2$tf, rule = 2)
    competitor0kd <- sigimp2
    
    mrnafirst <- (1*((active(0))^2)/(((active(0))^2)+((KDTF2*competitor0kd(0))^2)+1))/kdeg
    
    odeModel <- function (times, state, pars) {
      with(as.list(c(state, pars)),{
        
        dmrna <- (1*((active(times))^2)/(((active(times))^2)+((KDTF2*competitor0kd(times))^2)+1))-(kdeg*mrna)
        
        return(list(c(dmrna)))
      })
    }
    
    state <- c(mrna=mrnafirst)
    
    out0 <- ode(y = state, times = seq(-3600,0,by=1),
                func = odeModel, parms = pars)
    out0 <- as.data.frame(out0)
    return(out0)
  }
  
  kdModel <- function(kdeg,KDTF2,tau) { 
    
    
    signal2 <- data.frame(time = kdBase0(k1,k2,KD1,KD2,tau)[,1],C = kdBase0(k1,k2,KD1,KD2,tau)[,2],O = kdBase0(k1,k2,KD1,KD2,tau)[,3])
    signal2$A <- 1-(signal2$C+signal2$O)
    sigimp2 <- approxfun(signal2$time, signal2$A, rule = 2)
    active <- sigimp2
    
    signal2 <- data.frame(time = D4kd(kdegTF,KDTF)[,1],tf = D4kd(kdegTF,KDTF)[,2])
    sigimp2 <- approxfun(signal2$time, signal2$tf, rule = 2)
    competitor1kd <- sigimp2
    
    mrnafirst <- kdModel0(kdeg,KDTF2)[3601,2]
    
    odeModel <- function (times, state, pars) {
      with(as.list(c(state, pars)),{
        
        dmrna <- (1*((active(times-tau))^2)/(((active(times-tau))^2)+((KDTF2*competitor1kd(times-7200))^2)+1))-(kdeg*mrna)
        
        
        return(list(c(dmrna)))
      })
    }
    
    state <- c(mrna=mrnafirst)
    
    out0 <- ode(y = state, times = seq(0,10800,by=1),
                func = odeModel, parms = pars)
    out0 <- as.data.frame(out0)
    return(out0)
  }
  
  
  dftarget <- kdModel(kdeg,KDTF2,tau)
  
  x <- dftarget[which(dftarget[,1] %in% seq(0,10800,by=1)),1]
  y <- dftarget[which(dftarget[,1] %in% seq(0,10800,by=1)),2]
  df <- data.frame(time=x,mrna=y)
  
  return(df)
}

wtrmsdvalue <- c()
kdrmsdvalue <- c()
refall <- c()
for(r in 1:length(names)[1]){
  d <- names[r]
  
  #stackkari <- c()
  #stack <- c()
  #dfde <- c()
  #file <- paste("../Optimized_parameters/Model_v4/ERGsubcluster2/rep1/2_Cycle_",d,"_rep1",".txt",sep="")
  #datafile <- read.table(file,sep="\t",skip=1) 
  #datafile <- as.data.frame(datafile)
  #dfde <- datafile[,2:9]
  #dfde <- as.data.frame(dfde)
  #colnames(dfde) <- c("RMSD","k1","k2","KD1","KD2","KDTF2","kdeg","tau")
  #rep1dfde <- dfde
  #
  #k1 <- dfde[smallest[r,1],2]
  #k2 <- dfde[smallest[r,1],3]
  #KD1 <- dfde[smallest[r,1],4]
  #KD2 <- dfde[smallest[r,1],5]
  #KDTF2 <- dfde[smallest[r,1],6]
  #kdeg <- dfde[smallest[r,1],7]
  #tau <- dfde[smallest[r,1],8]
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
  #normwt <- data.frame(x=seq(0,10800,1),y=wtModels(k1,k2,KD1,KD2,KDTF2,kdeg,tau)[,2])
  #minnorm <- normwt[1,2]
  #for(i in 1:dim(normwt)[1]){
  #  normwt[i,2] <- normwt[i,2]/minnorm
  #}
  #
  #normkd <- data.frame(x=seq(0,10800,1),y=kdModels(k1,k2,KD1,KD2,KDTF2,kdeg,tau)[,2])
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
  file <- paste("../Optimized_parameters/Model_v4/ERGsubcluster2/rep2/2_Cycle_",d,"_rep2",".txt",sep="")
  datafile <- read.table(file,sep="\t",skip=1) 
  datafile <- as.data.frame(datafile)
  dfde <- datafile[,2:9]
  dfde <- as.data.frame(dfde)
  colnames(dfde) <- c("RMSD","k1","k2","KD1","KD2","KDTF2","kdeg","tau")
  rep2dfde <- dfde
  
  k1 <- dfde[smallest[r,2],2]
  k2 <- dfde[smallest[r,2],3]
  KD1 <- dfde[smallest[r,2],4]
  KD2 <- dfde[smallest[r,2],5]
  KDTF2 <- dfde[smallest[r,2],6]
  kdeg <- dfde[smallest[r,2],7]
  tau <- dfde[smallest[r,2],8]
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
  
  normwt <- data.frame(x=seq(0,10800,1),y=wtModels(k1,k2,KD1,KD2,KDTF2,kdeg,tau)[,2])
  minnorm <- normwt[1,2]
  for(i in 1:dim(normwt)[1]){
    normwt[i,2] <- normwt[i,2]/minnorm
  }
  
  normkd <- data.frame(x=seq(0,10800,1),y=kdModels(k1,k2,KD1,KD2,KDTF2,kdeg,tau)[,2])
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

Combination_refall2 <- data.frame(matrix(unlist(refall),ncol=26,byrow=T))

Combination_refall3 <- Combination_refall2[,1:13] 
Cycle_refall3 <- Cycle_refall2[,1:13] 
IFFL_refall3 <- IFFL_refall2[,1:13] 
Simple_refall3 <- Simple_refall2[,1:13] 

Combination_refall3 <- Combination_refall2[,14:26]
Cycle_refall3 <- Cycle_refall2[,14:26] 
IFFL_refall3 <- IFFL_refall2[,14:26] 
Simple_refall3 <- Simple_refall2[,14:26]

refall3 <- rbind(Simple_refall3,IFFL_refall3,Cycle_refall3,Combination_refall3)
refall3 <- as.data.frame(refall3)
Simple_ERG <- as.character(Simple_ERG)
IFFL_ERG <- as.character(IFFL_ERG)
Cycle_ERG <- as.character(Cycle_ERG)
Combination_ERG <- as.character(Combination_ERG)
rownames(refall3) <- c(Simple_ERG,IFFL_ERG,Cycle_ERG,Combination_ERG)
refall3 <- refall3[order(rownames(refall3)),]
refall3 <- as.data.frame(refall3)
names <- rownames(refall3)

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
  scale_fill_gradient(low="white",high="red",limits=c(0,1.0))+
  #scale_fill_gradient(low="white",high="mediumseagreen",limits=c(0,1.0))+
  labs(x = "",y = "") + 
  theme(panel.grid = element_blank())+
  scale_x_discrete(labels=c("0","15","30","45","60","75","90","105","120","135","150","165","180"),breaks=c("0","15","30","45","60","75","90","105","120","135","150","165","180"),expand = c(0, 0))+ 
  scale_y_discrete(breaks=sort(names,decreasing = T),labels=waiver(),limits=sort(names,decreasing = T),expand = c(0, 0))+
  theme(legend.position = "bottom",axis.ticks = element_blank(), axis.text.x = element_text(size = 7, angle = 0, hjust = 0.5, colour = "black"),axis.text.y = element_text(size = 6.0, angle = 0, hjust = 1.0, colour = "black"));g


#Draw a heatmap of nRMSD
Other_ERG <- c("ARRDC3","DLC1","ELF3","FAM46A","IER3","IER5","KLF3","LIF","MAP3K1","MXD1","NRIP1","PIM1","PLK2","PPP1R3C","SALL4","TNF")
Other_IRG <- c("ADGRF4","AKR1C2","AREG","BAG1","DMXL2","ETS2","FNBP4","ICAM1","IL17C","ITGAV","ITPKC","NKX3_1","NMD3","PCNA","PLAU","PLPP6","RAB9A","REL","SPTSSB","TET2","TXNRD1","VMP1")
Other_DRG <- c("HMGCR","IER5L","LINC00052","SLC6A14")

Simple_ERG <- c("MAP3K8","PLEKHF2")
Simple_IRG <- c("ANKRD18B","BCL3","FAM107B","LYPD6B","NFKBIE","S100A9","SDC4","SDCBP")
Simple_DRG <- c("APOL2","ARID5B","CHST15","CREBZF","CRIM1","DAXX","DRAM1","FAM117A","FGD6","ITGB8","KYNU","LTB","OPTN","PHF23","PLEKHG3","PPP1R18","RELB","RFX5","SEMA4B","SGPP2","STAT5A","TRAF3")

IFFL_ERG <- c("EFNA1","IRF1","SPRY4")
IFFL_IRG <- c("ANKRD33B")
IFFL_DRG <- c("EFEMP1","IL6ST","RAB3IP","RBM25","SEMA3C","TRIP10")

Cycle_ERG <- c("ZC3H12A")
Cycle_IRG <- c("ADGRG6","AMIGO2","ATP1B1","BCL2L11","BID","BIRC3","CLIC4","RBM3","RHOV","TNFRSF11B")
Cycle_DRG <- c("ABAT","ABTB2","AOX1","FOXO1","IFNGR2","IRF2","LINC02015","MGAT4A","NFATC2","NFKB1","PRKAR2B","RUNX2","TAP1","TBC1D9","TNFAIP2","VDR")

Combination_ERG <- c("GADD45A","IRS2","JUNB","KLF10","NUAK2","PHLDA1","PPP1R15A","PTGER4","RND1","SOX9","SPSB1","TNFAIP3")
Combination_IRG <- c("ANXA8","CDKN2B","EIF5","MAFF","NCOA7","NEDD9","TICAM1")
Combination_DRG <- c("B4GALT5","BAZ1A","CLK4","DDX58","HIVEP2","KLHL5","LAMC2","NFKB2","PHLDB2","RSRC2","SERPINB8","TUBD1")

Other_ERG <- data.frame(name=Other_ERG,group=rep("grey85",length(Other_ERG)))
Simple_ERG <- data.frame(name=Simple_ERG,group=rep("deepskyblue2",length(Simple_ERG)))
IFFL_ERG <- data.frame(name=IFFL_ERG,group=rep("red3",length(IFFL_ERG)))
Cycle_ERG <- data.frame(name=Cycle_ERG,group=rep("yellowgreen",length(Cycle_ERG)))
Combination_ERG <- data.frame(name=Combination_ERG,group=rep("mediumpurple1",length(Combination_ERG)))
ERG_df <- rbind(Other_ERG,Simple_ERG,IFFL_ERG,Cycle_ERG,Combination_ERG)
ERG_df <- as.data.frame(ERG_df)
ERG_df$name <- as.character(ERG_df$name)
ERG_df$group <- as.character(ERG_df$group)
ERG_df <- ERG_df[order(ERG_df[,1]),]
ERG_df$x <- rep("all",dim(ERG_df)[1])

Other_IRG <- data.frame(name=Other_IRG,group=rep("grey85",length(Other_IRG)))
Simple_IRG <- data.frame(name=Simple_IRG,group=rep("deepskyblue2",length(Simple_IRG)))
IFFL_IRG <- data.frame(name=IFFL_IRG,group=rep("red3",length(IFFL_IRG)))
Cycle_IRG <- data.frame(name=Cycle_IRG,group=rep("yellowgreen",length(Cycle_IRG)))
Combination_IRG <- data.frame(name=Combination_IRG,group=rep("mediumpurple1",length(Combination_IRG)))
IRG_df <- rbind(Other_IRG,Simple_IRG,IFFL_IRG,Cycle_IRG,Combination_IRG)
IRG_df <- as.data.frame(IRG_df)
IRG_df$name <- as.character(IRG_df$name)
IRG_df$group <- as.character(IRG_df$group)
IRG_df <- IRG_df[order(IRG_df[,1]),]
IRG_df$x <- rep("all",dim(IRG_df)[1])

Other_DRG <- data.frame(name=Other_DRG,group=rep("grey85",length(Other_DRG)))
Simple_DRG <- data.frame(name=Simple_DRG,group=rep("deepskyblue2",length(Simple_DRG)))
IFFL_DRG <- data.frame(name=IFFL_DRG,group=rep("red3",length(IFFL_DRG)))
Cycle_DRG <- data.frame(name=Cycle_DRG,group=rep("yellowgreen",length(Cycle_DRG)))
Combination_DRG <- data.frame(name=Combination_DRG,group=rep("mediumpurple1",length(Combination_DRG)))
DRG_df <- rbind(Other_DRG,Simple_DRG,IFFL_DRG,Cycle_DRG,Combination_DRG)
DRG_df <- as.data.frame(DRG_df)
DRG_df$name <- as.character(DRG_df$name)
DRG_df$group <- as.character(DRG_df$group)
DRG_df <- DRG_df[order(DRG_df[,1]),]
DRG_df$x <- rep("all",dim(DRG_df)[1])

df <- ERG_df

g <- ggplot(df,aes(as.factor(x),as.factor(name)))+
  geom_tile(aes(fill=group),color="grey50",stat="identity",na.rm = TRUE)+
  scale_fill_manual(values=c(grey85="grey85", deepskyblue2="deepskyblue2",red3="firebrick2",yellowgreen="yellowgreen",mediumpurple1="mediumpurple1"))+
  labs(x = "",y = "") + 
  theme(panel.grid = element_blank())+
  scale_x_discrete(labels=waiver(),breaks=group,expand = c(0, 0))+ 
  scale_y_discrete(breaks=sort(df$name,decreasing=T),labels=waiver(),limits=sort(df$name,decreasing=T),expand = c(0, 0))+
  theme(legend.position = "right",axis.ticks = element_blank(), axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, colour = "black"),axis.text.y = element_text(size = 12, angle = 0, hjust = 0.5, colour = "black"));g

#Draw a bar plot of the number of good fit genes from each model
simple <- c((2/18)*100)
IFFL <- c((3/18)*100)
cycle <- c((1/18)*100)
combination <- c((12/18)*100)
#other <- c(22)

df <- data.frame(simple,IFFL,cycle,combination)
df <- t(df)
df <- as.data.frame(df)
colnames(df) <- c("value")

g <- ggplot(df, aes(x = rownames(df), y = value, fill=rownames(df)))+
  geom_bar(stat="identity", color="black",width=0.4, size=0.7)+
  scale_fill_manual(values = c(simple="#66CCFF",IFFL="firebrick1",cycle="olivedrab2",combination="#9966FF"))+
  labs(x = "",y = "") + 
  theme_bw(base_size = 16)+
  theme(panel.grid = element_blank())+
  scale_x_discrete(breaks=c("simple","IFFL","cycle","combination"),labels=waiver(),limits=c("simple","IFFL","cycle","combination"),expand = c(0.12, 0.03))+
  scale_y_continuous(breaks=seq(0,100,20),labels=waiver(),limits=c(0,100),expand = c(0.035, 0.0))+
  theme(legend.position = "right",panel.border=element_rect(fill=NA,color="black", size=0.7),axis.ticks = element_blank(), panel.grid.major.x = element_line(colour="grey90", size =0.001),panel.grid.major.y = element_line(colour="grey90", size =0.001),axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, colour = "black"),axis.text.y=element_text(size = 10.0, angle = 0, hjust = 1, colour = "black",face="italic"));g

