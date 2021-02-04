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
refall <- refall[1:13,]     #siCtrl
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
  scale_fill_gradient(low="white",high="red",limits=c(0,1.0))+
  #scale_fill_gradient(low="white",high="mediumseagreen",limits=c(0,1.0))+
  labs(x = "",y = "") + 
  theme(panel.grid = element_blank())+#
  scale_x_discrete(labels=c("0","15","30","45","60","75","90","105","120","135","150","165","180"),breaks=c("0","15","30","45","60","75","90","105","120","135","150","165","180"),expand = c(0, 0))+ 
  scale_y_discrete(breaks=sort(names,decreasing = T),labels=waiver(),limits=sort(names,decreasing = T),expand = c(0, 0))+
  theme(legend.position = "bottom",axis.ticks = element_blank(), axis.text.x = element_text(size = 7.0, angle = 0, hjust = 0.5, colour = "black"),axis.text.y=element_text(size = 5.6, angle = 0, hjust = 1, colour = "black",face="italic"));g

#Draw a heatmap of simulation from simple model
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

smallest <- read.table("../Data/Simple_ERGsmallest.txt", sep = "\t")
smallest <- smallest[-which(rownames(smallest)=="NFKBIA"),]
smallest <- read.table("../Data/Simple_IRGsmallest.txt", sep = "\t")
smallest <- read.table("../Data/Simple_DRGsmallest.txt", sep = "\t")
names <- rownames(smallest)

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
  #
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

refall2 <- data.frame(matrix(unlist(refall),ncol=26,byrow=T))

refall3 <- refall2[,1:13]   #siCtrl
refall3 <- refall2[,14:26]  #siIkBa

refall3 <- t(refall3)
refall3 <- as.data.frame(refall3)
colnames(refall3) <- names

p <- colnames(refall3)
t <- colnames(refall3)
refall2$t <- seq(0,180,15)
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
df[which(df$name=="GADD45A"),2]="yellow"
df[which(df$name=="MAP3K8"),2]="yellow"
df[which(df$name=="PLEKHF2"),2]="yellow"
df[which(df$name=="PTGER4"),2]="yellow"
df[which(df$name=="RND1"),2]="yellow"
#df[which(df$name=="ADGRG6"),2]="yellow"
#df[which(df$name=="ANKRD18B"),2]="yellow"
#df[which(df$name=="ATP1B1"),2]="yellow"
#df[which(df$name=="BCL3"),2]="yellow"
#df[which(df$name=="BID"),2]="yellow"
#df[which(df$name=="BIRC3"),2]="yellow"
#df[which(df$name=="CDKN2B"),2]="yellow"
#df[which(df$name=="CLIC4"),2]="yellow"
#df[which(df$name=="EIF5"),2]="yellow"
#df[which(df$name=="FAM107B"),2]="yellow"
#df[which(df$name=="LYPD6B"),2]="yellow"
#df[which(df$name=="MAFF"),2]="yellow"
#df[which(df$name=="NCOA7"),2]="yellow"
#df[which(df$name=="NFKBIE"),2]="yellow"
#df[which(df$name=="RHOV"),2]="yellow"
#df[which(df$name=="S100A9"),2]="yellow"
#df[which(df$name=="SDC4"),2]="yellow"
#df[which(df$name=="SDCBP"),2]="yellow"
#df[which(df$name=="TICAM1"),2]="yellow"
simple_genes <- c("ABAT","ABTB2","AOX1","APOL2","ARID5B","BAZ1A","CHST15","CREBZF","CRIM1","DAXX","DRAM1","EFEMP1","FAM117A","FGD6","ITGB8","KYNU","LTB","OPTN","PHF23","PLEKHG3","PPP1R18","RELB","RFX5","SEMA4B","SGPP2","STAT5A","TRAF3")
except_genes <- c("CLK4","DDX58","HMGCR","IER5L","IL6ST","IRF2","LINC00052","RSRC2","SLC6A14","TUBD1")
df[-which(df$name %in% except_genes),2]="yellow"
df$x <- rep("all",dim(df)[1])

g <- ggplot(df,aes(as.factor(x),as.factor(name)))+
  geom_tile(aes(fill=group),color="grey50",stat="identity",na.rm = TRUE)+
  scale_fill_manual(values=c(black="grey85", yellow="yellow"))+
  labs(x = "",y = "") + 
  theme(panel.grid = element_blank())+#
  scale_x_discrete(labels=waiver(),breaks=group,expand = c(0, 0))+ 
  scale_y_discrete(breaks=sort(df$name,decreasing=T),labels=waiver(),limits=sort(df$name,decreasing=T),expand = c(0, 0))+
  theme(legend.position = "right",axis.ticks = element_blank(), axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, colour = "black"),axis.text.y = element_text(size = 12, angle = 0, hjust = 0.5, colour = "black"));g

#Draw lineplots of example gene
smallest <- read.table("../Data/Simple_ERGsmallest.txt", sep = "\t")
smallest <- smallest[-which(rownames(smallest)=="NFKBIA"),]
smallest <- read.table("../Data/Simple_IRGsmallest.txt", sep = "\t")
smallest <- read.table("../Data/Simple_DRGsmallest.txt", sep = "\t")
names <- rownames(smallest)

d <- names[r]

stackkari <- c()
stack <- c()
dfde <- c()
file <- paste("../Optimized_parameters/Simple_model/ERGsubcluster2/rep1/2_",d,"_subplex_rep1",".txt",sep="")
datafile <- read.table(file,sep="\t",skip=1) 
dfde <- datafile[,2:5]
dfde <- as.data.frame(dfde)
colnames(dfde) <- c("RMSD","kdeg","KD","tau")

kdeg <- dfde[smallest[r,1],2]
kd <- dfde[smallest[r,1],3]
tau <- dfde[smallest[r,1],4]
rmsd <- dfde[smallest[r,1],1]

#stackkari <- c()
#stack <- c()
#dfde <- c()
#file <- paste("../Optimized_parameters/Simple_model/ERGsubcluster2/rep2/2_",d,"_subplex_rep2",".txt",sep="")
#datafile <- read.table(file,sep="\t",skip=1) 
#dfde <- datafile[,2:5]
#dfde <- as.data.frame(dfde)
#colnames(dfde) <- c("RMSD","kdeg","KD","tau")
#
#kdeg <- dfde[smallest[r,2],2]
#kd <- dfde[smallest[r,2],3]
#tau <- dfde[smallest[r,2],4]
#rmsd <- dfde[smallest[r,2],1]

count <- read.csv("../Data/160712_ikbakd.csv")
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


bothsignal <- rbind(wtsignal,kdsignal)
bothsignal <- as.data.frame(bothsignal)

df <- rescale(bothsignal[,2], to=c(2,100))
df <- data.frame(time=rep(seq(0,10800,900),2),mrna=df)
df <- df[1:13,]                                         

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

bothsignal <- rbind(wtsignal,kdsignal)
bothsignal <- as.data.frame(bothsignal)

df <- rescale(bothsignal[,2], to=c(2,100))
df <- data.frame(time=rep(seq(0,10800,900),2),mrna=df)
df <- df[14:26,]                                         

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



#count <- read.csv("../Data/160617_ikbakd.csv")
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
#df <- rescale(wtsignal[1:4,2], to=c(2,100))
#unit <- (wtsignal[4,2]-wtsignal[1,2])/(df[4]-df[1])
#wts <- c()
#for(i in 5:13){
#  kari <- wtsignal[i,2]-wtsignal[1,2]
#  final <- kari/unit+2
#  wts[i] <- final
#}
#wts <- na.omit(wts)
#wts <- c(df,wts)
#
#kds <- c()
#for(i in 1:dim(kdsignal)[1]){
#  kari <- kdsignal[i,2]-wtsignal[1,2]
#  final <- kari/unit+2
#  kds[i] <- final
#}
#
#df <- data.frame(time=seq(0,10800,900),mrna=wts)
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
#df <- data.frame(time=seq(0,10800,900),mrna=kds)
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
#times <- seq(-3600,0,1)
#df <- rep(nfkb_kd[1,2],length(times))
#fp <- pchipfun(times, df) 
#sigimp <- approxfun(times, df, rule = 2)
#steadykd <- sigimp





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
  scale_colour_manual(values = c(wt="red",kd="green4"))+
  scale_x_continuous(breaks=seq(0,10800,1800),labels=seq(0,180,30))+
  scale_y_continuous(breaks=c(0.0,0.2,0.4,0.6,0.8,1.0),labels=waiver(),limits=c(0.0,1.0))+
  scale_linetype_manual(values=c(solid="solid", dashed="dotted"))+
  scale_alpha_manual(values=c(solid=0.7, dashed=1.0))+
  theme(panel.border=element_rect(fill=NA,color="black", size=0.6),plot.title = element_text(hjust = 0.5,size=15.0,margin=margin(t = 0, r = 0, b = 0.1, l = 0)),panel.grid.major.x = element_line(colour="grey90", size =0.001),panel.grid.major.y = element_line(colour="grey90", size =0.001),legend.position = "none",title = element_text(size=7),legend.text = element_text(size=7),axis.text.x = element_text(size=10),axis.text.y = element_text(size=10))+
  geom_line(lwd=2.0,aes(color=condition));g
