library(ggplot2)

    count <- read.csv("../Data/160712_ikbakd.csv")
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
    wtsignal$group <- rep("none",dim(wtsignal)[1])
    wtsignal$time <- as.numeric(wtsignal$time)
    atac <- c(0,30,75,120)
    wtsignal[which(wtsignal$time %in% atac),3] <- c("atac")
    wtrep1signal <- wtsignal
    
    count <- read.csv("../Data/Ando_160617_ikbakd.csv")
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
    wtsignal$group <- rep("none",dim(wtsignal)[1])
    wtsignal$time <- as.numeric(wtsignal$time)
    atac <- c(0,30,75,120)
    wtsignal[which(wtsignal$time %in% atac),3] <- c("atac")
    wtrep2signal <- wtsignal
    
    wtsignal <- wtrep1signal
    two <- cbind(wtrep1signal$nfkb,wtrep2signal$nfkb)
    two <- as.data.frame(two)
    two$mean <- apply(two,1,mean)
    wtsignal$nfkb <- two$mean
    
    count <- read.csv("../Data/160712_ikbakd.csv")
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
    kdsignal$group <- rep("none",dim(kdsignal)[1])
    kdsignal$time <- as.numeric(kdsignal$time)
    atac <- c(0,30,75,120)
    kdsignal[which(kdsignal$time %in% atac),3] <- c("atac")
    kdrep1signal <- kdsignal
    
    count <- read.csv("../Data/Ando_160617_ikbakd.csv")
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
    kdsignal$group <- rep("none",dim(kdsignal)[1])
    kdsignal$time <- as.numeric(kdsignal$time)
    atac <- c(0,30,75,120)
    kdsignal[which(kdsignal$time %in% atac),3] <- c("atac")
    kdrep2signal <- kdsignal
    
    kdsignal <- kdrep1signal
    two <- cbind(kdrep1signal$nfkb,kdrep2signal$nfkb)
    two <- as.data.frame(two)
    two$mean <- apply(two,1,mean)
    kdsignal$nfkb <- two$mean
    
    g <- ggplot(kdsignal,aes(x=time,y=nfkb,color=group,size=group))+
      theme_bw(base_size = 16)+
      theme(panel.grid = element_blank())+
      labs(x="",y="",title="",colour="")+
      #geom_line(lwd=2.2,color="#DC143C")+
      geom_line(lwd=2.2,color="seagreen")+
      scale_x_continuous(breaks=seq(0,180,15),labels=waiver(),limits=c(0,180))+
      scale_y_continuous(breaks=seq(0.25,0.5,0.05),labels=waiver(),limits=c(0.25,0.5))+
      geom_point(aes(shape=group,color=group))+
      scale_color_manual(values=c(none="seagreen",atac="black"))+
      scale_shape_manual(values = c(none=20,atac=21))+
      scale_size_manual(values=c(none=0.0,atac=10.0))+ 
      theme(panel.border=element_rect(fill=NA,color="black", size=1.3),panel.grid.major.x = element_line(colour="grey", size =0.7),panel.grid.major.y = element_blank(),legend.position = "none",title = element_blank(),legend.text = element_blank(),axis.text.x = element_text(size=12),axis.text.y = element_text(size=12));g

    
