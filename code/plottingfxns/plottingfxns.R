library(reshape2)
library(viridis)


MatrixPrep <- function(v){
  # Melt a matrix and display upright, 1,1 in the corner
  numClusters = nrow(v)
  v <- melt(t(v)[,numClusters:1])
  v
}

MakeColormap <- function(v){
  make0 <- seq(min(v$value),max(v$value),length = 1000)
  make0 <- which(diff(sign(make0)) != 0)
  customcmap <- viridis(1000,option = 'plasma')
  customcmap[(make0-1):(make0+1)] <- "grey50"
  customcmap
}

fliplr <- function(x){
  x <- x[length(x):1]
  return(x)
}

TP.beta.plot <- function(v,clusterColors,title){
	numClusters = nrow(v)
	v <- melt(t(v)[,numClusters:1])

	make0 <- seq(min(v$value),max(v$value),length = 1000)
	make0 <- which(diff(sign(make0)) != 0)
	customcmap <- viridis(1000,option = 'plasma')
	customcmap[(make0-1):(make0+1)] <- "grey50"
	p <- ggplot() + geom_tile(aes(x = v$Var1,y = v$Var2, fill = v$value)) +   scale_fill_gradientn(colours = customcmap,name ="") +
	theme_bw() + xlab('') + ylab('') +
	theme(axis.text.x = element_text(size=8,angle = 90,hjust = 0.95,vjust = 0.5)) +
  theme(axis.title.x = element_blank()) + theme(axis.ticks.x = element_blank()) + #size = 0.2,colour = "black")
  scale_x_discrete(limits = 1:numClusters, breaks = 1:numClusters, labels = list(clusterNames), expand = c(0,0)) + 
  scale_y_discrete(limits = numClusters:1, breaks = numClusters:1, labels = list(clusterNames), expand = c(0,0)) +
  theme(axis.title.y = element_blank()) + theme(axis.ticks.y = element_blank()) + 
  theme(axis.text.y = element_text(size=8)) +
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 8), legend.position = "bottom",
        legend.background = element_rect(fill=0,size = 0.1),legend.key.width = unit(0.2,"in"), legend.key.height = unit(0.05,"in")) +
  	ggtitle(title) + theme(text = element_text(size = 8),plot.title = element_text(hjust = 0.5,size = 8)) + coord_fixed(ratio = 1)

if(numClusters == 5){
	p <- p + theme(axis.text.x = element_text(color = clusterColors),axis.text.y = element_text(color = clusterColors))
}

}

p.signif.matrix <- function(p){
  # take matrix of p values and make asterisks
  #  : p > 0.05 (blank)
  # *: p <= 0.05

  p.new <- matrix(data = '',nrow = nrow(p),ncol=ncol(p))
  p.new[p > 0.05] <- ''
  p.new[p < 0.05] <- '*'
  return(p.new)
  
}

plot.beta.matrix <- function(b.mat,p.mat,clusterColors,clusterNames,title,bmin=min(b.mat),bmax=max(b.mat)){
    # plots matrix of beta weights associated with transitions between clusters
    # b.mat: matrix of betas
    # p.mat: matrix of p-values
    # clusterColors: character vector of colors for each cluster
    # clusterNames: character vector of names for each cluster
    # title: string containing title
    # bmin and bmax: numbers signifying min/max values for color axis

    melted_betas <- melt(t(b.mat))
    melted_betas$Var2 <- fliplr(melted_betas$Var2)  
    p <- p.signif.matrix(p.mat)
    melted_p <- melt(t(p))
    melted_p$Var2 <- fliplr(melted_p$Var2)
    p1 <- ggplot() + 
      geom_tile(data = melted_betas, aes(x=Var1, y=Var2, fill=value)) + xlab("") + ylab("") +
      geom_text(data = melted_p, aes(x=Var1,y=Var2,label=value),color='white',size=2.5) +
      scale_fill_viridis(option='plasma',limits=c(bmin,bmax)) +
      ggtitle(title) +
      scale_y_discrete(limits=fliplr(clusterNames),labels = fliplr(clusterNames),expand=c(0,0)) + 
      scale_x_discrete(limits=clusterNames,labels = clusterNames,expand=c(0,0),position='bottom') + coord_fixed() +
      theme_classic() + theme(text=element_text(size = 6), legend.key.size = unit(0.1,'inches'),
            legend.position = 'right',
            axis.line = element_blank(),
            axis.text.x=element_text(color=clusterColors,angle=90,size=8),
            axis.text.y = element_text(color=fliplr(clusterColors),size=8),
            plot.title = element_text(size=8,hjust=0.5,face = 'bold'),
            plot.margin = unit(c(0, 0, 0, 0), "cm"),
            axis.ticks = element_blank())
    return(p1)
}

p.violin.nptest <- function(val,dist,xlab='',ylab='',col='blue'){
  # val: observed value
  # dist: null distribution
  # xlab: name of quantity for x-axis tick label

  # plots violin plot of distribution, with val as a point
  df <- data.frame(x=rep(xlab,length(dist)),y=dist)
  p <- ggplot(df) + geom_violin(aes(x=x,y=dist),fill=col,alpha=0.5) + geom_point(aes(x=xlab,y=val)) +
  xlab('') + theme_classic() + theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) + theme(text=element_text(size=8))
  return(p)
}

getClusterColors <- function(k){
  # retrieve cluster colors for pre-defined cluster orders at different values of k
  # current I have only defined specific colors for k = 5 and k = 6
  # all other k's will return vectors of all black

  if(k==5){
    clusterColors <- c("1"="#AB484F","2"="#591A23", "3"="#AA709F","4"="#527183","5"="#7E874B")
  } else if(k == 6){
    # 4 = 99628E
    clusterColors <- c("1"="#AB484F","2"="#591A23", "3"="#D98FCA","4"="#945F8A","5"="#C6CDF7","6"="#7294D4")
  } else { clusterColors <- rep('black',k)}

  return(clusterColors)
}

### Control plotting ###

p.Tsweep <- function(t.rng,E.brain,E.null,TP,ED,col,ttl,method='pearson',leg=TRUE){
  # plot t.rng vs. correlation between control energy and transition probability
  # also show correlation between trans probs and null energies, and trans probs and distance
  # t.rng: vector of t values tested
  # E.brain: n_transitions-by-T matrix of energies
  # E.null: list of n_transitions-by-T-by-nperms matrix of energies. names of list elements are legend labels of null
  # TP: transition probabilities
  # ED: inter-state euclidean distance

  nT <- dim(E.brain)[2]
  r.brain <- cor(E.brain,TP,method=method)
  r.ed <- rep(cor(ED,TP,method=method),nT)
  df1 <- data.frame(T=t.rng,r.brain=r.brain,r.ed=r.ed)

  p <- ggplot() + 
  geom_line(aes(x=df1$T,y=df1$r.brain),color=col,size=0.25) + 
  geom_line(aes(x=df1$T,y=df1$r.ed),color='grey50',linetype='dashed',size=0.25)
  if(!is.null(E.null)){
      n.nulls <- length(E.null)      
      r.null <- lapply(E.null, function(E.null.i) matrix(sapply(1:nT, function(T) cor(E.null.i[,,T],TP,method=method)),nrow=dim(E.null.i)[2],ncol=nT))
      r.null.min <- lapply(r.null, function(X) sapply(1:nT, function(T) min(X[,T])))
      r.null.max <- lapply(r.null, function(X) sapply(1:nT, function(T) max(X[,T])))
      df2 <- data.frame(r.null=unlist(lapply(r.null,colMeans)),
        r.null.min=unlist(r.null.min),r.null.max=unlist(r.null.max),
        T.nulls=rep(t.rng,n.nulls),null.fill=rep(names(E.null),each=nT))
  p <- p + geom_ribbon(data=df2,aes(x=T.nulls,ymin=r.null.min,ymax=r.null.max,fill=null.fill),alpha=0.5)+ 
  geom_line(data=df2,aes(x=T.nulls,y=r.null,color=null.fill),size=0.25) +
  scale_fill_manual(limits = names(E.null),values= brewer.pal('Set3',n=10)[1:n.nulls])+
  scale_color_manual(limits = names(E.null),values= brewer.pal('Set3',n=10)[1:n.nulls])
  }
  p <- p + scale_y_continuous(limits=c(-1,1),breaks=seq(-1,1,length.out = 5)) +
  xlab('T') + ylab('r(TP,E)') + ggtitle(ttl) + 
    theme_classic() + theme(text = element_text(size = 8)) + 
    theme(plot.title = element_text(size=8,hjust=0.5,face = "bold")) +
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
  if(!leg){p <- p + theme(legend.position = 'none')}
return(p)
}

p.xy <- function(x,y,xlab,ylab,ttl='',col='black'){
  df <- data.frame(x=x,y=y)
  r.text <- paste('r_p =',signif(cor(x,y),2),'r_s =',signif(cor(x,y,method='spearman'),2))
  p <- ggplot(df) + geom_point(aes(x=x,y=y),color=col) + geom_smooth(aes(x=x,y=y),fill=col,color=col,method='lm') + 
    xlab(xlab) + ylab(ylab) + ggtitle(ttl) + 
    annotate("text",size = 2, x = min(df$x),y = 0, label = r.text,hjust=0) +
      theme_classic() + theme(text = element_text(size = 8)) + 
      theme(plot.title = element_text(size=8,hjust=0.5,face = "bold")) +
      theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) + theme(legend.position = 'none')
  return(p)
}

scatter.ETP <- function(SC,TP,col,ttl='',plabel.df=NULL,xlabel='Transition Energy (rank)',ylabel='Transition Probability (rank)',method='pearson'){
  # SC: x -values of plot, here control energies
  # TP: y-values of plot, here transition probabilities
  # plabel.df: data frame with p-values for different permutation tests, with name of test in column named test, label in column named plabel

  TPmin <- min(TP)
  TPmax <- max(TP)
  SCmin <- min(SC)
  SCmax <- max(SC)
  r.text <- paste('r =',signif(cor(SC,TP,method=method),2))
  if(!is.null(plabel.df)){
  r.text <- paste(c(r.text,
    sapply(1:nrow(plabel.df), function(R) paste(plabel.df$test[R],plabel.df$plabel[R]))),collapse='\n')
  }
  p <- ggplot() + geom_point(aes(x = SC,y = TP),color = col,size = 1, alpha = 0.75,stroke = 0) + 
      geom_smooth(aes(x = SC,y = TP),color = col, fill = col,method = 'lm',size=0.5) +
      scale_y_continuous(breaks= seq(0,TPmax,length.out=5)) + scale_x_continuous() +
      annotate("text",size = 2, x = -Inf,y = -Inf, label = bquote(.(r.text)),hjust=-0.1,vjust=-0.15) +
      xlab(xlabel) + ylab(ylabel) + ggtitle(ttl) + 
      theme_classic() + theme(text = element_text(size = 8),axis.title=element_text(size=6)) + 
      theme(plot.title = element_text(size=8,hjust=0.5,face = "bold")) +
      theme(legend.position = 'none') + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
  return(p) 
}

collapse.columns <- function(df,cnames=colnames(df),groupby=NULL){
  # INPUTS:
  # df: dataframe
  # cnames: column names to perform operation on, default to all columns
  # groupby: column name to group variables by, treated separately from variables in cnames
  
  # OUTPUTS:
  # df.new: dataframe with 2 columns:
  # values: all columns in cnames vertically concatenated. 
  # names: elements of cnames corresponding to rows
  # group: groups of observations in df for each variable in cnames
  
  df.names <- do.call('cbind',lapply(cnames, function(n) rep(n,nrow(as.matrix(df[,cnames])))))  
  df.new <- data.frame(values = as.vector(as.matrix(df[,cnames])),names=as.vector(df.names))
  if(!is.null(groupby)){
    df.grp <- do.call('cbind',lapply(cnames,function(n) df[,groupby]))
    df.new$group <- as.vector(df.grp)
  }
  return(df.new)
}

matrix.to.df <- function(m,dnames){
  # INPUTS:
  # m: matrix
  # dnames: 2 element list specifying dimnames for m
  #
  # OUTPUTS:
  # df:
  dimnames(m) <- dnames
  return(as.data.frame(m))

}