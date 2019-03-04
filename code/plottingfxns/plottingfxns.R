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

