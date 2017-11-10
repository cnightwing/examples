

coords <- read.table("strains/alignment/coords/C10-C11.coords",stringsAsFactors=F)
coords <- cbind(coords,coords[,4]<coords[,3]) #Reverse-complement?
n = nrow(coords)

geneTable <- read.table("strains/clustering/new-groups/geneTable.txt",header=T,stringsAsFactors=F)

if(0){
#Reorder the alignments starting with the largest and putting many-to-one matches together THIS IS HARDER THAN I THOUGHT IT WOULD BE ARGH
newCoords <- coords[0,]
unused <- rep(TRUE,n)
while(sum(unused)>0){
	longest <- which(coords[,5]==max(coords[unused,5]))
	topDone = FALSE
	botDone = FALSE
	#Check the top contig
	tempCoords <- coords[unused & coords[,8]==coords[longest,8],]
	if(nrow(tempCoords)==1){
		topDone = TRUE
		tempCoords <- coords[unused & coords[,9]==coords[longest,9],]
		if(nrow(tempCoords)==1){
			botDone=TRUE
		}
	}
	while(!topDone | !botDone){

	}
	
	
		
	
	topCount <- sum(coords[unused,8]==tempCoords[,8])
	botCount <- sum(coords[unused,9]==tempCoords[,9])
	cat(longest,topCount,botCount,"\n")
	if(topCount == 1 & botCount == 1){
		unused[longest] = FALSE
	}
	while(topCount > 1 | botCount > 1){
		cat("In the loop\n")
		if(topCount > 1 & botCount > 1){
			cat("WTF?")
		}
		if(topCount>botCount){
			tempCoords <- coords[unused & coords[,8]%in%tempCoords[,8],]
			tempCoords <- tempCoords[order(tempCoords[,1]),]
			unused[coords[,8]==tempCoords[,8]] = FALSE
		}else{
			tempCoords <- coords[unused & coords[,9]%in%tempCoords[,9],]
			tempCoords <- tempCoords[order(tempCoords[,3]),]
			unused[coords[,9]==tempCoords[,9]] = FALSE
		}
		topCount <- sum(coords[unused,8]%in%tempCoords[c(1,nrow(tempCoords)),8])
		botCount <- sum(coords[unused,9]%in%tempCoords[c(1,nrow(tempCoords)),9])
	}
	newCoords <- rbind(newCoords,tempCoords)
	line <- readline()
}
}

l1 <- read.table("strains/alignment/C10-lengths.txt",row.names=1)	
n1 <- nrow(l1)
s1 <- l1
s1[,1] <- 0
e1 <- l1
e1[,1] <- 0
r1 <- l1
r1[,1] <- FALSE

l2 <- read.table("strains/alignment/C11-lengths.txt",row.names=1)
n2 <- nrow(l2)
s2 <- l2
s2[,1] <- 0
e2 <- l2
e2[,1] <- 0
r2 <- l2
r2[,1] <- FALSE

if(0){
#Correct for "<" characters
for(i in 1:n){
	if(grepl("<",coords[i,8])){
		coords[i,8] <- substr(coords[i,8],1,nchar(coords[i,8])-1)
	}
	if(grepl("<",coords[i,9])){
		coords[i,9] <- substr(coords[i,9],1,nchar(coords[i,9])-1)
	}
}
}

#Redo RC coordinates
for(i in 1:n){
	if(coords[i,3]>coords[i,4]){
		coords[i,3] <- l2[coords[i,9],1]+1-coords[i,3]
		coords[i,4] <- l2[coords[i,9],1]+1-coords[i,4]
	}
}

lastTop <- ""
lastBot <- ""

for(i in 1:n){
	if(coords[i,8]==lastTop){
		anchor = s1[coords[i,8],1]
		offset = coords[i,1]-coords[i,3]
		s2[coords[i,9],1] <- anchor+offset
		e2[coords[i,9],1] <- s2[coords[i,9],1]+l2[coords[i,9],1]
		r2[coords[i,9],1] <- coords[i,10]
	}else if(coords[i,9]==lastBot){
		anchor = s2[coords[i,9],1]
		offset = coords[i,3]-coords[i,1]
		s1[coords[i,8],1] <- anchor+offset
                e1[coords[i,8],1] <- s1[coords[i,8],1]+l1[coords[i,8],1]
		r1[coords[i,8],1] <- coords[i,10]
	}else{
		anchor = max(0,e1[lastTop,1],e2[lastBot,2],na.rm=T)+1000
		offset = coords[i,1]-coords[i,3]
		if(offset<0){
			anchor = anchor-offset
		}
		s1[coords[i,8],1] <- anchor
		e1[coords[i,8],1] <- s1[coords[i,8],1]+l1[coords[i,8],1]
		s2[coords[i,9],1] <- anchor+offset
		e2[coords[i,9],1] <- s2[coords[i,9],1]+l2[coords[i,9],1]
		r2[coords[i,9],1] <- coords[i,10]
	}
	lastTop = coords[i,8]
	lastBot = coords[i,9]
}

pdf("strains/alignment/C10-C11.pdf",width=200,height=1)
par(mar=c(0,0,0,0),xaxs="i",yaxs="i")
plot(0,type="n",xlab="",ylab="",xlim=c(0,6e6),ylim=c(0,2))

for(i in 1:nrow(s1)){
	rect(s1[i,1],1,e1[i,1],1.3,border=1,col=2+r1[i,1],lwd=0.1)
	text(s1[i,1]+200,1.3,paste(i,rownames(s1)[i]),pos=4,srt=90,offset=0,cex=0.1)
	genes <- geneTable[geneTable[,4]==rownames(s1)[i],]
	if(nrow(genes)>0){
		arrows(s1[i,1]+genes$Start,1.05,s1[i,1]+genes$End,1.05,length=0.005,lwd=0.1)
		text(s1[i,1]+apply(genes[,c("Start","End")],1,mean),1.1,genes$Group,pos=4,srt=90,offset=0,cex=0.1)
	}
}

for(i in 1:nrow(s2)){
	rect(s2[i,1],0.7,e2[i,1],1,border=1,col=2+r2[i,1],lwd=0.1)
	text(s2[i,1]+200,0.7,paste(i,rownames(s2)[i]),pos=2,srt=90,offset=0,cex=0.1)
	genes <- geneTable[geneTable[,4]==rownames(s2)[i],]
	if(nrow(genes)>0){
		arrows(s2[i,1]+genes$Start,0.95,s2[i,1]+genes$End,0.95,length=0.005,lwd=0.1)
		text(s2[i,1]+apply(genes[,c("Start","End")],1,mean),0.9,genes$Group,pos=2,srt=90,offset=0,cex=0.1)
	}
}

if(0){
R = 0
anchor = 0
for(i in 1:n){
	if(!p1[coords[i,8],1] & !p2[coords[i,9],1]){
		anchor = R+1000
	}
	cat(i,":",anchor)
	if(!p1[coords[i,8],1]){
		topStart <- anchor-coords[i,1]
		topEnd <- anchor-coords[i,1]+l1[coords[i,8],1]
		rect(topStart,1,topEnd,1.3,border=NA,col=2)
		text(topStart+200,1.3,paste(i,coords[i,8]),pos=4,srt=90,offset=0,cex=0.1)

		genes <- geneTable[geneTable[,4]==coords[i,8],]
		if(nrow(genes)>0){
			arrows(topStart+genes$Start,1.05,topStart+genes$End,1.05,length=0.005,lwd=0.1)
			text(topStart+apply(genes[,c("Start","End")],1,mean),1.1,genes$Group,pos=4,srt=90,offset=0,cex=0.1)
		}

		p1[coords[i,8],1] = TRUE
		cat(" Top",topStart-anchor,":",topEnd-anchor)
	}
	if(!p2[coords[i,9],1]){
		botStart <- anchor+coords[i,1]-coords[i,3]
		botEnd <- anchor+coords[i,1]-coords[i,3]+l2[coords[i,9],1]
		rect(botStart,0.7,botEnd,1,border=NA,col=3-coords[i,10])
		text(botStart,0.7,paste(i,coords[i,9]),pos=2,srt=90,offset=0,cex=0.1)

		p2[coords[i,9],1] = TRUE
		cat(" Bottom",botStart-anchor,":",botEnd-anchor)
	}
	R = max(topEnd,botEnd)
	cat("\n")
}
}
segments(0,1,6e6,1)
dev.off()
