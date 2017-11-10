#So the idea here is to work out how rapidly synteny falls off as we consider further and further apart genes
#Using entropy to score the consistency of genes at a fixed distance
#Build a connected graph for each strain, then a distance matrix in which pairs can easily be checked.

library(igraph)

geneMatrix <- read.table("strains/clustering/new-groups/geneMatrix.txt",header=T,row.names=1)
geneTable <- read.table("strains/clustering/new-groups/geneTable.txt",header=T,stringsAsFactors=F)
geneFreq <- apply(geneMatrix,1,function(x) sum(x>0))
strains <- read.table("strains/sc1_strains.txt",header=F,stringsAsFactors=F)[,1]

if(import){
graphs <- list()
dms <- list()
for(strain in strains){
	cat(strain,"\n")

	strainGenes <- geneTable[geneTable$Strain==strain,]
	strainGenes[strainGenes$Group=="unique",]$Group <- strainGenes[strainGenes$Group=="unique",]$Gene
	strainParas <- names(which(table(strainGenes$Group)>1))
	replacements <- c()
	for(para in strainParas){
		newNames <- paste(strainGenes[strainGenes$Group==para,]$Group,1:sum(strainGenes$Group==para),sep="")
		strainGenes[strainGenes$Group==para,]$Group <- newNames
		replacements <- c(replacements,newNames)
	}

	strainNeighbours <- data.frame(gene1=character(),gene2=character())
	for(contig in unique(strainGenes$Contig)){
		contigGenes <- strainGenes[strainGenes$Contig==contig,]
		leftEnds <- apply(contigGenes[,5:6],1,min)
                contigGenes <- contigGenes[order(leftEnds),]
		contigNeighbours <- cbind(contigGenes$Group[-length(contigGenes$Group)],contigGenes$Group[-1])
		strainNeighbours <- rbind(strainNeighbours,contigNeighbours)
	}
	
	strainGraph <- graph_from_edgelist(as.matrix(strainNeighbours),directed=F)
	replacementNodes <- which(V(strainGraph)$name%in%replacements)
	V(strainGraph)$name[replacementNodes] <- sub("\\d+$","",V(strainGraph)$name[replacementNodes])
	graphs[[strain]] <- strainGraph

	strainDM <- distances(strainGraph)
	dms[[strain]] <- strainDM
}
}

#Scoring function
entropy <- function(a){
#Function to calculate the entropy of a vector
        at = table(a)
        at = at/sum(at)
        H = -sum(at*log2(at))
        return(H)
}

#Go through each ogg and score neighbour-pairs at different distances
dScores <- list()
for(d in 1:50){
	cat("Distance",d)
	scores = list()
	for(ogg in rownames(geneMatrix)[geneFreq>1]){
		oggNbs = list()
		for(strain in colnames(geneMatrix)[geneMatrix[ogg,]>0]){
#			nbRows <- dms[[strain]][ogg,]
#			nbs = c()
#			for(nbRow in nbRows){
#				rowNbs = c(names(nbRow==d),"Unknown")[1:2]
#				nbs = c(nbs,rowNbs)
#			}

			vs = which(V(graphs[[strain]])$name==ogg)
			if(length(vs)>0){
				nbs = names(unlist(ego(graphs[[strain]],d,vs,mindist=d)))
				if(!is.null(nbs)){
					nbs = c(unlist(nbs),rep("Unknown",length(vs)))[1:(length(vs)*2)]
					oggNbs[[strain]] = nbs
				}
			}
		}
		scores[[ogg]] = entropy(unlist(oggNbs))
	}
	dScores[[d]] <- unlist(scores)
	cat("\n")
}

scoresMatrix <- do.call(cbind,dScores)
write.table(scoresMatrix,"strains/synteny/syntenyScores.txt")
