#Script to analyse trees produced from the output of OrthoMCL; this version has had redundant functions removed and others modified

library(ape)

loadTrees <- function(n){
	#Load amino-acid and nucleotide trees into R (not very flexible)
	aaTrees <- list()
	ntTrees <- list()
	#8955 groups
	for(i in 1:n){
		if(i%%45==0){
			cat(".")
		}
		if(file.exists(paste("strains/clustering/trees/SC1",i+999,".faa.tree",sep=""))){
			aaTrees[[i]] <- read.tree(paste("strains/clustering/trees/SC1",i+999,".faa.tree",sep=""))
			ntTrees[[i]] <- read.tree(paste("strains/clustering/trees/SC1",i+999,".ffn.tree",sep=""))
		}
	}
	cat("\n")
	return(list(aaTrees=aaTrees,ntTrees=ntTrees))
}

ttt <- function(tree){
	#Calculate the largest tip-to-tip distance in a tree
	if(max(tree$edge.length)>0){
		distance <- max(cophenetic(tree))
	}else{
		distance <- 0
	}
	return(distance)
}

ntips <- function(tree){
	#Count the number of tips in a tree
	return(length(tree$tip.label))
}

cutTree <- function(tree,edge){
	#Reformulation to cut a tree at a given edge by taking the 'lower' node of that edge and separating the tree of all tips beneath it from the tree of all other tips; this version stores no information about the cut
	node <- tree$edge[edge,2]
	Ntips <- ntips(tree)
	
	subTrees <- list()
	if(node<=Ntips){
		tips1 <- node
	}else{
		tips1 <- prop.part(tree)[[node-Ntips]]
	}
	tips2 <- setdiff(1:Ntips,tips1)

	if(length(tips1)>1){
		subTrees[[1]] <- drop.tip(tree,tips2,rooted=F)
		if(length(subTrees[[1]]$edge.length)>2){
			subTrees[[1]] <- unroot(subTrees[[1]])
		}
	}else{
		subTrees[[1]] <- list(tip.label=tree$tip.label[tips1],edge.length=0)
	}

	if(length(tips2)>1){
		subTrees[[2]] <- drop.tip(tree,tips1,rooted=F)
		if(length(subTrees[[2]]$edge.length)>2){
			subTrees[[2]] <- unroot(subTrees[[2]])
		}
	}else{
		subTrees[[2]] <- list(tip.label=tree$tip.label[tips2],edge.length=0)
	}

	return(subTrees)
}

identifyStrains <- function(tree){
	#Identifies the strains associated with each tip in a tree
	strains <- unlist(regmatches(tree$tip.label,gregexpr("[A-Z]\\d\\d?",tree$tip.label)))
	return(strains)
}

countPresentStrains <- function(tree){
	#Counts the frequency of presents strains in a tree
	strains <- identifyStrains(tree)
	count <- table(strains)
	return(count)
}

findParaTips <- function(tree){
	#Determines which tips are paralogous within a tree
	strains <- identifyStrains(tree)
	count <- countPresentStrains(tree)
	paralog <- as.numeric(count[strains])>1
	return(paralog)
}

findParaGroups <- function(tree){
	#Determine the paralogous groups of a tree
	strains <- identifyStrains(tree)
        count <- table(strains)
        paraStrains <- names(count[count>1])
        paraGroups <- list()
	if(length(paraStrains)>0){
	        for(i in 1:length(paraStrains)){
	                paraGroups[[i]] <- which(strains==paraStrains[i])
	        }
	}
        return(paraGroups)
}

findParaEdges <- function(tree){
	#Determine
        paraGroups <- findParaGroups(tree)
        paraPaths <- list()
	if(length(paraGroups)>0){
	        for(i in 1:length(paraGroups)){
	                paraPaths[[i]] <- which.edge(tree,paraGroups[[i]])
	        }
	}
	paraPaths <- table(unlist(paraPaths))
	paraEdges <- rep(0,length(tree$edge.length))
	paraEdges[as.numeric(names(paraPaths))] <- paraPaths
        return(paraEdges)
}

treeCensus <- function(tree){
	#Count the number of strains present, with paralogs, and their average representation
	count <- countPresentStrains(tree)
	present <- sum(count>0)
	plural <- sum(count>1)
	representation <- mean(count)
	return(list(present=present,plural=plural,representation=representation))
}

cutLongestPPBranch <- function(tree){
	#Cut a tree at the longest branch on a paralog path
	paraEdges <- findParaEdges(tree)
	if(max(paraEdges)>0){
		longestParaEdges <- which(tree$edge.length==max(tree$edge.length[paraEdges>0]))
		longestParaEdge <- longestParaEdges[order(paraEdges[longestParaEdges],decreasing=T)][1]
		length <- tree$edge.length[longestParaEdge]
		ordinal <- which(order(tree$edge.length,decreasing=T)==longestParaEdge)
		subTrees <- cutTree(tree,longestParaEdge)
		return(list(subTrees=subTrees,length=length,ordinal=ordinal))
	}else{
		return(list(subTrees=list(tree),length=0,ordinal=0))
	}
}

cutMostPPBranch <- function(tree){
	#Cut a tree at the branch on the most paralog paths UNUSED
	paraEdges <- findParaEdges(tree)
	if(max(paraEdges)>0){
		mostParaEdges <- which(paraEdges==max(paraEdges))
		mostParaEdge <- mostParaEdges[order(tree$edge.length[mostParaEdges],decreasing=T)][1]
		length <- tree$edge.length[mostParaEdge]
		ordinal <- which(order(tree$edge.length,decreasing=T)==mostParaEdge)
		subTrees <- cutTree(tree,mostParaEdge)
		RdR <- findRdR(tree,mostParaEdge)
		return(list(subTrees=subTrees,length=length,ordinal=ordinal,RdR=RdR))
	}else{
		return(list(subTrees=list(tree),length=0,ordinal=0,RdR=list(q=0,R=0,dR=0,iter=0)))
	}
}

removeParalogs <- function(tree,mode="long",cutoff=0,plot=F){
	#Remove all paralogs from the tree by repeatedly cutting the longest branch that lies on a paralog path, recording information about each cut
	trees <- list(tree)
	lengths <- vector()
	ordinals <- vector()
	presents <- vector()
	plurals <- vector()
	representations <- vector()
	losses <- vector()
	resolutions <- vector()

	while(1){
		noTrees <- length(trees)
		longestParaEdgesByTree <- vector()
		for(tree in trees){
			paraEdges <- findParaEdges(tree)
			if(max(paraEdges)>0){
				longestParaEdges <- which(tree$edge.length==max(tree$edge.length[paraEdges>0]))
		                longestParaEdge <- longestParaEdges[order(paraEdges[longestParaEdges],decreasing=T)][1]
				longestParaEdgesByTree <- c(longestParaEdgesByTree,tree$edge.length[longestParaEdge])
			}else{
				longestParaEdgesByTree <- c(longestParaEdgesByTree,0)
			}
		}
		if(max(longestParaEdgesByTree)>cutoff){
			treeToCut <- which.max(longestParaEdgesByTree)
			census <- treeCensus(trees[[treeToCut]])
			if(mode=="long"){
				results <- cutLongestPPBranch(trees[[treeToCut]])
			}
			if(mode=="most"){
				results <- cutMostPPBranch(trees[[treeToCut]])
			}
			lengths <- c(lengths,results$length)
			ordinals <- c(ordinals,results$ordinal)
			presents <- c(presents,census$present)
			plurals <- c(plurals,census$plural)
			representations <- c(representations,census$representation)
			census1 <- treeCensus(results$subTrees[[1]])
			census2 <- treeCensus(results$subTrees[[2]])
			losses <- c(losses,census$present-max(census1$present,census2$present))
			resolutions <- c(resolutions,(census$present*(census$representation-1))-sum(census1$present*(census1$representation-1),census2$present*(census2$representation-1)))
			if(plot){
				par(mfrow=c(1,3))
				par(oma=c(4,0,0,0))
				plotCutTree(trees[[treeToCut]],results$subTrees,highlight=which(trees[[treeToCut]]$edge.length==sort(trees[[treeToCut]]$edge.length,decreasing=T)[results$ordinal]))
				plotTrees(results$subTrees)
				mtext(paste(format(lengths[length(lengths)],3)," (",format(ordinals[length(ordinals)],1),"): ",format(losses[length(losses)],1)," lost, ",format(resolutions[length(resolutions)],1)," resolved",sep=""),1,outer=T,cex=2)
                        }
			trees <- c(trees[setdiff(1:length(trees),treeToCut)],results$subTrees)
		}
		if(length(trees)==noTrees){
			break
		}
	}
	return(list(subTrees=trees,lengths=lengths,ordinals=ordinals,presents=presents,plurals=plurals,representations=representations,losses=losses,resolutions=resolutions))
}

calculateP <- function(s,t,n){
	#Calculate the probability of n genes ending up in both halves of a pool of s strains when divided into pools of t strains and s-t strains UNUSED
	pbar <- (choose(t,n)+choose(s-t,n))/choose(s,n)
	return(1-pbar)
}

calculateRdR <- function(tree,edge,q){
	#Calculate the probability of seeing a given division of strains if a tree is cut on a certain edge, with q as the probability of orthology UNUSED
	tips <- ntips(tree)
	count <- countPresentStrains(tree)

        subTrees <- cutTree(tree,edge)
	subTips <- unlist(lapply(subTrees,ntips))
	subCounts <- lapply(subTrees,countPresentStrains)
	bstrains <- intersect(names(subCounts[[1]]),names(subCounts[[2]]))
	bbarstrains <- setdiff(union(names(subCounts[[1]]),names(subCounts[[2]])),bstrains)
	#bbarstrains <- setdiff(names(count[count>1]),bstrains)
	bbar <- length(bbarstrains)

	ps <- sapply(count,function(n) calculateP(tips,ntips(subTrees[[1]]),n))
	bsum <- sum(log(1+q*(1-ps[bstrains])/ps[bstrains]))
	bbarsum <- bbar*log(1-q)
	R <- bsum+bbarsum

	bsum <- sum(1/(q+(ps[bstrains]/(1-ps[bstrains]))))
	bbarsum <- bbar/(1-q)
	dR <- bsum-bbarsum

	return(list(R=R,dR=dR))
}

findRdR <- function(tree,edge){
	#Find q at which dR equals zero UNUSED
	RdR <- calculateRdR(tree,edge,0)
	if(RdR$dR < 0){
		return(list(q=0,R=RdR$R,dR=RdR$dR,iter=0))
	}

	low <- 0
	high <- 1-10^-6
	i <- 0

	while(abs(RdR$dR)>1e-3 & i<100){
		i <- i+1
		q <- (high+low)/2
		RdR <- calculateRdR(tree,edge,q)
		if(RdR$dR<0){
			high <- q
		}else{
			low <- q
		}
	}

	return(list(q=q,R=RdR$R,dR=RdR$dR,iter=i))
}

plotTrees <- function(trees,live=F){
	#Plot a list of trees
	for(i in 1:length(trees)){
		if(length(trees[[i]]$tip.label)==1){
			plot(0,type="n",axes=F,xlab="",ylab="",main=trees[[i]]$tip.label,col.main=rainbow(length(trees))[i])
		}else{
			census <- treeCensus(trees[[i]])
			paraEdges <- findParaEdges(trees[[i]])
			paraTips <- findParaTips(trees[[i]])
			heat <- colorRampPalette(c("yellow","red"))
			plot(trees[[i]],type="u",lab4ut="axial",cex=0.5,tip.color=1+paraTips,edge.color=c("black",heat(max(paraEdges)),"green")[1+paraEdges+ifelse(paraEdges==census$plural/2 & paraEdges>0,1,0)])
			title(paste("No. diff strains:",census$present,"; No. repeat strains:",census$plural,";\nTip-to-tip:",ttt(trees[[i]])))
			edgelabels(paste(trees[[i]]$edge.length,paraEdges),cex=0.5,col=4,bg="white")
		}
		if(live){
			line <- readline()
		}
	}
}

plotCutTree <- function(tree,trees,use.edge.length=T,edgelabels=F,highlight=0){
	#Plots starting tree with colours according to resultant trees
	census <- treeCensus(tree)
	tipIndex <- rep(0,ntips(tree))
	names(tipIndex) <- tree$tip.label
	for(i in 1:length(trees)){
		tipIndex[trees[[i]]$tip.label] <- i
	}
	edgecolors <- rep(1,length(tree$edge.length))
	if(highlight>0){
		edgecolors[highlight]  <- 2
	}
	plot(tree,type="u",lab4ut="axial",cex=0.5,tip.color=rainbow(length(trees))[tipIndex[tree$tip.label]],edge.color=edgecolors,use.edge.length=use.edge.length,lwd=2)
	title(paste("No. diff strains:",census$present,"; No. repeat strains:",census$plural,"; Tip-to-tip:",ttt(tree)))
	if(edgelabels){
		edgelabels(tree$edge.length,cex=0.5,col=4,bg="white")
	}
}

outputGroups <- function(cutTrees){
	#Output ortholog groups based on a list of lists of trees
	letters <- c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z","aa","bb","cc","dd","ee","ff","gg","hh","ii","jj","kk","ll","mm","nn","oo","pp","qq","rr","ss","tt","uu","vv","ww","xx","yy","zz")
	for(i in 1:length(cutTrees)){
		for(j in 1:length(cutTrees[[i]])){
			cat("SC1",i+999,letters[j]," ",sep="")
			for(k in 1:length(cutTrees[[i]][[j]]$tip.label)){
				cat(cutTrees[[i]][[j]]$tip.label[k],"")
			}
			cat("\n")
		}
	}
}

if(0){
	#LOAD TREES
	aaTrees <- loadTrees(8955)$aaTrees
}
if(0){
	#CUT TREES, PLOT TREES AND OUTPUT NEW GROUPS
	plot=T
	cutTrees <- list()
	lengths <- list()
	ordinals <- list()
	presents <- list()
	plurals <- list()
	representations <- list()
	losses <- list()
	resolutions <- list()
	treenumbers <- list()
	for(i in 1:length(aaTrees)){
		pdf(paste("strains/clustering/pdf/SC1",i+999,"-cut.pdf",sep=""),width=30,height=10)
		if(i%%45==0){
			cat(".")
		}
		results <- removeParalogs(aaTrees[[i]],mode="long",plot=plot,cutoff=0.05)
		cutTrees[[i]] <- results$subTrees
		lengths[[i]] <- results$lengths
		ordinals[[i]] <- results$ordinals
		presents[[i]] <- results$presents
		plurals[[i]] <- results$plurals
		representations[[i]] <- results$representations
		losses[[i]] <- results$losses
		resolutions[[i]] <- results$resolutions
		treenumbers[[i]] <- rep(i,length(results$subTrees))
		if(plot){
			plotCutTree(aaTrees[[i]],results$subTrees,edgelabels=T)
			plotCutTree(aaTrees[[i]],results$subTrees,use.edge.length=F)
			plotTrees(aaTrees[i],live=F)
		}
		dev.off()
	}
	cat("\n")
	sink("strains/clustering/new-groups.txt")
	outputGroups(cutTrees)
	sink()

	finalTrees <- vector()
	for(i in 1:length(lengths)){
		finalTrees <- c(finalTrees,rep(i,length(lengths[[i]])))
	}
	finalLengths <- unlist(lengths)
	finalOrdinals <- unlist(ordinals)
	finalPresents <- unlist(presents)
	finalPlurals <- unlist(plurals)
	finalRepresentations <- unlist(representations)
	finalLosses <- unlist(losses)
	finalResolutions <- unlist(resolutions)

	finalDF <- as.data.frame(cbind(finalTrees,finalLengths,finalOrdinals,finalPresents,finalPlurals,finalRepresentations,finalLosses,finalResolutions))
	colnames(finalDF) <- c("Tree","Length","Ordinal","Present","Plural","Representation","Losses","Resolutions")

	finalCheck <- unlist(lapply(cutTrees,function(x) unlist(lapply(x,function(y) treeCensus(y)$plural))))
}
