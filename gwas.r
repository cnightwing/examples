library(plyr)
source("~/phenotyping/phenotype-functions.r")

strains <- unlist(read.table("~/strains/sc1_strains.txt",stringsAsFactors=F))
clades <- read.table("~/strains/clades.txt",row.names=1)

entropy <- function(a){
#Function to calculate the entropy of a vector
	at = table(a)
	at = at/sum(at)
	H = -sum(at*log2(at))
	return(H)
}

mi <- function(a,b){
#Function to calculate the mutual information between two vectors, unless one of them does not vary (NAs are ignored)
	valid = !is.na(a) & !is.na(b)
	at = table(a[valid])
	at = at/sum(at)
	bt = table(b[valid])
	bt = bt/sum(bt)
	if(length(at)<2|length(bt)<2){
		return(0)
	}
	ct = table(a[valid],b[valid])
	ct = ct/sum(ct)
	
	I = 0
	for(j in 1:length(at)){
		for(k in 1:length(bt)){
			if(ct[j,k]>0){
				i = ct[j,k] * log2(ct[j,k]/(at[j]*bt[k]))
				I = I+i
			}
		}
	}
	return(as.vector(I))
}

rmi <- function(a,b){
#Function to calculate the mutual information between two vectors as a fraction of the maximum possible information
	I <- mi(a,b)
	minH <- min(entropy(a[!is.na(b)]),entropy(b[!is.na(a)]))
	if(minH==0){
		return(0)
	}
	return(I/minH)
}

wmi <- function(a,b){
#Function to calculate the mutual information between a discrete vector and a vector of probabilities, as a fraction of the maximum possible information
	valid = !is.na(a) & !is.na(b)
	a = a[valid]
	b = b[valid]

	at = table(a)
	at = at/sum(at)

	bt = c(sum(b),sum(1-b))
	bt = bt/sum(bt)

        if(length(at)<2|bt[1]==0|bt[2]==0){
                return(0)
        }

        I = 0
        for(aind in 1:length(at)){
		#b+
		aval = names(at)[aind]
        	bval = sum(b[a==aval])/length(b)
		if(bval>0){
			i = bval * log2(bval/(at[aind]*bt[1]))
			I = I+i
		}
		#b-
		aval = names(at)[aind]
		bval = sum((1-b)[a==aval])/length(b)
		if(bval>0){
			i = bval * log2(bval/(at[aind]*bt[2]))
			I = I+i
		}
	}

	minH <- min(-sum(at*log2(at)),-sum(bt*log2(bt)))
        return(as.vector(I/minH))
}

identifyStrains <- function(dm){
        strains <- unlist(regmatches(rownames(dm),gregexpr("[A-Z]\\d\\d?",rownames(dm))))
        return(strains)
}

countAlleles <- function(dm){
        alleles <- nrow(unique(dm))
        return(alleles)
}

whichAlleles <- function(dm){
	dmStrains <- identifyStrains(dm)
        rows <- apply(dm,1,v2str)
        alleles <- apply(unique(dm),1,v2str)
	results <- match(rows,alleles)
	names(results) <- dmStrains
        return(results)
}

assignAlleles <- function(dm){
	dmAlleles <- whichAlleles(dm)
	alleles <- dmAlleles[strains]
	alleles[is.na(alleles)] <- 0
	names(alleles) <- strains
	return(alleles)	
}

#Four datasets: phenotypes, genotypes, alleles, snps
if(import){
cat("Importing data\n")
grandTable <- read.table("~/phenotyping/grandTable.txt")
phenoGood <- apply(grandTable,1,function(x) as.logical(length(table(x))-1))

geneMatrix <- read.table("~/strains/clustering/new-groups/geneMatrix.txt",header=T,row.names=1)
prabMatrix <- geneMatrix
prabMatrix[prabMatrix>1] <- 1
geneGood <- apply(prabMatrix,1,function(x) sum(x)<92)
aaGood <- apply(geneMatrix,1,function(x) as.logical(sum(x>1)==0 & sum(x)>1))

aaDMs <- list()
for(i in rownames(prabMatrix)[aaGood]){
        if(!is.na(file.info(paste("~/strains/post-clustering/aadm/",i,".dm",sep=""))$size)){
                aaDMs[[length(aaDMs)+1]] <- read.table(paste("~/strains/post-clustering/aadm/",i,".dm",sep=""))
        }
}
aaAlleles <- lapply(aaDMs,assignAlleles)
aaAlleleMatrix <- do.call(rbind,aaAlleles)
rownames(aaAlleleMatrix) <- rownames(prabMatrix)[aaGood]
colnames(aaAlleleMatrix) <- strains

snpMatrix <- read.table("snpmatrix.txt",header=T)
names(snpMatrix) <- sapply(names(snpMatrix),function(x) sub(".scaffolds","",x))
snpMatrix <- snpMatrix[,strains]

redSnpMatrix <- count(snpMatrix)[,1:92]

#snpCount <- read.table("~/strains/dnds/reducedSnps.txt",colClasses=c("character","integer"))
#snpMatrix <- t(apply(snpCount,1,function(x) str2v(x[1])))
#colnames(snpMatrix) <- sort(colnames(prabMatrix))
#snpMatrix <- t(apply(snpMatrix,1,function(x) x[colnames(prabMatrix)]))
}

if(geneCalc){
#Phenotype-Genotype MI
cat("Calculating MI between phenotypes and gene presence/absence\n")
PhenoGeneMI <- apply(grandTable[phenoGood,],1,function(x) apply(prabMatrix[geneGood,],1,function(y) rmi(x,y)))
}
if(alleleCalc){
#Phenotype-Allele MI
cat("Calculating MI between phenotypes and alleles\n")
PhenoAlleleMI <- apply(grandTable[phenoGood,],1,function(x) apply(aaAlleleMatrix,1,function(y) rmi(x,y)))
}
if(snpCalc){
#Phenotype-Snp MI
cat("Calculating MI between phenotypes and snps\n")
PhenoSnpMI <- apply(grandTable[phenoGood,],1,function(x) apply(snpMatrix,1,function(y) rmi(x,y)))
}

geneTable <- read.table("~/strains/clustering/new-groups/geneTable.txt",header=T,row.names=1)

#Function to backreference an OGG
oggInfo <- function(ogg){
	info <- geneTable[geneTable$Group==ogg,]
	if(nrow(info)==0){
		info <- geneTable[ogg,]
	}
	strains <- paste(sort(unique(info$Strain)),collapse=",")
	names <- paste(unique(info$Name),collapse=" | ")
	descs <- paste(unique(info$Description),collapse=" | ")
	
	return(data.frame(OGG=ogg,Strains=strains,Names=names,Descriptions=descs,stringsAsFactors=F))
}

#Function to find the gene corresponding to a snp
coreTable <- read.table("~/strains/core-genomes/codingCoreStarts.txt",header=F,row.names=1)
coreStarts <- coreTable[,3]
coreEnds <- c(coreTable[2:nrow(coreTable),3],coreTable[nrow(coreTable),3]+10000)
geneInfo <- function(pos){
	return(coreTable[which(pos>coreStarts & pos<coreEnds),])
}

#Function to backreference a Snp
snpIndex <- read.table("~/strains/dnds/reducedSnpsIndex.txt",header=F,colClasses=c("integer","character"))

snpInfo <- function(pattern){
	#Rearrange pattern back to alphabetical strain order (patterns from snp analysis are always in alphabetical order, whereas genes and phenotypes are in alphanumerical order)
	pattern = str2v(pattern)
	pattern = pattern[order(strains)]
	pattern = v2str(pattern)

	snps = which(snpIndex[,2]==pattern)
	poss = snpIndex[snps,1]
	genes = do.call("rbind",lapply(poss,geneInfo))
	return(genes)
}

#Produce reports on the best hits for each phenotype, and comparisons between the three tables



if(hits){
#Get the hits that have the best hits from the three tables
CombinedMI <- rbind(PhenoGeneMI,PhenoAlleleMI,PhenoSnpMI)
CombinedMI[is.na(CombinedMI)] <- 0

bestFive <- function(PhenoMI){
	PhenoMI <- round(PhenoMI,8)
	scores <- sort(unique(PhenoMI),decreasing=T)[1:5]
	hits <- vector()
	ids <- vector()
	for(i in 1:5){
		if(length(hits)<5){
			hits <- c(hits,PhenoMI[PhenoMI==scores[i]])
			ids <- c(ids,which(PhenoMI==scores[i]))
		}
	}

	types <- c(rep("Gene",nrow(PhenoGeneMI)),rep("Allele",nrow(PhenoAlleleMI)),rep("Snp",nrow(PhenoSnpMI)))[ids]
	bestScores <- vector()
	bestTypes <- vector()
	bestOggs <- vector()
	for(i in 1:length(types)){
		if(types[i]=="Snp"){
			oggs <- rownames(snpInfo(v2str(snpMatrix[ids[i]-(nrow(PhenoGeneMI)+nrow(PhenoAlleleMI)),])))
			bestScores <- c(bestScores,rep(hits[i],length(oggs)))
			bestTypes <- c(bestTypes,rep(types[i],length(oggs)))
			bestOggs <- c(bestOggs,oggs)
		}else{
			bestScores <- c(bestScores,hits[i])
			bestTypes <- c(bestTypes,types[i])
			bestOggs <- c(bestOggs,names(ids[i]))
		}
	}

	results <- cbind(bestScores,bestTypes,bestOggs)
	rownames(results) <- 1:nrow(results)
	return(results)
}

MIResults <- list()
for(i in 1:ncol(CombinedMI)){
	results <- bestFive(CombinedMI[,i])
	results <- cbind(rep(colnames(CombinedMI)[i],nrow(results)),results)
	MIResults[[length(MIResults)+1]] <- results
}
MIResults <- do.call("rbind",MIResults)
MIResults <- cbind(MIResults,sapply(MIResults[,4],function(x) oggInfo(x)[4]))
}

#Produce some graphs to show how useless this whole exercise was
library(Hmisc)
pdf("~/paper-phenotypes/gwas.pdf")
plot(as.vector(apply(PhenoGeneMI,2,max)),as.vector(apply(PhenoSnpMI,2,max)),xlab="OGG Mutual Information",ylab="Snp Mutual Information",pch=20,col=rgb(1,0,0,0.2),xlim=c(0,1),ylim=c(0,1))
abline(0,1)
plot(as.vector(apply(PhenoGeneMI,2,max)),as.vector(apply(PhenoSnpMI,2,max)),xlab="OGG Mutual Information",ylab="Snp Mutual Information",pch=20,col=rgb(1,0,0,0.2),xlim=c(0,1),ylim=c(0,1))
abline(0,1)
dev.off()

if(cladeCalc){
#Gwas by clade
cat("Looking at phenotypes within clades..\n")
polyClades <- which(table(clades)>1)
cladeGenoMI <- list()
cladeGenoWMI <- list()
cladeSNPMI <- list()
cladeAlleleMI <- list()
genoMI <- list()
snpMI <- list()
alleleMI <- list()
for(clade in polyClades){
	cladeStrains <- rownames(clades)[which(clades==clade)]

	cladePheno <- phenoMatrix[,cladeStrains]
	cladePhenoAna <- phenoMatrixAna[,cladeStrains]
	cladeFilter <- apply(cladePheno,1,function(x) length(table(x))>1)
	cladePheno <- cladePheno[cladeFilter,]
	cladePhenoAna <- cladePhenoAna[cladeFilter,]

	pdf(paste("clade",clade,"phenotypes.pdf",sep=""),width=10,height=15)
	for(pheno in rownames(cladePheno)){
		compareClade(pheno,strains[clades==clade])
	}
	dev.off()

	cladeGeno <- prabMatrix[,cladeStrains]
	cladeGeno <- cladeGeno[apply(cladeGeno,1,function(x) length(table(x))>1),]

	cladeSNP <- snpMatrix[,cladeStrains]
	cladeSNP <- as.data.frame(cladeSNP)[apply(cladeSNP,1,function(x) length(table(x))>1),]

	cladeAllele <- aaAlleleMatrix[,cladeStrains]
	cladeAllele <- cladeAllele[apply(cladeAllele,1,function(x) length(table(x))>1),]

	cladeGenoMI[[clade]] <- apply(cladePheno,1,function(x) apply(cladeGeno,1,function(y) rmi(x,y)))
	cladeGenoWMI[[clade]] <- apply(cladePhenoAna,1,function(x) apply(cladeGeno,1,function(y) wmi(y,x)))
	cladeSNPMI[[clade]] <- apply(cladePheno,1,function(x) apply(cladeSNP,1,function(y) rmi(x,y)))
	cladeAlleleMI[[clade]] <- apply(cladePheno,1,function(x) apply(cladeAllele,1,function(y) rmi(x,y)))

	genoMI[[clade]] <- PhenoGeneMI[rownames(cladeGeno),rownames(cladePheno)]
	snpMI[[clade]] <- PhenoSnpMI[as.numeric(rownames(cladeSNP)),rownames(cladePheno)]
	alleleMI[[clade]] <- PhenoAlleleMI[rownames(cladeAllele),rownames(cladePheno)]
}

save.image()
}


