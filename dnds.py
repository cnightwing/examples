#Script to calculate dn/ds statistics for a nucleotide alignment created from an amino acid alignment
import sys,re
import numpy as np
from Bio import AlignIO,SeqIO
from Bio.Seq import Seq,translate
from Bio.Alphabet import generic_dna

if len(sys.argv)<2:
    print 'Usage: python dnds.py <nt alignment file> [ungapped]'
    sys.exit(0)
aln = AlignIO.read(sys.argv[1],"clustal")
if len(sys.argv)>2:
	ungapped = True
else:
	ungapped = False

def catAln(alns):
        alphabet = alns[0][0,:].seq.alphabet
        catSeqs = list()
        for i in range(len(alns[0])):
                catSeq = list()
                for j in range(len(alns)):
                        catSeq.append(str(alns[j][i].seq))
                catSeq = ''.join(catSeq)
                catSeqs.append(catSeq)
        result = AlignIO.MultipleSeqAlignment(SeqIO.SeqRecord(Seq(catSeqs[x],alphabet=alphabet),id=alns[0][x,:].id) for x in range(0,len(catSeqs)))
        return(result)

def stripGaps(align):
	goodCols = [align[:,x:x+1] for x in range(0,align.get_alignment_length()) if "-" not in align[:,x]]
	return(catAln(goodCols))		

def nearCodons(codon):
	codon = codon.upper()
	ncs = list()
	for pos in range(3):
		for nt in 'ACGT':
			new = list(codon)
			new[pos] = nt
			ncs.append("".join(new))
	ncs = [nc for nc in ncs if nc!=codon]
	return(ncs)

def codonChange(codon):
	aa = translate(codon)
	ncs = nearCodons(codon)
	aas = [translate(nc) for nc in ncs]
	ns = [("N","S")[int(naa==aa)] for naa in aas]
	return(ns)

#Create NS map for all codons
allCodons = [x+y+z for x in 'ACGT' for y in 'ACGT' for z in 'ACGT']
NSDict = dict()
for codon in allCodons:
	NSDict[codon] = codonChange(codon)

#Remove gapped columns
completeAln = aln
if not ungapped:
	aln = stripGaps(aln)

codonList = list()
aaList = list()
codonCount = list()
snpCount = list()
codonStatus = list()
codonStatusMN = list()
codonStatusMS = list()
snpStatus = list()
nsList = list()
majorCount = list()
snpPosition = list()

for i in range(0,aln.get_alignment_length(),3):
	if i%1000==0:
		print >> sys.stderr, i
	
	codons = aln[:,i:i+3]
	codonSet = list(set([str(codon.seq) for codon in codons]))
	codonList.append(codonSet)
	aaList.append([translate(codon) for codon in codonSet])
	codonCount.append(len(codonSet))
	nsList.append(np.array([NSDict[str(codon.seq)] for codon in codons]))

	#Conserved codons
	if len(codonSet)==1:
		codonStatus.append("C")
		snpStatus.extend("CCC")

	#Multiple codons have NS status determined by all pairs within a single mutation of each other
	if len(codonSet)>2:
		#For each codon in the set, calculate the distance to each other codon and determine NS status if 1 snp
		changeSet = list()
		for j in range(0,len(codonSet)):
			aa = translate(codonSet[j])
			for k in range(j,len(codonSet)):
				if sum([int(codonSet[j][x]!=codonSet[k][x]) for x in range(0,3)])==1:
					if translate(codonSet[k])==aa:
						changeSet.append("S")
					else:
						changeSet.append("N")
		#Minimal number of changes is codonSet-1, so final count is scaled to sum to this value
		codonStatusMN.append(changeSet.count("N")*(float(len(codonSet)-1)/len(changeSet)))
		codonStatusMS.append(changeSet.count("S")*(float(len(codonSet)-1)/len(changeSet)))
                codonStatus.append("M")
		for j in range(0,3):
			snpStatus.append("C") if len(set(codons[:,j]))==1 else snpStatus.append("M")

	#Exactly two codons
	if len(codonSet)==2:
		majorCount.append(max([sum(str(codon.seq)==x for codon in codons) for x in codonSet]))
		diff = [codonSet[0][i]!=codonSet[1][i] for i in range(0,3)]
		snpCount.append(sum(diff))
		if(sum(diff)==1):
			snpPosition.extend([i+1 for i in range(0,3) if diff[i]])
		else:
			snpPosition.append(0)
		aaSet = set([translate(codon) for codon in codonSet])
		if len(aaSet)==2:
			codonStatus.append("N")
			for j in range(0,3):
				if diff[j] & (sum(diff)==1):
					snpStatus.append("N")
				elif diff[j]:
					snpStatus.append("A")
				else:
					snpStatus.append("C")
		else:
			codonStatus.append("S")
			for j in range(0,3):
				if diff[j] & (sum(diff)==1):
                                        snpStatus.append("S")
                                elif diff[j]:
                                        snpStatus.append("A")
                                else:
                                        snpStatus.append("C")
	else:
		snpCount.append(0)
		snpPosition.append(0)

bgN = sum([list(ns.flatten()).count("N") for ns in nsList])
bgS = sum([list(ns.flatten()).count("S") for ns in nsList])
fgN = codonStatus.count("N")
fgS = codonStatus.count("S")
fgMN = sum(codonStatusMN)
fgMS = sum(codonStatusMS)

posbgN = list()
posbgS = list()
posfgN = list()
posfgS = list()
for pos in [0,2]:
	posbgN.append(sum([list(ns[:,0+(pos*3):3+(pos*3)].flatten()).count("N") for ns in nsList]))
	posbgS.append(sum([list(ns[:,0+(pos*3):3+(pos*3)].flatten()).count("S") for ns in nsList]))
	fg = [codonStatus[x] for x in range(0,len(codonStatus)) if snpPosition[x]==pos+1]
	posfgN.append(fg.count("N"))
	posfgS.append(fg.count("S"))

ntAlleles = len(set([str(seq.seq) for seq in aln]))
aaAlleles = len(set([translate(seq) for seq in aln]))

#SeqCount ntAlleles aaAlleles AlnLength GoodCodons BiCodons Multicodons w1snp w2snp w3snp bgN bgS fgN fgS fgMN fgMS (pos1bgN pos1bgS pos1fgN pos1fgS pos3bgN pos3bgS pos3fgN pos3fgS)
if ungapped:
	print len(completeAln),ntAlleles,aaAlleles,completeAln.get_alignment_length(),len(codonList),codonCount.count(2),len(codonStatusMN),snpCount.count(1),snpCount.count(2),snpCount.count(3),bgN,bgS,fgN,fgS,fgMN,fgMS#,posbgN[0],posbgS[0],posfgN[0],posfgS[0],posbgN[1],posbgS[1],posfgN[1],posfgS[1]
#OGG SeqCount ntAlleles aaAlleles AlnLength GoodCodons Bicodons Multicodons w1snp w2snp w3snp bgN bgS fgN fgS fgMN fgMS M A N S (pos1bgN pos1bgS pos1fgN pos1fgS pos3bgN pos3bgS pos3fgN pos3fgS)
else:
	print re.findall("SC1\d+\w",sys.argv[1])[0],len(completeAln),ntAlleles,aaAlleles,completeAln.get_alignment_length(),len(codonList),codonCount.count(2),len(codonStatusMN),snpCount.count(1),snpCount.count(2),snpCount.count(3),bgN,bgS,fgN,fgS,fgMN,fgMS,snpStatus.count("M"),snpStatus.count("A"),snpStatus.count("N"),snpStatus.count("S")#,posbgN[0],posbgS[0],posfgN[0],posfgS[0],posbgN[1],posbgS[1],posfgN[1],posfgS[1]

if ungapped:
	#Print snp status
	fo = open("strains/dnds/snpStatus.txt","w")
	for n,snp in enumerate(snpStatus):
		fo.write(snp)
	fo.close()
	
	#Extract only snp columns
	snps = []
	poss = []
	for n,snp in enumerate(snpStatus):
		if snp!="C":
			snps.append(aln[:,n])
			poss.append(n)
		if not n%1000:
			print >> sys.stderr, n
	
#	fo = open("strains/dnds/onlySnps.txt","w")
#	for pos in poss:
#		fo.write(str(pos)+" ")
#	for s in range(0,len(aln)):
#		fo.write(aln[s].id+" ")
#		for n in range(0,len(snps)):
#			fo.write(snps[n][s]+" ")
#		fo.write("\n")
#	fo.close()
	

#print "Alignment Length:",len(codonList)*3
#print "Background N/S:",bgdNdS
#print "Foreground N/S:",fgdNdS
#print "Codon Count:",{x:codonCount.count(x) for x in set(codonCount)}
#print "SNP Count:",{x:snpCount.count(x) for x in set(snpCount)}
#print "SNP Position:",{x:snpPosition.count(x) for x in set(snpPosition)}
