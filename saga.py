#SAGA: Search for Alien Genes Algorithm
#Blast a protein sequence against the NCBI database
#First against Escherichia only
#Second against all other bacteria
#Look for hits that score better in an organism other than Escherichia

import sys,os
import pandas as pd
from Bio.Blast import NCBIWWW as blast
from Bio.Blast import NCBIXML as parse
from Bio import SeqIO as seqio

#Input file
if len(sys.argv)<2:
    print('Usage: python <fasta file>')
    sys.exit(0)

seqs = list(seqio.parse(sys.argv[1],'fasta'))
results = pd.read_table("strains/post-clustering/saga-results-nt.txt",header=-1,delimiter=" ")
completed = list(results[0])
seqs = [seq for seq in seqs if seq.id not in completed]

escores = []
enames = []
oscores = []
onames = []
for seq in seqs:
	sys.stdout.write("\""+seq.name+"\" ")

	escher = blast.qblast('blastn','nr',seq.seq,entrez_query='txid561[Organism]')
	other = blast.qblast('blastn','nr',seq.seq,entrez_query='NOT txid561[Organism]') #txid2[Organism]

	eresults = parse.parse(escher)
	oresults = parse.parse(other)

	eresult = next(eresults)
	oresult = next(oresults)

	if(len(eresult.alignments)>0):
		ebest = eresult.alignments[0]
		sys.stdout.write(str(len(seq.seq))+" \""+ebest.title+"\" "+str(ebest.hsps[0].align_length)+" "+str(ebest.hsps[0].identities))
	else:
		sys.stdout.write(str(len(seq.seq))+" \"none\" "+"NA NA")
	if(len(oresult.alignments)>0):
		obest = oresult.alignments[0]
		sys.stdout.write(" \""+obest.title+"\" "+str(obest.hsps[0].align_length)+" "+str(obest.hsps[0].identities))
	else:
		sys.stdout.write(" \"none\" "+"NA NA")

	sys.stdout.write("\n")
	sys.stdout.flush()
