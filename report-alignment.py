#Program to produce a report/statistics on a mummer alignment of two sets of contigs
#Requires coords and snps files
#Represent each strain as a list of contigs, each contig a list of letters
#Start with the longest alignments and don't overwrite once a base has an assignment

import sys,time
from Bio import SeqIO
from Bio.Seq import Seq,translate
from Bio.Alphabet import generic_dna
import pandas as pd
import numpy as np
from itertools import groupby
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages

pd.options.mode.chained_assignment = None

def vcount(v):
	return[[name,len(list(group))] for name,group in groupby(v)]

#Some fixed values for the graphics
hl = 2
ll = -2
minsize = 100000

#Get strains from command line argument so the program can be parallelised
if len(sys.argv)<4:
    print 'Usage: python report-alignment.py <strain 1> <strain 2> <plot 1/0>'
    sys.exit(0)

s1 = sys.argv[1]
s2 = sys.argv[2]
vis = sys.argv[3]

#Import contigs and make lists of contig coverage for later
contigs1 = list(SeqIO.parse("strains/alignment/arranged/"+s1+"-arranged.fasta",'fasta'))
contigs2 = list(SeqIO.parse("strains/alignment/arranged/"+s2+"-arranged.fasta",'fasta'))

cov1 = {contig.id:np.zeros(len(str(contig.seq))) for contig in contigs1}
cov2 = {contig.id:np.zeros(len(str(contig.seq))) for contig in contigs2}

#Import the gene table to label where the snps are
geneTable = pd.read_table("strains/clustering/new-groups/geneTable.txt",sep=" ")
genes1 = geneTable[geneTable['Strain']==s1]
genes2 = geneTable[geneTable['Strain']==s2]
genes1['Shared'] = [gene in genes2.Group.tolist() for gene in genes1.Group.tolist()]
genes2['Shared'] = [gene in genes1.Group.tolist() for gene in genes2.Group.tolist()]
snpnotes1 = {contig.id:[] for contig in contigs1}
snpnotes2 = {contig.id:[] for contig in contigs1}

#Import and sort coords file by alignment length
aligns = pd.read_table("strains/alignment/coords/"+s1+"-"+s2+".coords",header=-1,index_col=False)
aligns[9] = aligns[4]*aligns[6]
aligns.sort_values(by=9,axis=0,ascending=False,inplace=True)

#Import snps
snps = pd.read_table("strains/alignment/snps/"+s1+"-"+s2+".snps",header=None,index_col=False,sep='[\t ]',engine='python')

#Create a figure for each contig
if vis:
	voffsets = {contig.id:0 for contig in contigs1}
	for contig in [contig for contig in contigs1 if len(contig)>minsize]:
		print contig.id
		cgenes = genes1[genes1['Contig']==contig.id]
		plt.figure(contig.id,figsize=(len(contig)/2000,2))
		plt.xlim(0,len(contig))
		plt.ylim(-10,10)
		plt.title(contig.id)
		ax = plt.gca()
		ax.axis('on')
		plt.xticks(np.arange(0,len(contig),5000))
		plt.yticks([0])
		ax.xaxis.grid(True)
		ax.plot([0,len(contig)],[hl,hl],'b-')
		ax.plot([0,len(contig)],[ll,ll],'b--')
		for i,row in cgenes.iterrows():
			if isinstance(row.Name,float):
				row.Name = ""
       		 	ax.add_patch(patches.Arrow(row.Start,hl,row.End-row.Start,0,1))
		        ax.text((row.Start+row.End)/2,hl+0.5,row.Group+" "+str(row.Name),rotation=90,va='bottom',color=('b' if row.Shared else 'r'),size='xx-small')

#wtf = list()
#Loop through the alignments and mark coverage and snps
for n,align in aligns.iterrows():
	#if n>0:
	#	continue
	#sys.stderr.write(str(i)+" ")
	contig1 = align[7]
	contig2 = align[8]

	ctg1 = [contig for contig in contigs1 if contig.id==contig1][0]

	if align[0] < align[1]:
		region1 = np.arange(align[0]-1,align[1])
	else:
		region1 = np.arange(align[0]-1,align[1]-2,-1)
	if align[2] < align[3]:
		dir = "fwd"
	        region2 = np.arange(align[2]-1,align[3])
	else:
		dir = "rev"
		region2 = np.arange(align[2]-1,align[3]-2,-1)

        index1 = set([i for i,x in enumerate(region1) if cov1[contig1][x]==0])
        index2 = set([i for i,x in enumerate(region2) if cov2[contig2][x]==0])
	index = list(index1.intersection(index2))
	region1 = region1[index]
	region2 = region2[index]

	if (dir=="fwd") & (len(region1)>0) &(len(region2)>0):
		offset = min(region1)-min(region2)
	elif (dir=="rev") & (len(region1)>0) &(len(region2)>0):
		offset = min(region1)-(len(ctg1)-max(region2))
	else:
		offset = 0

	cov1[contig1][region1] = 1
	cov2[contig2][region2] = 1

	#Plotting the aligned regions
	if len(ctg1)>minsize:
		cgenes = genes2[(genes2['Contig']==contig2) & (([start in region2 for start in genes2['Start']]) or ([end in region2 for end in genes2['End']]))]
		if vis:
			plt.figure(contig1)
			print contig1
			voffset = voffsets[contig1]
			ax = plt.gca()
			ax.plot(align[0:2],[ll-voffset,ll-voffset],'b-')
			ax.plot(min(align[0:2]),ll-voffset,'ko')
			ax.plot(max(align[0:2]),ll-voffset,'ko')
			#voffsets[contig1] = voffsets[contig1]+0.25
			if dir=="fwd":
				for i,row in cgenes.iterrows():
					if isinstance(row.Name,float):
						row.Name = ""
					ax.add_patch(patches.Arrow(offset+row.Start,ll,row.End-row.Start,0,1))
		       	 	        ax.text(offset+(row.Start+row.End)/2,ll-0.5,str(row.Name)+" "+row.Group,rotation=90,va='top',color=('b' if row.Shared else 'r'),size='xx-small')
			if dir=="rev":
				ctglen = plt.xlim()[1]
				for i,row in cgenes.iterrows():
					if isinstance(row.Name,float):
						row.Name = ""
					ax.add_patch(patches.Arrow(offset+ctglen-row.Start,ll,row.Start-row.End,0,1))
					ax.text(offset+ctglen-(row.Start+row.End)/2,ll-0.5,str(row.Name)+" "+row.Group,rotation=90,va='top',color=('b' if row.Shared else 'r'),size='xx-small')
			ax.text((align[0]+align[1])/2,-16,contig2,va='bottom',ha='center',size='xx-small')
	asnps = snps[snps[0]==n+1]
	for m,snp in asnps.iterrows():
		if (snp[1] in region1) & (snp[4] in region2):
			if (snp[2]==".") | (snp[3]=="."):
				cov1[contig1][snp[1]] = -2
				cov2[contig2][snp[4]] = -2
				if (len(ctg1)>minsize) & bool(vis):
					ax.plot([snp[1],snp[1]],[ll+0.25,ll-0.25],'b-',linewidth=0.1)
			else:
				cov1[contig1][snp[1]] = -1
				cov2[contig2][snp[4]] = -1
				if (len(ctg1)>minsize) & bool(vis):
					ax.plot([snp[1],snp[1]],[ll+0.25,ll-0.25],'r-',linewidth=0.1)
			cgenes1 = genes1[genes1['Contig']==contig1]
			gene1 = cgenes1[((cgenes1['Start']>=snp[1]) & (cgenes1['End']<=snp[1])) | ((cgenes1['Start']<=snp[1]) & (cgenes1['End']>=snp[1]))]
			if len(gene1):
				snpnotes1[contig1].append(gene1.index.tolist()[0])
			else:
				snpnotes1[contig1].append("None")
			cgenes2 = genes2[genes2['Contig']==contig2]
			gene2 = cgenes2[((cgenes2['Start']>=snp[4]) & (cgenes2['End']<=snp[4])) | ((cgenes2['Start']<=snp[4]) & (cgenes2['End']>=snp[4]))]
			if len(gene2):
				snpnotes2[contig1].append(gene2.index.tolist()[0])
			else:
				snpnotes2[contig1].append("None")

        #sys.stderr.write(str(abs(align[0]-align[1]))+" "+str(len(region1))+" "+str(list(cov1[contig1][region1]).count(-1))+"\n")
        #sys.stderr.write(str(abs(align[2]-align[3]))+" "+str(len(region2))+" "+str(list(cov2[contig2][region2]).count(-1))+"\n\n")
	#wtf.append(cov1[contigs1[0].id])

if vis:
	for contig in [contig for contig in contigs1 if len(contig)>minsize]:
		id = contig.id
		snppos = [x for x in np.arange(0,len(contig)) if cov1[id][x]==-1]
		hist = np.histogram(snppos,np.arange(0,1+(len(contig)/1000))*1000)
		if np.sum(hist[0])>0:
			plt.figure(id)
			ax=plt.gca()
			ax.plot(hist[1][1:]-500,[(float(x)/max(hist[0])*2)-1 for x in hist[0]],'r-')

#Summary statistics
stats = open("strains/alignment/stats/"+s1+"-"+s2+".stats",'w')
stats.write(s1+" ")
stats.write(s2+" ")
stats.write(str(sum([len(c) for c in contigs1]))+" ")
stats.write(str(sum([len(c) for c in contigs2]))+" ")
stats.write(str([x for c in cov1.values() for x in c].count(1))+" ")
stats.write(str([x for c in cov2.values() for x in c].count(1))+" ")
stats.write(str([x for c in cov1.values() for x in c].count(0))+" ")
stats.write(str([x for c in cov2.values() for x in c].count(0))+" ")
stats.write(str([x for c in cov1.values() for x in c].count(-1))+" ")
stats.write(str([x for c in cov2.values() for x in c].count(-1))+" ")
stats.write(str([x for c in cov1.values() for x in c].count(-2))+" ")
stats.write(str([x for c in cov2.values() for x in c].count(-2))+" ")
stats.close()

rep1 = open("strains/alignment/coverage/"+s1+"-"+s2+".coverage",'w')
for contig in contigs1:
	cov = (vcount(cov1[contig.id]))
	i = 0
	for part in cov:
		rep1.write(contig.id+" "+str(len(contig.seq))+" "+str(int(part[0]))+" "+str(part[1]))
		if part[0]==-1:
			if snpnotes1[contig.id][i]=="None":
				rep1.write(" NoGene NoDesc")
			else:
				rep1.write(" "+str(geneTable.iloc[snpnotes1[contig.id][i],6])+" "+str(geneTable.iloc[snpnotes1[contig.id][i],7]))
			if snpnotes2[contig.id][i]=="None":
				rep1.write(" NoGene NoDesc")
			else:
				rep1.write(" "+str(geneTable.iloc[snpnotes2[contig.id][i],6])+" "+str(geneTable.iloc[snpnotes2[contig.id][i],7]))
			i+=1
		rep1.write("\n")
rep1.close()

rep2 = open("strains/alignment/blocks/"+s1+"-"+s2+".blocks",'w')
cov1flat = np.concatenate(list(cov1.values()))
blocks1 = np.zeros(len(cov1flat/100))
for i in range(0,len(blocks1)):
	blocks1[i] = list(cov1flat[(i*100):(i*100)+1000]).count(-1)
	rep2.write(str(blocks1[i])+"\n")
rep2.close()

rep3 = open("strains/alignment/gaps/"+s1+"-"+s2+".gaps",'w')
allparts = []
for contig in contigs1:
	cov = cov1[contig.id]
	bincov = [x if x==0 else 1 for x in cov]
	covparts = vcount(bincov)
	allparts.append(covparts)
	for part in covparts:
		rep3.write(contig.id+" "+str(len(contig.seq))+" "+str(int(part[0]))+" "+str(part[1])+"\n")
rep3.close()

rep4 = open("strains/alignment/gapsums/"+s1+"-"+s2+".gaps",'w')
allparts2 = []
for parts in allparts:
	if len(parts)>1:
		if parts[0][0]==1:
			parts.pop(0)
		if parts[-1][0]==1:
			parts.pop(-1)
	if len(parts)>0:
		allparts2.append(parts)
allparts2 = np.concatenate(allparts2)
gaps = [x[1] for x in allparts2 if x[0]==0]
fills = [x[1] for x in allparts2 if x[0]==1]
rep4.write(s1+" "+s2+" "+str(len(gaps))+" "+str(np.mean(gaps))+" "+str(np.std(gaps))+" "+str(len(fills))+" "+str(np.mean(fills))+" "+str(np.std(fills))+"\n")
rep4.close()

if vis:
	pp = PdfPages("strains/alignment/pdf/"+s1+"-"+s2+".pdf")
	for f in plt.get_fignums():
		plt.figure(f)
		pp.savefig(bbox_inches='tight')
	pp.close()

#Try to plot some of this madness
#contig = contigs1[0]
#cgenes1 = genes1[genes1['Contig']==contig.id]
#fig = plt.figure(figsize=(len(contig)/1000,2))
#plt.xlim(0,len(contig))
#plt.ylim(-2,2)
#ax = plt.gca()
#for i,row in cgenes1.iterrows():
#	ax.add_patch(patches.Arrow(row.Start,0.5,row.End-row.Start,0,1))
#	ax.text((row.Start+row.End)/2,0.5,row.Group+" "+str(row.Name),rotation=90,va='bottom')
#pp = PdfPages('pytest.pdf')
#pp.savefig()
##pp.close()
