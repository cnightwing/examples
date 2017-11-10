#Program to parse the coverage.txt output from realphy with the reference genbank file

import os,sys,time
from Bio import SeqIO
from Bio.Seq import Seq,translate
from Bio.Alphabet import generic_dna
import pandas as pd
import numpy as np
import gc

if len(sys.argv)<3:
	print 'Usage: python realphyCoreParse.py <Coverage File> <Genbank File>'
	sys.exit(0)

covFile = sys.argv[1]
genFile = sys.argv[2]
covDir = os.path.dirname(covFile)

index = pd.read_table(covFile,header=None,names=["contig","coverage"])
genome = list(SeqIO.parse(genFile,'genbank'))
genLen = sum([len(x.seq) for x in genome])

if len(index)!=genLen:
	exit("Incompatible Files")

#index.loc[:,'coding'] = False
#index.loc[:,'position'] = 0
#index.loc[:,'feature'] = "None"
#index.loc[:,'name'] = "None"

#for contig in genome:
#	for feature in [x for x in contig.features if (x.type!='source') & (x.type!='gene')][1:10]:
#		for part in feature.location.parts:
#			if feature.type=='CDS':
#				index.ix[part.start:part.end-1,'coding'] = True
#				index.ix[part.start:part.end-1,'position'] = [[2-feature.strand,2,2+feature.strand][(x-part.start)%3] for x in range(part.start,part.end)]
#				index.ix[part.start:part.end-1,'name'] = feature.qualifiers['gene']
#			index.ix[[x for x in range(part.start,part.end) if index.ix[x,'feature']=='None'],'feature'] = feature.type
#		gc.collect()

#index = index.ix[index['coverage']==1,:]
#index.to_csv(covDir+"/coreIndex.txt")

#List version runs CRAZY faster. Pandas sucks.
location = list(index.contig)
coverage = list(index.coverage)
coding = np.repeat(0,genLen)
position = np.repeat(0,genLen)
feat = np.repeat('None',genLen).astype(np.dtype('S32'),copy=False)
name = np.repeat('None',genLen).astype(np.dtype('S32'),copy=False)

for contig in genome:
	for feature in [x for x in contig.features if (x.type!='source') & (x.type!='gene')]:
		for part in feature.location.parts:
			if feature.type=='CDS':
				coding[part.start:part.end] = 1
				position[part.start:part.end] = [[2-feature.strand,2,2+feature.strand][(x-part.start)%3] for x in range(part.start,part.end)]
				name[part.start:part.end] = feature.qualifiers['gene']
			feat[[x for x in range(part.start,part.end) if feat[x]=='None']] = feature.type

index = pd.DataFrame({'contig':location,'coverage':coverage,'coding':coding,'position':position,'feature':feat,'name':name},columns=['contig','coverage','coding','position','feature','name'])
index = index.ix[index['coverage']==1,:]
index.to_csv(covDir+"/coreIndex.txt")
