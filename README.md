# examples

PYTHON

realphyCoreParser.py : realphy produces a core genome from multiple genomes or readsets, this program parses which bases of a given reference genome are included, and which genes they are part of if in a coding region
saga.py : interfaces with NCBI Blast to search their entire nucleotide database for hits within E. coli and other bacteria to find genes that likely come from non-E. coli species, in serial to conform to Blast usage rules
dnds.py : calculates the dN and dS scores for orthologous gene groups without assuming a particular phylogeny for those sequences
report-alignment.py : parses a pairwise alignment produced by nucmer (part of MUMmer) and produces a contigs-on-contigs plot of the sub-alignments complete with gene and snp labelling, primarily to visualise alignments and look for recombination events

R

synteny-ld.r : takes the gene annotations for multiple strains and calculates synteny across all of them based on orthologous gene groups, then estimates how this decays with distance, ie: as you look at further apart pairs of genes, is synteny maintained?
plot-alignment.r : similar to report-alignment.py, but produces an abstract view of the pairwise alignments by turning them into graph structures, where nodes represent genes and edges synteny.
analyse-trees.r : the result of lots of different attempts to decluster orthologous gene groups produced by orthoMCL, where it was clear that some groups contained too many paralogs, it includes functions for rearranging and cutting trees according to various criteria, labelling and plotting them to make it visually clear what the declustering process is doing, and reporting statistics on each tree and subtree.
