# plant-mito-sharing
Simulation of molecular content and sharing in plant mitochondria


Say we have NM mitochondria in the cell, and NG mtDNA-encoded genes.

An individual mitochondrion i is described by:
* DNA_ij, how many copies of DNA encoding gene j it contains;
* mRNA_ij, how many transcripts for gene j it contains;
* protein_ij, how many copies of the protein product of gene j it contains.

Over time, the following processes occur within a mitochondrion i:
* Transcription, 		 mRNA_ij = mRNA_ij + 1, for every DNA_ij 
* Translation, 		 protein_ij = protein_ij + b, for every mRNA_ij
* RNA degradation,	 mRNA_ij = mRNA_ij - 1, for every mRNA_ij

The proteins within a mitochondrion are assigned to complexes as follows. The number of complete complexes is given by the minimum protein copy number across the different genes. All protein copy numbers in excess of this are assumed to be "lone" proteins, not yet embedded in a complex.
e.g. for protein counts 2, 5, 3, 10 we say there are 2 intact complexes and the remaining 0, 3, 1, 8 are unbound

Then we also have these processes in a mitochondrion i:
* Lone protein degradation,		 protein_ij -> protein_ij - 1, for every lone protein_ij
* Protein degradation within complex,	 protein_ij -> protein_ij - 1, for every embedded protein_ij

We also have interactions between mitochondria k and l:
* DNA exchange, 		  DNA_ki = DNA_li and DNA_li = DNA_ki
* RNA mixing, 		  mRNA_ki = r (mRNA_ki+mRNA_li) and mRNA_li = (1-r) (mRNA_ki+mRNA_li)
* Protein mixing,		  protein_ki = r(protein_ki+protein_li) and protein_li = (1-r)(protein_ki+protein_li)

We constrain DNA_ij to be either 0 or 1.

Shortcomings:
* translational, but not transcriptional, bursts are modelled (size b)
* uniform distributions for RNA and protein mixing feel wrong. binomial instead?
* proteins are immediately embedded in a complex if the stoichiometry is right. should be a delay?
* restrictive DNA structure -- everything encodes only one gene -- but this could be viewed as gene "regions" instead
* protein mixing on encounters doesn't respect unbound / bound structure

Here are some toy parameters:

* Transcription rate = 2 / hour
* Translation rate = 10 / hour
* Burst size b = 1
* RNA degradation rate = 0.25 / hour
* Lone protein degradation rate = 1 / hour
* Bound protein degradation rate = 0.01 / hour
* Mitochondrial encounter rate = (free to vary)
* DNA exchange probability per encounter = 0.
* RNA mixing probability per encounter = 1
* Protein mixing probability per encounter = 1
