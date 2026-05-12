# plant-mito-sharing
Simulation of molecular sharing and mtDNA repair in plant mitochondria

<img width="2426" height="1310" alt="image" src="https://github.com/user-attachments/assets/d67966dd-9213-4a80-8bdb-1275a7ee183a" />

MtDNA repair requires the action of sets of protein machinery, which are often present at small copy numbers per mitochondrion. How do we ensure that all mtDNAs have access to the repair machinery required?

Simulation code is `damage-loop.c`, which takes a command-line parameter determining the specific experiment to run. These are passed by `damage-loop.sh`. `sim-plot-damage-both.R` plots the results and some probability calculations related to the problem.

Notes
---

Broadly, our model consists of a collection of mitochondria, which can contain different numbers of different proteins (which turn over with time) and mtDNA (which becomes randomly damaged over time). MtDNA damage can be repaired by the action of proteins. Different targetting, fusion, and exchange rules – the target of research in the study – respectively determine how proteins are imported to mitochondria, which mitochondria undergo fusion, and which biomolecules are exchanged upon fusion-fission events.
 
Model basics
A cell in our model (Fig. 1A) is pictured as containing n_M mitochondria. We consider a protein complex (which will act as a metabolon) that consists of n subunits. The copy number of each of these subunits, and mtDNA molecules, in mitochondrion i is described by a vector c^((i)) (where elements 1 through n describe copy numbers of each protein) and a scalar describing the mitochondrion’s copy number of mtDNA. Initially, n_M/3 mtDNA molecules (Fuchs et al., 2020; Preuten et al., 2010) are randomly assigned across the mitochondria in the cell. 
MtDNA may be intact or damaged. Intact mtDNA becomes damaged randomly with rate μ. Once damaged, mtDNA can be returned to the intact state through the action of repair protein complexes.

Repair protein dynamics
The cell produces proteins at rate λ. Each protein is imported by a mitochondrion with rate ρ according to a targetting rule (Fig. 1Ai). These targetting rules are: Damaged (only target to mitochondria containing damaged mtDNA); DNA-bearing (only target to mitochondria containing mtDNA); DNA-bearing, incomplete (only target to mitochondria with mtDNA but without a full set of repair proteins); Incomplete (only target to mitochondria without a combination of mtDNA and a full set of repair complexes); or Random (target each protein to a random mitochondrion). Proteins turn over with rate ν, reflecting their damage, removal, and replacement.
In our model, repair protein complexes can convert damaged mtDNA to intact mtDNA through two pathways. The first, base excision repair (BER)-like pathway, involves n_rep steps, each accomplished by a different protein complex acting on a single damaged mtDNA. The second, template repair (TR)-like pathway, involves one step, accomplished by a protein complex acting on a damaged mtDNA in the same compartment as another mtDNA (Fig. 1A). 
BER-like: 
Damaged –[complex 1]-> Repair state a –[complex 2]-> Repair step b –[…]-> Intact
TR-like:
Damaged + template –[complex 1]-> Intact + template
Although each step is generally catalysed by a protein complex, we consider only the lowest-copy-number component of each complex as the stoichiometric limiting factor in these processes, assuming that other subunits are present in excess.

Mitochondrial fusion and exchange dynamics
Mitochondrial fusion occurs with rate κ. Two rules determine behaviour upon fusion. The social rule determines which mitochondria fuse – this is either All mitochondria, or DNA-bearing (only mitochondria contain mtDNA fuse). The exchange rule (Fig. 1Aii) determines which biomolecules are mixed upon fusion. This can be Metabolons (complete sets of repair proteins are mixed), DNA + Subunits (individual mtDNA molecules and protein subunits are mixed); DNA-metabolons (complexes of mtDNA and complete sets of repair proteins are mixed); DNA-metabolons + Metabolons (DNA-metabolon complexes and metabolons without mtDNA are mixed); None (no mixing); Subunits (individual protein subunits are mixed). After mixing, fused mitochondria immediately fission, modelling “kiss-and-run” behaviour (Arimura, 2018; Logan, 2010).

Parameter values
Experimental evidence guides choices of some parameters. We take n_M=300, broadly matching observations from microscopy (Chustecki et al., 2021, 2025; Logan, 2010). Mixing of mitochondrial matrix content through the cellular population (through fusion) occurs on the timescale of 1-2 h in onion epidermis (Arimura et al., 2004). This corresponds to every mitochondrion having undergone at least one fusion event on average. The average number of fusion events per mitochondrion is 2 κ t /〖 n〗_M, so κ=n_M/2=150 τ^(-1) achieves this mixing rate. Across a collection of functions and families, plant proteins have a mean turnover rate around 0.1 per day, so we use this turnover rate for our model proteins (Nelson et al., 2013). 
DNA damage rate estimates involve a variety of estimates, sometimes propagated through the scientific literature from early sources (see Supplementary Information). At an order-of-magnitude level (Johnston et al., 2014), one single-strand break per mtDNA per day is consistent with several of these arguments. We adopt this damage rate by default for our model mtDNA, and will allow it to vary in different instances of our simulations.
The BER pathway reported in (Ferrando et al., 2019; Jeppesen et al., 2011) involves a collection of five steps catalysed by protein complexes, the least numerous subunits of which are often present at low copy numbers (Supplementary Information) around 4 per mitochondrion in (Fuchs et al., 2020). A key player in TR, MSH1, is present at a copy number around 2 per mitochondrion and likely acts in dimer form (Fuchs et al., 2020). We set the production and import rates λ and ρ in our model to produce these average per-mitochondrion copy numbers.
We report results after 1000 timesteps in our simulations, having found that system behaviour is stabilised at this time (Supp. Fig. 1).

