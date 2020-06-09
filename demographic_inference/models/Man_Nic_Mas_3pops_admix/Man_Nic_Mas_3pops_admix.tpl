//Parameters for the coalescence simulation program : simcoal.exe
3 samples to simulate :
//Population effective sizes (number of genes)
NPOPman
NPOPnic
NPOPmas
//Samples sizes and samples age 
20
20
20
//Growth rates: negative growth implies population expansion
Rman
Rnic
Rmas
//Number of migration matrices : 0 implies no migration between demes
1
//Migration matrix 0
0 MIG_gl 0
MIG_gl 0 0
0 0 0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
8 historical event
TADMIX 2 1 MIGRANTS_nic 1 keep 0
TADMIX 2 0 MIGRANTS_man 1 keep 0
TDIV_col 2 2 0 1 0 0
TDIV_col 2 0 CONT_man 1 keep 0
TDIV_col 2 1 1 1 keep 0
TDIV_gl 0 0 0 1 0 nomig
TDIV_gl 0 1 1 1 keep nomig
TCHANGE_nic 1 1 0 RESIZE_nic 0 nomig
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0 3.5e-9
