//Parameters for the coalescence simulation program : simcoal.exe
3 samples to simulate :
//Population effective sizes (number of genes)
NPOPcit
NPOPlab
NPOPasl
//Samples sizes and samples age 
20
20
20
//Growth rates: negative growth implies population expansion
Rcit
Rlab
Rasl
//Number of migration matrices : 0 implies no migration between demes
1
//Migration matrix 0
0 MIG_gl 0
MIG_gl 0 0
0 0 0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
8 historical event
TADMIX 2 0 MIGRANTS_cit 1 keep 0
TADMIX 2 1 MIGRANTS_lab 1 keep 0
TDIV_col 2 2 0 1 0 0
TDIV_col 2 1 CONT_lab 1 keep 0
TDIV_col 2 0 1 1 keep 0
TDIV_gl 1 1 0 1 0 nomig
TDIV_gl 1 0 1 1 keep nomig
TCHANGE_cit 0 0 0 RESIZE_cit 0 nomig
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0 3.5e-9
