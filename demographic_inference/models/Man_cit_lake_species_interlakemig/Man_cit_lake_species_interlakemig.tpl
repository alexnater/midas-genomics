//Parameters for the coalescence simulation program : simcoal.exe
4 samples to simulate :
//Population effective sizes (number of genes)
N_Man_cit
N_Man_lab
N_Nic_cit
N_Nic_lab
//Samples sizes and samples age 
20
20
20
20
//Growth rates: negative growth implies population expansion
R_Man_cit
R_Man_lab
R_Nic_cit
R_Nic_lab
//Number of migration matrices : 0 implies no migration between demes
4
//Migration matrix 0
0.00000 MIG_man MIG_cit 0.00000
MIG_man 0.00000 0.00000 MIG_lab
MIG_cit 0.00000 0.00000 MIG_nic
0.00000 MIG_lab MIG_nic 0.00000
//Migration matrix 1
0.00000 0.00000 MIG_cit 0.00000
0.00000 0.00000 0.00000 0.00000
MIG_cit 0.00000 0.00000 MIG_nic
0.00000 0.00000 MIG_nic 0.00000
//Migration matrix 2
0.00000 MIG_man MIG_cit 0.00000
MIG_man 0.00000 0.00000 0.00000
MIG_cit 0.00000 0.00000 0.00000
0.00000 0.00000 0.00000 0.00000
//Migration matrix 3
0.00000 0.00000 MIG_cit 0.00000
0.00000 0.00000 0.00000 0.00000
MIG_cit 0.00000 0.00000 0.00000
0.00000 0.00000 0.00000 0.00000
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
8 historical event
TDIV_man 1 1 0 1 0 1
TDIV_man 1 0 1 1 keep 1
TDIV_nic 3 3 0 1 0 2
TDIV_nic 3 2 1 1 keep 2
TMIGSTOP 0 0 0 1 keep 3
TDIV_cit 2 2 0 1 0 nomig
TDIV_cit 2 0 1 1 keep nomig
TCHANGE_cit 0 0 0 RESIZE_cit 0 nomig
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0 3.5e-9
