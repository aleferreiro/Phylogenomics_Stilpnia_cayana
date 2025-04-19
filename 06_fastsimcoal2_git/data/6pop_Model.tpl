//Parameters for the coalescence simulation program : fastsimcoal.exe
6 samples to simulate :
//Population effective sizes (number of genes)
$NGUY$
$NPAR$
$NMJO$
$NCAA$
$NEST$
$NWST$
//Samples sizes and samples age 
10
10
6
10
10
10
//Growth rates : negative growth implies population expansion
0
0
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
6
//Migration matrix 0
0 $GUYPAR$ $GUYMJO$ $GUYCAA$ $GUYEST$ $GUYWST$
$GUYPAR$ 0 $PARMJO$ $PARCAA$ $PAREST$ $PARWST$
$GUYMJO$ $PARMJO$ 0 $MJOCAA$ $MJOEST$ $MJOWST$
$GUYCAA$ $PARCAA$ $MJOCAA$ 0 $CAAEST$ $CAAWST$
$GUYEST$ $PAREST$ $MJOEST$ $CAAEST$ 0 $ESTWST$
$GUYWST$ $PARWST$ $MJOWST$ $CAAWST$ $ESTWST$ 0
//Migration matrix 1 
0 0 $GUYMJO1$ $GUYCAA1$ $GUYEST1$ $GUYWST1$
0 0 0 0 0 0
$GUYMJO1$ 0 0 $MJOCAA1$ $MJOEST1$ $MJOWST1$
$GUYCAA1$ 0 $MJOCAA1$ 0 $CAAEST1$ $CAAWST1$
$GUYEST1$ 0 $MJOEST1$ $CAAEST1$ 0 $ESTWST1$
$GUYWST1$ 0 $MJOWST1$ $CAAWST1$ $ESTWST1$ 0
//Migration matrix 2
0 0 $GUYMJO2$ 0 $GUYEST2$ $GUYWST2$
0 0 0 0 0 0 
$GUYMJO2$ 0 0 0 $MJOEST2$ $MJOWST2$
0 0 0 0 0 0
$GUYEST2$ 0 $MJOEST2$ 0 0 $ESTWST2$
$GUYWST2$ 0 $MJOWST2$ 0 $ESTWST2$ 0
//Migration matrix 3 
0 0 0 0 $GUYEST3$ $GUYWST3$
0 0 0 0 0 0
0 0 0 0 0 0
0 0 0 0 0 0
$GUYEST3$ 0 0 0 0 $ESTWST3$
$GUYWST3$ 0 0 0 $ESTWST3$ 0
//Migration matrix 4 
0 0 0 0 $GUYEST4$ 0
0 0 0 0 0 0
0 0 0 0 0 0 
0 0 0 0 0 0
$GUYEST4$ 0 0 0 0 0
0 0 0 0 0 0
//Migration matrix 5
0 0 0 0 0 0 
0 0 0 0 0 0
0 0 0 0 0 0 
0 0 0 0 0 0
0 0 0 0 0 0 
0 0 0 0 0 0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
10 historical event
$TDIV1$ 1 0 1 1 0 1 // Coalescence GUYPAR
$TDIV1$ 1 1 0 0 0 1 // Kill PAR deme
$TDIV2$ 3 4 1 1 0 2 // Coalescence ESTCAA
$TDIV2$ 3 3 0 0 0 2 // Kill CAA deme
$TDIV3$ 2 0 1 1 0 3 // Coalescence GUYMJO
$TDIV3$ 2 2 0 0 0 3 // Kill MJO deme
$TDIV4$ 5 4 1 1 0 4 // Coalescence ESTWST
$TDIV4$ 5 5 0 0 0 4 // Kill WST deme
$TDIV5$ 0 4 1 1 0 5 // Coalescence GUYEST
$TDIV5$ 0 0 0 0 0 5 // Kill GUY deme
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0 4.6e-9 OUTEXP
