RDMSV
=====

Drug Markov Singular Values of Transition Probabilities in R

Cristian R Munteanu, BiGCaT, Maastricht University (UM), Netherlands

Haralambos Sarimveis, National Technical University of Athens (NTUA), Greece 

Contact: muntisa [at] gmail [dot] com

This tool is using the same flow as the RDMMP or MInD-Prot but it introduces a new class of descriptors such as the Markov Singular Values for the transition probabilities of the molecular graph (the minimum and the maximum values).

The inputs are SMILES formulas and and 4 atom physical - chemical properties such as Mulliken electronegativity (EM), Kang-Jhon polarizability (PKJ), van der Waals area (vdWA) and atom contribution to partition coefficient (AC2P). There are 6 types of atom types: All (all atoms), Csat (saturated C), Cinst (insaturated C), Hal (halogen), Het (heteroatoms) and HetNoX (heteroatoms but not halogens. The list of atom properties could be extended.

RDMSV can calculate:
- 48 descriptors - if the flag fAllKs = 0 (default), that includes only the averaged descriptors over all ks for each atomic property and atomic type 
- 336 descriptors - if the flag fAllKs = 1, that includes all the descriptors for all k power values, atomic properties and atom type + the averaged values over all ks for each atomic property and atomic type.


Script steps:
- Read the inputs: SMILES formulas and atom properties
- Get connectivity matrix (CM), nodes = atoms, edges = chemical bonds
- Get weights vector (w) for each atom property
- Calculate weighted matrix (W) using CM and w
- Calculate transition probability (P) based on W
- Calculate k powers of P -> Pk matrices
- Calculate Markov Singurla Values of Transition Probability Matrix for all power, each type of atom property and atom type (288 descriptors)
- Calculate the average values over all K values (48 descriptors)
- The total number of descriptor is 336
 -If the flag fAllKs is 1 -> 336 descriptors; if 0 -> only 48 (averaged values)
- Descriptors for all k values -> 288 = (4 properties * 6 atom types * 6 powers) * 2 for Min and Max
- Averages descriptors         -> 48  = (4 properties * 6 atom types) * 2 for Min and Max
