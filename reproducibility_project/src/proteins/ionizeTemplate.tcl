package require cionize
namespace import ::cionize::* 

source HELPER_METHODS

set numCatAn [getNumIons PATH_2_SOLVATED_PSF PATH_2_SOLVATED_PDB CATION_NAME CATION_VAL ANION_NAME ANION_VAL SALC_CONC_VALUE]

puts "$numCatAn"

# load solute molecule
set solute [mol new PATH_2_PROT_PSF waitfor all]
mol addfile PATH_2_PROT_PDB mol $solute waitfor all

cionize -mol 0 -np 4 -mg -prefix default -ions "{CATION_NAME [lindex $numCatAn 0] CATION_VAL} {ANION_NAME [lindex $numCatAn 1] ANION_VAL}" 
