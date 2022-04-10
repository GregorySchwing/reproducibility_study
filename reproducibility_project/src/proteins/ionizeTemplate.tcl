package require cionize
namespace import ::cionize::* 

source HELPER_METHODS

set numCatAn [getNumIons PATH_2_SOLVATED_PSF PATH_2_SOLVATED_PDB CATION_NAME CATION_VAL ANION_NAME ANION_VAL SALC_CONC_VALUE]

puts "$numCatAn"

# load solute molecule
set solute [mol new PATH_2_PROT_PSF waitfor all]
mol addfile PATH_2_PROT_PDB mol $solute waitfor all

set all [atomselect $solute all] 
set minmax [measure minmax $all]
set vals "[lindex [lindex $minmax 0] 0] [lindex [lindex $minmax 0] 1] [lindex [lindex $minmax 0] 2] [lindex [lindex $minmax 1] 0] [lindex [lindex $minmax 1] 1] [lindex [lindex $minmax 1] 2]"
puts "$vals"
cionize -mol 0 -np 4 -mg -prefix default -border 0 -ga test -ions "{CATION_NAME [lindex $numCatAn 0] CATION_VAL} {ANION_NAME [lindex $numCatAn 1] ANION_VAL}" 
