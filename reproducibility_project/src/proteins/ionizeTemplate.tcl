package require cionize
namespace import ::cionize::* 

source HELPER_METHODS

set numCatAn [getNumIons PATH_2_SOLVATED_PSF PATH_2_SOLVATED_PDB CATION_NAME CATION_VAL ANION_NAME ANION_VAL SALC_CONC_VALUE]

puts "$numCatAn"
return "$numCatAn"

# load solute molecule
set solute [mol new $solvatedProtein.psf waitfor all]
mol addfile $solvatedProtein.pdb mol $solute waitfor all

cionize -mol 0 -np 4 -mg -prefix default -ions "{SOD 85 1} {CLA 104 -1}" 
