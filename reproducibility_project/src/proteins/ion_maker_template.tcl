package require psfgen

topology PATH_2_IONS_TOPOLOGY

segment CAT {
    pdb SINGLE_CATION_PDB.pdb
    first none
    last none
}

coordpdb CATION_NAME.pdb CAT

writepsf CATION_NAME.psf

resetpsf

segment ANI {
    pdb SINGLE_ANION_PDB.pdb
    first none
    last none
}

coordpdb ANION_NAME.pdb ANI

writepsf ANION_NAME.psf

