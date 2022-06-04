package require psfgen

topology PATH_2_IONS_TOPOLOGY

segment CAT {
    pdb CATION_NAME.pdb
    first none
    last none
}

writepsf ./ANION_NAME.psf

segment ANI {
    pdb ANION_NAME.pdb
    first none
    last none
}

writepsf ./ANION_NAME.psf

