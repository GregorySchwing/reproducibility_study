package require psfgen

topology PATH_2_IONS_TOPOLOGY

segment CATION_NAME {
    pdb CATION_NAME.pdb
    first none
    last none
}

writepsf ./ANION_NAME.psf

segment ANION_NAME {
    pdb ANION_NAME.pdb
    first none
    last none
}

writepsf ./ANION_NAME.psf

