package require psfgen

topology ./ions.str

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

