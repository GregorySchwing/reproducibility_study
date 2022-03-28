from numpy import zeros
import numpy as np
import pandas as pd
from parmed.charmm import CharmmPsfFile
from parmed.formats.registry import load_file

#sel = Select()
#ResiduesToStrip = sel.getBoolArray(pos2comb, "same residue as (exwithin 2.4 of protein)")

solvbox2 = CharmmPsfFile("mosdef_box_0.psf")
solvbox2_pos = load_file("mosdef_box_0.pdb")
#, positions="mosdef_box_0.psf
pro2_parm = CharmmPsfFile("6g6k.psf")

pro2_parm_pos = load_file("6g6k_aligned.pdb")
parmPosComb = pro2_parm_pos + solvbox2_pos
# Doesnt work because prody isnt smart about cmap
#systemPSF = CharmmPsfFile.from_structure("intermediate.psf")
systemPSF = CharmmPsfFile.from_structure(pro2_parm + solvbox2)
systemPSF.write_psf("parmedComb.psf")
parmPosComb.write_pdb("parmedComb.pdb", renumber=False)
#systemPSF.write_pdb("parmed.pdb")
#strippedPSF=systemPSF.strip(ResiduesToStrip)
systemPSFComb = CharmmPsfFile("parmedComb.psf")
systemPSFComb_pos = load_file("parmedComb.pdb")

df = systemPSFComb_pos.to_dataframe()
df2 = systemPSFComb.to_dataframe()

print(df2.head())

endOfProtein = df2.segid.eq('SYS').idxmax()
endOfSolvent = len(df2.index)
startOfSolvent = endOfProtein+1
prot = "@0-{end}".format(end=endOfProtein)
solv = "@{start}-{end}".format(start=startOfSolvent, end=endOfSolvent)
print(prot)
print(solv)


amberMaskRes = "{protS}<:2.4&{solvS}".format(protS=prot, solvS=solv)
amberMaskRes = ":1-10"
badWatersRes = systemPSFComb_pos[amberMaskRes]


bwDF = badWatersRes.to_dataframe()

JustResFull = df2[["name","resname","resnum"]]
JustResPart = bwDF[["name","resname","resnum"]]

print(JustResFull.head())
print(len(JustResFull.index))
print(JustResPart.head())
print(len(JustResPart.index))
df = pd.merge(JustResFull, JustResPart, on=["name","resname","resnum"], how='left', indicator='Exist')
df['Exist'] = np.where(df.Exist == 'both', True, False)
stripInput = df['Exist'].to_numpy()
print(len(stripInput))
#maybe = stripInput.to_list()

splitRes = systemPSFComb.strip(stripInput)
print(splitRes)
systemPSFComb.write_psf("parmedComb2.psf")

splitRes = systemPSFComb_pos.strip(stripInput)
print(splitRes)
systemPSFComb_pos.write_pdb("parmedComb2.pdb")
# Select all residues which are less than 2.4 from a protein atom
# Also, solvent residue atoms must be contained in solvent atoms.
# For a residue to be included in this mask, only one atom needs to be below the cutoff!
# This is the desired behavior.

# To see this uncomment these lines and compare
#amberMask = "{protS}<@2.4&{solvS}".format(protS=prot, solvS=solv)
#badWaters = systemPSFComb_pos[amberMask]
#print(badWatersRes.to_dataframe())
#print(badWaters.to_dataframe())




