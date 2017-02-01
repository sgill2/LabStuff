import pdbfixer
from pdbfixer import PDBFixer
from simtk import unit
from simtk.openmm.app import PDBFile
output_file = 't_h.pdb'
fixer = PDBFixer(filename='VER_apo.pdb')
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(pH=7.5)
fixer.addSolvent(padding=11*unit.angstrom, ionicStrength=0.050*unit.molar)
#PDBFile.writeHeader(fixer.topology, open(output_file, 'w'))                     
PDBFile.writeFile(fixer.topology, fixer.positions, open(output_file, 'w'))
#PDBFile.writeFooter(fixer.topology, open(output_file, 'a'))      



