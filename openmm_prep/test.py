import simtk.unit as u
from simtk.openmm import app
import simtk.openmm as mm
import mdtraj as md
import numpy as np

ligand_name = "sustiva"
pdb_filename = "t_h.pdb"
#mol2_filename = "./chemicals/%s/%s.mol2" % (ligand_name, ligand_name)
mol2_filename = "lig2.mol2"
xml_filename = "lig2.xml"
output_name = 'combined.pdb'

temperature = 300 * u.kelvin
friction = 0.3 / u.picosecond
timestep = 2.0 * u.femtosecond

traj = md.load(pdb_filename)
pbc_lengths = traj.unitcell_lengths
pbc_vectors = traj.unitcell_vectors
pbc_angles = traj.unitcell_angles
protein_traj = md.load(pdb_filename)
#protein_traj.center_coordinates()

protein_top = protein_traj.top.to_openmm()
protein_xyz = protein_traj.openmm_positions(0)

ligand_traj = md.load(mol2_filename)
#ligand_traj.center_coordinates()

#Move the pre-centered ligand sufficiently far away from the protein to avoid a clash.  
#min_atom_pair_distance = ((ligand_traj.xyz[0] ** 2.).sum(1) ** 0.5).max() + ((protein_traj.xyz[0] ** 2.).sum(1) ** 0.5).max() + 0.3
#ligand_traj.xyz += np.array([1.0, 0.0, 0.0]) * min_atom_pair_distance
forcefield = app.ForceField("amber99sbildn.xml", xml_filename, "tip3p.xml")
pdb = app.PDBFile(pdb_filename)
system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0 * u.nanometers, constraints=app.HAngles)
ligand_xyz = ligand_traj.openmm_positions(0)
ligand_top = ligand_traj.top.to_openmm()

model = app.modeller.Modeller(protein_top, protein_xyz)
model.add(ligand_top, ligand_xyz)
app.PDBFile.writeFile(model.topology, model.positions, open(output_name, 'w'))
a = md.load(output_name)
a = a.center_coordinates()
a.unitcell_lengths = pbc_lengths
a.unitcell_vectors = pbc_vectors
a.unitcell_angles = pbc_angles
edit_length = a.unitcell_lengths *(2./4.)
#edit_length = a.unitcell_lengths

atom_coords = a.xyz
print(atom_coords)
added_length = np.tile( edit_length, (a.n_atoms, 1))
atom_coords = atom_coords + added_length
print(atom_coords)
a.xyz = atom_coords
a.save(output_name)


#model.addSolvent(forcefield, padding=0.4 * u.nanometer)

#system = forcefield.createSystem(model.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0 * u.nanometers, constraints=app.HAngles)
if 0:
    integrator = mm.LangevinIntegrator(temperature, friction, timestep)

    simulation = app.Simulation(model.topology, system, integrator)
    simulation.context.setPositions(model.positions)

    #simulation.minimizeEnergy(tolerance=2.0)

    simulation.reporters.append(app.PDBReporter("simulation.pdb", 250))
    print("running")
    simulation.step(2500)

