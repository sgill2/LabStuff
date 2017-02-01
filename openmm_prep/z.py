from simtk import openmm, unit
from simtk.openmm import app
import mdtraj as md
from openmmtools import testsystems
import math
import numpy as np
from mdtraj.reporters import HDF5Reporter
from simtk import openmm, unit

#parameters
temperature = 300.0 * unit.kelvin
stepsNPT = 100000
stepsNVT = 250000

coord_file = 'eqToluene.inpcrd'
top_file =   'eqToluene.prmtop'
if 1:
#    pdb_file = 'Hsp_h.pdb'
#    pdb_file = 'pdbfix_water.pdb'
    pdb_file = 'combined2.pdb'


    pdb = app.PDBFile(pdb_file)
    forcefield = app.ForceField("amber99sbildn.xml", "lig2.xml", "tip3p.xml")
#    forcefield = app.ForceField("amber99sbildn.xml", "tip3p.xml")

    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, nonbondedCutoff=1.0 * unit.nanometers, constraints=app.HBonds)
#if getting strucutre from pdb


md_integrator = openmm.openmm.LangevinIntegrator(temperature, 1/unit.picosecond, 0.002*unit.picoseconds)
print('minizmizing')
#select if min is not presnet, otherwise 0
if 1:
    md_simulation = openmm.app.simulation.Simulation(topology=pdb.topology, system=system, integrator=md_integrator)
    md_simulation.context.setPositions(pdb.positions)
    openmm.LocalEnergyMinimizer.minimize(md_simulation.context)
    md_simulation.saveState('min.xml')
else:
    md_simulation = openmm.app.simulation.Simulation(topology=pdb.topology, system=system, integrator=md_integrator, state='min.xml')

md_simulation.context.setVelocitiesToTemperature(temperature)
md_simulation.reporters.append(openmm.app.dcdreporter.DCDReporter('npt.dcd', 100000))
md_simulation.reporters.append(HDF5Reporter('npt.h5', 100000))
md_simulation.reporters.append(openmm.app.statedatareporter.StateDataReporter('nptinfo.csv', 100000, step=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True))

if 1:
    print('Equilibrating Volume')
    systemFromContext = md_simulation.context.getSystem()
    print systemFromContext.getDefaultPeriodicBoxVectors()

    md_simulation.step(stepsNPT)
    temp_xml = 'temp.xml'
    barostat_index = systemFromContext.getNumForces() -1
    md_simulation.system.removeForce(barostat_index)

#    md_simulation.saveState('temp.xml')
    stateinfo = md_simulation.context.getState(True, True, True, True, True, True)
    statepos = stateinfo.getPositions(asNumpy=True)[:]
    statevel = stateinfo.getVelocities(asNumpy=True)[:]
    boxVectors = stateinfo.getPeriodicBoxVectors(asNumpy=True)[:]
    md_simulation.context.reinitialize()
    md_simulation.context.setPositions(statepos)
    md_simulation.context.setVelocities(statevel)
    md_simulation.context.setPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2])
    md_simulation.saveState(temp_xml)

if 1:
    print('Equilbrating NVT')
    print('this selected')
#    stateinfo = md_simulation.context.getState(True, True, True, True, True)
#    statepos = stateinfo.getPositions(asNumpy=True)[:]
#    boxVectors = stateinfo.getPeriodicBoxVectors()[:]
    nvt_integrator = openmm.openmm.LangevinIntegrator(temperature, 1/unit.picosecond, 0.002*unit.picoseconds)

    nvt_simulation = openmm.app.simulation.Simulation(topology=pdb.topology, system=system, integrator=nvt_integrator, state=temp_xml)
    nvt_simulation.reporters.append(openmm.app.dcdreporter.DCDReporter('nvt.dcd', 100000))
    nvt_simulation.reporters.append(HDF5Reporter('nvt.h5', 100000))
    nvt_simulation.reporters.append(openmm.app.statedatareporter.StateDataReporter('nvtinfo.csv', 50000, step=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True))


    nvt_simulation.step(stepsNVT)
