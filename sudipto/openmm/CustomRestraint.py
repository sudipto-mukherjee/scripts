from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

pdb = PDBFile('input_exercise4.pdb')
forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
print "Creating System..."
system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff)
print "Creating external potential.."
force = CustomExternalForce('10*max(0, r-1)^2; r=sqrt(x*x+y*y+z*z)')
for i in range(system.getNumParticles()):
    force.addParticle(i, ())
system.addForce(force)
integrator = LangevinIntegrator(1000*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
print "Using Platform:", simulation.context.getPlatform().getName()
simulation.context.setPositions(pdb.positions)
simulation.reporters.append(PDBReporter('output_exercise4.pdb', 100))
print "Adding Reporters to report Potential Energy and Temperature every 100 steps"
simulation.reporters.append(StateDataReporter(stdout, 100, step=True,
potentialEnergy=True, temperature=True))
simulation.step(1000)
print "Finished Simulation."

