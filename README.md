# CQBP_Python_Package
This documentation provides an overview of the structural modules inclusive in the TSAMS' Classical Quantum Python Suite for Biophysics.


READ ME
Classical Quantum Biophysics (CQBP)
A comprehensive framework for modeling biological systems using principles from classical and quantum physics, with applications in drug discovery and biophysical modeling.
Overview
The Classical Quantum Biophysics (CQBP) package provides a unified framework for studying biological systems across scales, from quantum effects in enzyme catalysis to classical molecular dynamics and pharmacokinetic modeling. CQBP bridges the gap between quantum and classical descriptions of biological processes, enabling researchers to explore phenomena that span these domains.
Key Features
Core Framework 
• Field-based representations of molecules and their interactions 
• Quantum mechanical operators for modeling quantum effects in biological systems 
• Classical molecular systems for traditional molecular modeling 
• Hybrid quantum-classical approaches for multiscale modeling
Quantum Capabilities 
• Quantum dynamics simulations for studying time evolution of quantum systems 
• Path integral methods for incorporating nuclear quantum effects 
• Open quantum system dynamics using Lindblad master equation 
• Feynman path integral calculations for quantum tunneling and coherence
Molecular Modeling 
• Molecular structure representations with atoms, bonds, and properties 
• Protein structure modeling with residues, chains, and secondary structure 
• Membrane modeling for lipid bilayers and membrane proteins 
• Molecular field calculations for electrostatics, sterics, and hydrophobicity
Drug Discovery Applications 
• Drug-target interaction analysis using field-based approaches 
• Virtual screening of compound libraries against protein targets 
• Pharmacophore modeling for identifying key interaction features 
• PBPK modeling for simulating drug distribution in the body
Installation
Clone the repositories (temporarily individualized while under development)
git clone https://github.com/ctibedoJ/cqbp 
Install the package
pip install -e .
Dependencies • NumPy (>=1.19.0) • SciPy (>=1.5.0) • Matplotlib (>=3.3.0) • NetworkX (>=2.5.0) • scikit-learn (>=0.24.0)
Usage Examples
Quantum Dynamics Simulation
import numpy as np 
from cqbp.quantum import QuantumState, Hamiltonian 
from cqbp.simulations import QuantumDynamics
Define Pauli matrices
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex) sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)
Create a Hamiltonian
H_matrix = 1.0 * sigma_z + 0.5 * sigma_x hamiltonian = Hamiltonian(H_matrix)
Create an initial state
psi_0 = np.array([1.0, 0.0], dtype=complex) initial_state = QuantumState(psi_0)
Create a quantum dynamics simulator
simulator = QuantumDynamics(method="unitary", solver="exact") simulator.initialize((hamiltonian, initial_state))
Run the simulation
results = simulator.run( steps=100, dt=0.1, observe_every=1, observers=[lambda sys, t: np.real(simulator.get_expectation_value(sigma_z))] )
Drug-Target Interaction Analysis
from cqbp.molecular import Molecule, Protein 
from cqbp.applications import DrugTargetInteraction
Load molecules
drug = Molecule.from_file("drug.mol") protein = Protein.from_file("protein.pdb")
Create drug-target interaction object
interaction = DrugTargetInteraction(drug, protein)
Define binding site
binding_site = [45, 46, 49, 52, 73, 76] # Residue indices interaction.set_binding_site(binding_site)
Calculate molecular fields
interaction.calculate_fields(["electrostatic", "steric", "hydrophobic"])
Calculate field similarity
similarities = interaction.calculate_field_similarity("tanimoto")
Dock the drug to the target
poses = interaction.dock_drug(method="field_based", n_poses=10)
Calculate binding energy
binding_energy = interaction.calculate_binding_energy()
Analyze interactions
interactions = interaction.analyze_interactions()
Path Integral Molecular Dynamics
from cqbp.molecular import MolecularSystem 
from cqbp.simulations import PathIntegral
Create a molecular system
system = MolecularSystem.from_file("molecule.xyz")
Create a path integral simulator
pi_simulator = PathIntegral(name="PIMD", n_beads=32, temperature=300.0) pi_simulator.initialize(system)
Define observer functions
def position_observer(system, time): return pi_simulator._get_state()["centroid_positions"][0, 0]
Run the simulation
results = pi_simulator.run( steps=1000, dt=0.001, observe_every=10, observers=[position_observer] )
Module Structure • core: Base classes for fields, operators, and systems 
• fields: Molecular field representations and transformations 
• quantum: Quantum systems, density matrices, and Hamiltonians 
• molecular: Molecule and protein representations 
• membrane: Membrane and lipid bilayer models 
• simulations: Simulation methods including MD, QD, and PIMD 
• applications: Drug discovery and PBPK modeling applications
Examples
The examples directory contains Jupyter notebooks demonstrating various features of the CQBP package:
• quantum_dynamics_example.ipynb: Quantum dynamics simulations • drug_discovery_example.ipynb: Drug discovery applications • path_integral_example.ipynb: Path integral simulations for quantum effects
License
MIT License
Citation
If you use CQBP in your research, please cite:
ctibedo@gmail.com {cqbp2025, author = {Charles Tibedo}, title = {CQBP: Classical Quantum Biophysics}, year = {2025}, url = {https://github.com/ctibedoJ/cqbp} }
Contributing
Contributions are welcome! Please feel free to submit a Pull Reque
