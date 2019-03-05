This repository includes some C++ codes for variety of different purposes in the scope of data generation to derive Machine Learning Potentials. The tasks are briefly as follow:

lj_fit.cpp --> This piece of code takes lattice constant and atomic mass as input, generate POSCAR files of perturbed atomic structures for DFT calculations, and then fits Lennard-Jones 6-12 parameters using Nelder-Mead minimization over the data provided by OUTCAR files 

ae_atomic_forces.cpp --> Reads a .cfg data file and write the atomic forces in an output file for visualization

cfg_to_poscar.cpp --> Reads configurations from a .cfg file and write POSCAR files for each atomic environment configurations

data_read_try.cpp --> Find the max energy atom for a given LAMMPS dump output and returns its atom ID

mlip_energy_force_correlation.cpp --> Reads the given .cfg file and writes forces and energies for correlation plots

mlip_parameters_errors.cpp --> Reads force, energy and stress errors from an error file and generate data for heat-map visualization

select_ae.cpp --> Select the atoms from different regions of a given configuration(dump) file, and write the corresponding atomic environment configurations in a POSCAR file

select_generate_poscars.cpp --> This code takes lattice constant, atomic mass and number of desired simulations, generates rattled atomic structures and write relevant POSCARs to generate DFT data

write_outcar_mlip.cpp --> Writes .cfg database file based on given set of OUTCAR files


