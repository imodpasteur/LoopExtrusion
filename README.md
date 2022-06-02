# LoopExtrusion
A tool to simulate loop extrusion in chromatin

This repository contains four directories, covering the following cases:
 - **No_extrusion** used to model the dynamics of a polymer without loop extrusion.
 - **Loop300kb** used to model a loop of 300 kb with random loading of cohesin.
 - **Loop150kb** used to model a loop of 150 kb with random loading of cohesin.
 - **Loop150kb_Unloading** used to model a loop of 300 kb with random loading of cohesin and unloading of cohesin after a time spent in the closed state drawn from an exponential distribution.

Each directory contains the following:
 - **initial.py**: Generates a file "initial_conformation.txt" containing the atoms, bonds and angles used in the simulations as well as the initial coordinates of the beads and the size of the boundary.
 - **input_data_generate_f1.py**: Defines the length L of the polymer, the index of the bead "start" on which the extrusion complex will land (randomly chosen between the anchors), the anchors "anchor_left" and "anchor_right" bead indexes, the random seeds "seed_velocity" and "seed_langevin" and the radius "size_nuc" of the confinement sphere. Variables i and j are updated during the simulation (**lammps_input_file.in**) and represent the beads that get in close contact during extrusion.
 - **interactions**: Defines the potentials between the beads. Rigidity is modeled by the harmonic angle_style.
 - **lammps_input_file.in**: Defines the time course and events of the molecular dynamics simulation steps.
 - **run.sh**: Set the number N of independent simulations to run. Create N different folders and copy each file in each folder.

To model different sizes of loops, the "anchor_left" and "anchor-right" variables in the **input_data_generate_f1.py** script were modified.
To model unloading of cohesin, a variable "death" in the **input_data_generate_f1.py** script was drawn from an exponential distribution, which defined the length (in simulation steps) of the closed state (anchors kept close to each other).

All files were put in the same folder and the simulations were started in parallel using ./run.sh command on the Institut Pasteur cluster. The following packages were used to run the simulations: lammps/2016.11.17 ; gcc/9.2.0 ; openmpi/4.0.5 ; Python/3.8.3.

The resulting .dcd files containing beads coordinates were then imported in Python using the MDAnalysis package in python.
