# LoopExtrusion

A tool to simulate chromatin loop extrusion with LAMMPS and estimate loop state proportions from simulations (closed vs extruding vs open)

This repository contains three directories:
 - **simulations** used to simulate a polymer with or without loop extrusion.
 - **simulated_data_subsample** contains a few polymer simulations with and without loop extrusion.
 - **Anchor-anchor_analysis** used to analyze simulated loop anchor coordinates to estimate loop state proportions (closed vs extruding vs open) from simulated static imaging.

Licence: MIT

## Simulations

The directory **simulations** contains four directories, covering the following cases:
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



## Simulated data

The directory **simulated_data_subsample** contains the coordinates of 2 anchors (beads #275 and #324) encompassing a loop of 150kb. In **simulated_data_subsample/Free**, the polymer chain was free (no loop), while the polymer was submitted to loop extrusion in **simulated_data_subsample/Loop**. In simulations with extrusion, the polymer chain was simulated such as it went through 3 different states : i) Open state  (absence of loops) ii) Extruding state where the loop size increases with time (anchor-anchor distance decreases) and iii) closed state corresponding to a stable loop with the two anchors in contact.

In each .txt file, the first three columns correspond to the XYZ coordinates of the anchor, the fourth column indicates the state label (0=Open, 1=Extruding, 2=Closed). Each row corresponds to a simulation timepoint (2991 timepoints in polymers submitted to loop extrusion). The simulation ID is indicated at the end of each .txt file. 

Only 19 samples were provided here. The full dataset is available at : 10.5281/zenodo.7521386.

## Anchor-anchor analysis

The directory **Anchor-anchor_analysis** contains a python script to estimate the proportion of Closed, Extruding and Open states from simulated static microscopy imaging of loop anchors.

From the simulated anchor coordinates given in **simulated_data_subsample** directory, the script first computes the difference distribution of the two anchors. Assuming Normal distributions for Closed and Open states, the script then estimates the standard deviations for these 2 states ( $\sigma_{closed}$ and $\sigma_{free}$ ). No parameter is estimated for the Extruding state since it can be modeled as an integral of Normal distribution with standard deviation varying from $\sigma_{closed}$ to $\sigma_{free}$. We end up with a 3-state model.

Finally, after simulating anchor-anchor difference distribution with different proportions of Open, Extruding and Closed states (with additive Gaussian noise to simulate localizations errors), the script estimates these proportions by fitting to the distribution our 3-states model.

### Requirements
The analysis script has been run using python 3.9 with the following libraries: matplotlib, scipy, numpy, pandas.
