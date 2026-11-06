# Master-s-Project
Modelling Cardiac Electrophysiology for 4th Year Integrated Master's Project at Durham University

To initialse environment run:
conda activate ppp
in terminal to use correct Python environment for this project

File Types:
1. Function Files - these contain all functions needed to run a specific model and will be named after the model they are taken from 
2. Trial Files - testing the function files to produce graphs and check that all functions run as expected, any other forms of testing will be named as Trial Files 
3. Trial Files - produce the results expected from models and generate graphs that can be used in report.
4. Cellml Files - downloaded from online sources, translations stored in other files to allowing running in Python

Trovato Purkinje Model:
- Somehow there is inconsistencies between the supplementary material equations and the cellml code repository for this model
- I have formulated functions for both the cellml version and the supplementary material to perform a comparison and determine which version shoul be kept 
- Might be best to adapt this model for a type of animal cell so there is more accesible experimental data
- Supplementary equations don't seem to be correct upon further investigation so will continue with cellml file for now 


TWorld Model:
- taken cellml code from the open source online github repository
- translation from cellml performed by generative AI
- comparisons made using the simulated conditions in peer-reviewed journal for TWorld model to check implementation in Python 

Bidomain Model:
- constructing for the Purkinje fibre in 1D currently 
- used matrix formulation and built in solvers for matrix equations to evaluate variables from derivative approximations
- should be possible to use with both Purkinje and ventricular models
- bidomain function returns the new membrane voltage for all cellular positions when it is advanced one step in time


Purkinje_functions_tissues:
- currently works for 1D fibres 
- in the 'Trovato 1D Fibre.py' file I have defined a function that runs and plots the Purkinje fibre results of the simulation
- will need to adapt the code from the 1D fibre for plots and functions is other dimensions to account for vector propagation in 2D
