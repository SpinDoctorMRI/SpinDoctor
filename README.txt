# SpinDoctor NeuronModule

SpinDoctor is a software package that performs numerical simulations of diffusion magnetic resonance imaging (dMRI) for prototyping purposes.  
The user is advised to read the latest version from \url{https://github.com/jingrebeccali/SpinDoctor} 

The NeuronModule solves the Bloch-Torrey PDE on neurons. The PDEs are solved by P1 finite elements combined with built-in MATLAB routines for solving ordinary differential equations.  The finite element mesh generation is performed using an external package called TetGen that is included in the toolbox.

Paper about SpinDoctor can be found at https://arxiv.org/abs/1902.01025

Paper about SpinDoctor NeuronModule can be found at http://arxiv.org/abs/1910.07916

Website of SpinDoctor can be found at http://www.cmap.polytechnique.fr/~jingrebeccali/software.html

Software requirements: SpinDoctor NeuronModule is compatible with MATLAB (version R2017b or later).

Getting started
1) The DISTRIBUTE folder contains a commented general purpose driver called driver_spindoctor_neuronmodule_commented.m. It is highly recommended to read this driver to understand the workflow of SpinDoctor. 
2) driver_spindoctor_neuronmodule_commented.m does not use saved simulation data, all simulations are run from scratch.
3) Thirteen driver examples are given in the DISTRIBUTE folder, the simulation details are listed in the table below. These drivers save the simulated data or use previously saved simulation results if they are available.
4) The functions that the user is likely to call directly from the driver are located at the top level of DISTRIBUTE/SRC.
5) Other functions are stored in subfolders of DISTRIBUTE/SRC.
6) Some neuron finite element meshes are stored in DISTRIBUTE/msh_files.
7) The input files of the drivers are located at DISTRIBUTE/params_files.
8) The saved simulation data are saved in DISTRIBUTE/saved_simul.
9) Documentation of the top level functions can be found in DOC/.

| Driver name                         	| Simulation details                                                                       	|
|-------------------------------------	|------------------------------------------------------------------------------------------	|
| driver_btpde_multidir_soma_OGSE.m   	| BTPDE simulation of 03b_spindle4aACC_soma in 30 gradient directions with OGSE sequences. 	|
| driver_btpde_multidir_dendrites.m   	| BTPDE simulation of 03b_spindle7aACC_dendrites_1 in 30 gradient directions.              	|
| driver_btpde_multidir_soma.m        	| BTPDE simulation of 03b_spindle4aACC_soma in 60 gradient directions.                     	|
| driver_btpde_multidir_wholeneuron.m 	| BTPDE simulation of 25o_pyramidal18aFI in 80 gradient directions.                        	|
| driver_btpde_onedir_dendrites.m     	| BTPDE simulation of 03b_spindle7aACC_dendrites_1 in one gradient direction.              	|
| driver_btpde_onedir_soma.m          	| BTPDE simulation of 03b_spindle7aACC_soma in one gradient direction.                     	|
| driver_btpde_onedir_wholeneuron.m   	| BTPDE simulation of 25o_pyramidal18aFI in one gradient direction.                        	|
| driver_hadc_multidir_dendrites.m    	| HADC simulation of 03b_spindle7aACC_dendrites_1 in 60 gradient directions.               	|
| driver_hadc_multidir_soma.m         	| HADC simulation of 03b_spindle4aACC_soma in 30 gradient directions.                      	|
| driver_hadc_multidir_wholeneuron.m  	| HADC simulation of 25o_pyramidal18aFI in 30 gradient directions.                         	|
| driver_hadc_onedir_dendrites.m      	| HADC simulation of 03b_spindle7aACC_dendrites_1 in one gradient direction.               	|
| driver_hadc_onedir_soma.m           	| HADC simulation of 03b_spindle4aACC_soma in one gradient direction.                      	|
| driver_hadc_onedir_wholeneuron.m    	| HADC simulation of 25o_pyramidal18aFI in one gradient direction.                         	|

Authors: Jing-Rebecca Li, Van-Dang Nguyen, Chengran Fang, Try Nguyen Tran.

Copyright (C) 2019

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
