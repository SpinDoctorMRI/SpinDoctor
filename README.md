# SpinDoctor Neuron Module

==================================

This branch contains code and documentation to run some of the examples from the paper 

C. Fang, V.-D. Nguyen, D. Wassermann, J.-R. Li  
[Diffusion MRI simulation of realistic neurons with SpinDoctor and the Neuron Module](https://doi.org/10.1016/j.neuroimage.2020.117198)    
Neuroimage. 2020.  

==================================

SpinDoctor is a software package that performs numerical simulations of diffusion magnetic resonance imaging (dMRI) for prototyping purposes.  

The Neuron Module solves the Bloch-Torrey PDE on neurons. 

Software requirements

	SpinDoctor NeuronModule is compatible with MATLAB (version R2017b or later).

Getting started
1. The DISTRIBUTE folder contains a commented general purpose driver called driver_spindoctor_neuronmodule_commented.m. It is highly recommended to read this driver to understand the workflow of SpinDoctor. 
2. driver_spindoctor_neuronmodule_commented.m does not use saved simulation data, all simulations are run from scratch.
3. Thirteen driver examples are given in the DISTRIBUTE folder, the simulation details are listed in the table below. These drivers save the simulated data or use previously saved simulation results if they are available.
4. Six drivers that are used to plot the figures in the paper about SpinDoctor NeuronModule by utilizing the saved data are also given in the repository. The details can be found in the table below.
5. The functions that the user is likely to call directly from the driver are located at the top level of DISTRIBUTE/SRC.
6. Other functions are stored in subfolders of DISTRIBUTE/SRC.
7. Some neuron finite element meshes are stored in DISTRIBUTE/msh_files.
8. The input files of the drivers are located at DISTRIBUTE/params_files.
9. The saved simulation data are saved in DISTRIBUTE/saved_simul.
10. Documentation of the top level functions can be found in DOC/.
11. User Guide can be found [here](https://github.com/jingrebeccali/SpinDoctor/blob/NeuronModule/User%20Guide.pdf.

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
| driver_figure6.m                    	| Plot the Figure 6 in the paper about SpinDoctor NeuronModule.                            	|
| driver_figure7_8_9.m                	| Plot the Figure 7, 8 and Figure 9 in the paper about SpinDoctor NeuronModule.            	|
| driver_figure10.m                   	| Plot the Figure 10 in the paper about SpinDoctor NeuronModule.                           	|
| driver_figure11.m                   	| Plot the Figure 11 in the paper about SpinDoctor NeuronModule.                           	|
| driver_figureB14_15_16.m            	| Plot the Figure B14, B15 and Figure B16 in the paper about SpinDoctor NeuronModule.      	|
| driver_figureB17.m                  	| Plot the Figure B17 in the paper about SpinDoctor NeuronModule.                          	|

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
