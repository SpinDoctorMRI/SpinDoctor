# SpinDoctor Matrix Formalism Module

==================================

This branch contains code and documentation dating from 2019, which can be used to run some of the examples from the paper 

J.-R. Li, T. N. Tran, V.-D. Nguyen
[Practical computation of the diffusion MRI signal of realistic neurons based on Laplace eigenfunctions](https://doi.org/10.1002/nbm.4353)
NMR in Biomedicine. 2020


==================================

SpinDoctor is a software package that performs numerical simulations of diffusion magnetic resonance imaging (dMRI) for prototyping purposes.  

The Matrix Formalism Module computes a closed form representation of the diffusion MRI signal called Matrix Formalism, that is based on the eigendecomposition of the Laplace operator, defined in the relevant geometry.  

Currently, the Matrix Formalism Module allows the computation of the Matrix Formalism signal and the 
Matrix Formalism Gaussian Approximation signal for realistic neuron geometries (impermeable membranes) with the PGSE sequence.

Matrix Formalism for permeable membranes and for general diffusion-encoding sequences are under development 
and will be released in the future.  


Software requirements
 
	The SpinDoctor Toolbox and the Neuron Module have been developed in the MATLAB R2017b 
	and require no additional MATLAB Toolboxes.  
	
	*****However, the current version of the Matrix Formalism Module requires the MATLAB PDE Toolbox (2017 or later) 
	due to certain difficulties of implementing the matrix eigenvalue solution on a restricted eigenvalue interval.  
	This technical issue will be addressed in a future release.*****      


Getting started
		
1. The DISTRIBUTE folder contains drivers for typical simulations.
2. Some neuron finite element meshes are stored in DISTRIBUTE/msh_files.
3. The input files of the drivers are located at DISTRIBUTE/params_files.
4. The saved simulation data are saved in DISTRIBUTE/saved_simul.
5. The functions that the user is likely to call directly from the driver are located at the top level of DISTRIBUTE/SRC.
6. Other functions are stored in subfolders of DISTRIBUTE/SRC.
7. Documentation of the top level functions can be found in DOC/
8. User Guide is found [here](https://github.com/jingrebeccali/SpinDoctor/blob/MatrixFormalismModule/User%20Guide.pdf)


There are the following provided drivers

 Driver name                         	| Simulation details                                                                       	|
|-------------------------------------	|------------------------------------------------------------------------------------------	|
| driver_spindoctor_MF_1direction_commented.m   |   1 gradient direction, does not use saved simulation data, all simulations are run from scratch. |
| driver_spindoctor_MF_1direction_UseSavedData.m |  1 gradient direction, save the simulated data or use previously saved simulation results if they are available. The minimum length scale of the eigenvalue interval is made slightly larger so the eigendecomposition data take less than 1MB of space. |
| driver_spindoctor_MF_hardi_commented.m   |	HARDI, does not use saved simulation data, all simulations are run from scratch.|
| driver_spindoctor_MF_hardi_UseSavedData.m | HARDI, save the simulated data or use previously saved simulation results if they are available. |

Authors: Jing-Rebecca Li, Try Nguyen Tran, Van-Dang Nguyen. 

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
