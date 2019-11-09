# SpinDoctor Toolbox

SpinDoctor is a software package that performs numerical simulations of diffusion magnetic resonance imaging (dMRI) for prototyping purposes.

SpinDoctor can be used 

1) to solve the Bloch-Torrey PDE to obtain the dMRI signal 
(the toolbox provides a way of robustly fitting the dMRI signal to obtain the fitted Apparent Diffusion Coefficient (ADC)); 
2) to solve the diffusion equation of the H-ADC model to obtain the ADC;
3) a short-time approximation formula for the ADC is also included in the toolbox for comparison with the simulated ADC.

The PDEs are solved by P1 finite elements combined with built-in MATLAB routines for solving ordinary differential equations.
The finite element mesh generation is performed using an external package called TetGen that is included in the toolbox.

SpinDoctor provides built-in options of including 
1) spherical cells with a nucleus; 
2) cylindrical cells with a myelin layer; 
3) an extra-cellular space (ECS) enclosed either 
a) in a box or b) in a tight wrapping around the cells; 
4) deformation of canonical cells by bending and twisting.  
5) permeable membranes for the BT-PDE (the H-ADC assumes negligible permeabilty).
Built-in diffusion-encoding pulse sequences include the Pulsed Gradient Spin Echo and the Ocsilating Gradient Spin Echo.

# New Nov 2019: SpinDoctor has two modules called [Neuron Module](https://github.com/jingrebeccali/SpinDoctor/tree/NeuronModule) and [Matrix Formalism Module](https://github.com/jingrebeccali/SpinDoctor/tree/MatrixFormalismModule). Details about these modules can be found in the branches with those names.

Paper about SpinDoctor can be found at https://arxiv.org/abs/1902.01025

Software requirements: 

	The SpinDoctor Toolbox and the Neuron Module have been developed in the MATLAB R2017b 
	and require no additional MATLAB Toolboxes.  

	*****However, the current version of the Matrix Formalism Module requires the MATLAB PDE Toolbox (2017 or later) 
	due to certain difficulties of implementing the matrix eigenvalue solution on a restricted eigenvalue interval.  
	This technical issue will be addressed in a future release.*****   

Getting started
1) The DISTRIBUTE folder contains drivers for examples in the SpinDoctor paper as well as a general purpose driver called driver_spindoctor.m
2) The functions that the user is likely to call directly from the driver are located at the top level of DISTRIBUTE/SRC.
3) Other functions are stored in subfolders of DISTRIBUTE/SRC.
4) Documentation of the top level functions can be found in DOC/
5) User Guide is found [here](https://github.com/jingrebeccali/SpinDoctor/blob/master/User%20Guide.pdf)

Authors: Jing-Rebecca Li, Van-Dang Nguyen, Try Nguyen Tran.

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
