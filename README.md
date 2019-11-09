# SpinDoctor Matrix Formalism Module

SpinDoctor is a software package that performs numerical simulations of diffusion magnetic resonance imaging (dMRI) for prototyping purposes.  

The Matrix Formalism Module computes a closed form representation of the diffusion MRI signal called Matrix Formalism, that is based on the eigendecomposition of the Laplace operator, defined in the diffusion geometry.  

Currently, the Matrix Formalism Module allows the computation of the Matrix Formalism signal and the 
Matrix Formalism Gaussian Approximation signal for realistic neuron (impermeable membranes) with the PGSE sequence.

Matrix Formalism for permeable membranes and for general diffusion-encoding sequences are under development 
and will be released in the future.  

	The stable version of the MatrixFormalismModule will be released on November 11. 
	Do not download the MatrixFormalismModule until then.

Software requirements
 
	The SpinDoctor toolbox and the Neuron Module have been developed in the MATLAB R2017b and require no additional MATLAB toolboxes.  However, the current version of the Matrix Formalism Module requires the MATLAB PDE Toolbox (2017 or later) 
	due to certain difficulties of implementing the matrix eigenvalue solution on a restricted eigenvalue interval.  
	This technical issue will be addressed in a future release.      

Getting started
	
    The DISTRIBUTE folder contains a general purpose driver.
    Two sets of input files are provided, one for simulations in one gradient direction, and one for HARDI simulations.
    The functions that the user is likely to call directly from the driver are located at the top level of DISTRIBUTE/SRC.
    Other functions are stored in subfolders of DISTRIBUTE/SRC.
    Documentation of the top level functions can be found in DOC/

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
