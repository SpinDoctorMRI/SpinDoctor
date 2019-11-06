# SpinDoctor MatrixFormalismModule

SpinDoctor is a software package that performs numerical simulations of diffusion magnetic resonance imaging (dMRI) for prototyping purposes.  

The MatrixFormalismModule computes a closed form representation of the diffusion MRI signal called Matrix Formalism, that is based on the eigendecomposition of the Laplace operator, defined in the diffusion geometry.  Currently, it computes this representations in neurons subject to impermeable boundary conditions, for the PGSE sequence.  Permeable membranes and general diffusion-encoding sequences are subjects of future work.

	The stable version of the MatrixFormalismModule will be released on November 11. 
	Do not download the MatrixFormalismModule until then.

Paper about SpinDoctor can be found at https://arxiv.org/abs/1902.01025

Website of SpinDoctor can be found at http://www.cmap.polytechnique.fr/~jingrebeccali/software.html

Software requirements: SpinDoctor MatrixFormalismModule is compatible with MATLAB (version R2017b or later).

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
