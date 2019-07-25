# SpinDoctor NeuronModule

SpinDoctor is a software package that performs numerical simulations of diffusion magnetic resonance imaging (dMRI) for prototyping purposes.  

The NeuronModule solves the Bloch-Torrey PDE on neurons.  The PDEs are solved by P1 finite elements combined with built-in MATLAB routines for solving ordinary differential equations.  The finite element mesh generation is performed using an external package called TetGen that is included in the toolbox.

Paper about SpinDoctor can be found at https://arxiv.org/abs/1902.01025
Website of SpinDoctor can be found at http://www.cmap.polytechnique.fr/~jingrebeccali/software.html

Software requirements: SpinDoctor NeuronModule is compatible with MATLAB (version R2017b or later).

Getting started
1) The DISTRIBUTE folder contains a general purpose driver called driver_spindoctor_btpde_neuron_clean.m 
2) Two sets of input files are provided, one for simulations in one gradient direction, and one for HARDI simulations.
3) The functions that the user is likely to call directly from the driver are located at the top level of DISTRIBUTE/SRC.
4) Other functions are stored in subfolders of DISTRIBUTE/SRC.
5) Documentation of the top level functions can be found in DOC/

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
