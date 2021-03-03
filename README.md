# SpinDoctor Toolbox

SpinDoctor is a software package that performs numerical simulations of diffusion magnetic resonance imaging (dMRI) for prototyping purposes.

SpinDoctor can be used 

1) to solve the Bloch-Torrey partial differential equation (BTDPE) to obtain the dMRI signal (the toolbox provides a way of robustly fitting the dMRI signal to obtain the fitted Apparent Diffusion Coefficient (ADC));
2) to solve the diffusion equation for the homogenized ADC (HADC) model to obtain the ADC;
3) a short-time approximation formula for the ADC is also included in the toolbox for comparison with the simulated ADC;
4) Compute the dMRI signal using a matrix formalism (MF) analytical solution based Laplace eigenfunctions.

The PDEs and Laplace eigenvalue decompositions are solved by P1 finite elements combined with built-in MATLAB routines for solving ordinary differential equations.
The finite element mesh generation is performed using an external package called TetGen that is included in the toolbox.

SpinDoctor has support for the following features:
1. multiple compartments with compartment-wise constant
	* initial spin densities,
	* diffusion coefficients or diffusion tensors, and
	* T2-relaxation coefficients;
2. permeable membranes between compartments for the BTPDE and MF (the HADC assumes negligible permeability);
3. built-in diffusion-encoding pulse sequences, including
	* the Pulsed Gradient Spin Echo (PGSE) and double-PGSE,
	* the Ocsillating Gradient Spin Echo (OGSE), or
	* custom pulse sequences;
4. uniformly distributed gradient directions in 2D and 3D for high angular resolution diffusion imaging (HARDI)

SpinDoctor also comes with a geometry generation module, allowing for
1. spherical cells with a nucleus;
2. cylindrical cells with a myelin layer;
3. an extra-cellular space (ECS) enclosed in either
	* a box,
	* a convex hull, or
	* a tight wrapping around the cells;
4. deformation of canonical cells by bending and twisting.

In addition, a variety of neuron meshes is available, whose surface geometries were extracted from [NeuroMopho.org](http://neuromorpho.org). The neurons may also be enclosed in an extracellular space as described above.


## SpinDoctor article
The paper about SpinDoctor can be found at https://arxiv.org/abs/1902.01025.

## Software requirements

	The SpinDoctor Toolbox has been developed in the MATLAB R2020b
	and require no additional MATLAB Toolboxes.

	However, if the MATLAB Parallel Computing Toolbox is available,
	the simulations can be run in parallel.
 
 
## Getting started

1) The base folder contains a commented general purpose driver called `driver_spindoctor.m`. The other driver, `driver_save_load.m`, can save and load simulations.
2) The input files for the drivers are found in the folder `input_files`, and define the structures needed for the simulations.
3) Multiple neuron meshes are found in the folder `mesh_files`. These can be loaded in the `input_files/define_params_neuron.m` script.
4) The user guide is found [here](https://github.com/jingrebeccali/SpinDoctor/blob/master/user_guide.pdf).

Authors: Jing-Rebecca Li, Syver Døving Agdestein, Chengran Fang, Van-Dang Nguyen, Try Nguyen Tran.


## How to cite us

If you use our software for research, please consider citing us:
```bibtex
@article{Li2019,
title = {{SpinDoctor: A MATLAB toolbox for diffusion MRI simulation}},
journal = {NeuroImage},
volume = {202},
pages = {116120},
year = {2019},
issn = {1053-8119},
doi = {https://doi.org/10.1016/j.neuroimage.2019.116120},
url = {http://www.sciencedirect.com/science/article/pii/S1053811919307116},
author = {Jing-Rebecca Li and Van-Dang Nguyen and Try Nguyen Tran and Jan Valdman and Cong-Bang Trang and Khieu Van Nguyen and Duc Thach Son Vu and Hoang An Tran and Hoang Trong An Tran and Thi Minh Phuong Nguyen},
keywords = {Bloch-torrey equation, Diffusion magnetic resonance imaging, Finite elements, Simulation, Apparent diffusion coefficient}
}
```

Citations for the Neuron Module and Matrix Formalism module can be found in the `CITATION.bib` file.



## License

	Copyright (C) 2019-2021 Jing-Rebecca Li

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
