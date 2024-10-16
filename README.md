# Neuron Meshing paper

This code was used in a forthcoming publication by us to simulate diffusion magnetic resonance imaging (dMRI) experiments on realistic neuron and microlgia geometries. The primary purpose of the code is to test our meshes against the meshes produced by the Ultraliser tool (available here https://github.com/BlueBrain/Ultraliser).

To test for an individual cell, run 
```
driver_demo;
```
The sequences, cell geometry and experimental parameters can be adjusted in the setup file.

The simulations for the neuron 1-2-2.CNG can be done by running:
```
driver_run_neuron_simulations;
```

The simulations for the human microglia cells can be done by running:
```
driver_run_microglia_simulations; 
```
and comparing the signals of individual compartments is done in:
```
compare_microlgia_segmented;
```

We note that if the tetgen parameters are changed then the correspondance between branches may be altered. If branches are split/merged then it is necessary to manually assign the correspondance between branches. This is done in the switch statements in compare_microglia_segmented and also in the last section of driver_run_neuron_simulations.

The figures for the paper were plotted using plot_data.ipynb.


# SpinDoctor Toolbox

SpinDoctor is a software package that performs numerical simulations of diffusion magnetic resonance imaging (dMRI) for prototyping purposes. This branch contains the code and data used in a forthcoming publication describing our surface meshing algorithm https://github.com/AlexMcSD/mesh_swc.

## Software requirements

The SpinDoctor Toolbox has been developed in the MATLAB R2020b and tested with MATLAB R2018a-R2021a.

SpinDoctor requires no additional MATLAB Toolboxes. However, if the MATLAB Parallel Computing Toolbox is available,
the simulations can be run in parallel.


## Updates July 2024

New features:

* Simulation of general gradient sequences is now supported. These sequences can be read directly from a Camino .scheme file. Simulations can be run for solve_btpde, solve_btpde_midpoint, and solve_mf.
* Segmenting a finite element mesh into soma and processes/dendrites when the mesh comes from a swc file.
* Code refactoring, to incorporate the new features.



## Updates March 2021

All modules (master, T2, NeuronModule, MatrixFormalismModule) are now merged and found in the `master` branch. The previous versions of the modules are still found in their branches.

New features:

* Code refactoring, change of data structures
* New input file system, controlled by a structure `setup`
* Some speed and memory improvements (extensive vectorization, reduced memory usage in ODE solvers)
* Support for parallel outer loops (solve for multiple directions, amplitudes and sequences in parallel)
* Diffusion tensors as an alternative to scalar diffusion coefficients
* Matrix formalism with built in Matlab functions instead of PDE Toolbox
* Add convex hull ECS + tight wrap ECS for neurons (experimental feature)
* Simple interface for pulse sequences, the user can add custom sequences
	* by passing a function handle, for which the integral quantities are computed numerically,
	* or by manually defining a subclass of the `Sequence` class


## Authors

Jing-Rebecca Li, Syver Døving Agdestein, Chengran Fang, Van-Dang Nguyen, Try Nguyen Tran.

## How to cite us

The original paper about SpinDoctor can be found at https://arxiv.org/abs/1902.01025.

If you use our software for research, please cite us:

```bibtex
@article{Li2019,
  author  = {Jing-Rebecca Li and Van-Dang Nguyen and Try Nguyen Tran and Jan Valdman and Cong-Bang Trang and Khieu Van Nguyen and Duc Thach Son Vu and Hoang An Tran and Hoang Trong An Tran and Thi Minh Phuong Nguyen},
  doi     = {https://doi.org/10.1016/j.neuroimage.2019.116120},
  issn    = {1053-8119},
  journal = {NeuroImage},
  pages   = {116120},
  title   = {{SpinDoctor: A MATLAB toolbox for diffusion MRI simulation}},
  url     = {http://www.sciencedirect.com/science/article/pii/S1053811919307116},
  volume  = {202},
  year    = {2019}
}
```

Citations for further developements such as the Neuron Module and the Matrix Formalism module can be found in the `CITATION.bib` file.



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
