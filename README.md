# SpinDoctor Toolbox (developmental code, not documented)

SpinDoctor is a software package that performs numerical simulations of diffusion magnetic resonance imaging (dMRI) for prototyping purposes.

This is an extended version for simulating direction varying sequences for human microglia.

## Software requirements

The SpinDoctor Toolbox has been developed in the MATLAB R2020b and tested with MATLAB R2018a-R2021a.

SpinDoctor requires no additional MATLAB Toolboxes. However, if the MATLAB Parallel Computing Toolbox is available,
the simulations can be run in parallel.


## Getting started

1) The setup scripts in the setups folder can be used to modify the solver and types of gradient sequences.
2) To accommodate batch processing, an individual cell can be run using the driver_cell.m or driver_btpde.m functions.
    Here the user can specify the tetgen options used and length scales of the matrix formalism.
    An example of how to call these functions is given in example.m
3) driver_output_signals.m will save the signals into.mat files and also to plain text files.
4) driver_get_morez_errors.m and driver_get_pgse_errors.m compares the matrix formalism signals against the finite element ode15 solution.
5) The user guide is found [here](https://github.com/jingrebeccali/SpinDoctor/blob/master/user_guide.pdf).

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
