# SpinDoctor

SpinDoctor is a software package that performs numerical simulations for diffusion magnetic resonance imaging (dMRI) for prototyping purposes.

SpinDoctor can be used 

1) to solve the Bloch-Torrey PDE to obtain the dMRI signal 
(the toolbox provides a way of robustly fitting the dMRI signal to obtain the fitted Apparent Diffusion Coefficient (ADC)); 
2) to solve the diffusion equation (DE) of the H-ADC model to obtain the ADC;
3) a short-time approximation formula for the ADC is also included in the toolbox for comparison with the simulated ADC.

The PDEs are solved by P1 finite elements combined with build-in Matlab routines for solving ordinary differential equations.
The finite element mesh generation is performed using an external package called Tetgen that is included in the toolbox.

SpinDoctor provides built-in options of including 
1) spherical cells with a nucleus; 
2) cylindrical cells with a myelin layer; 
3) an extra-cellular space (ECS) enclosed either 
a) in a box or b) in a tight wrapping around the cells; 
4) deformation of canonical cells by bending and twisting.  
5) permeable membranes for the BT-PDE (the H-ADC assumes negligible permeabilty).
Built-in diffusion-encoding pulse sequences include the Pulsed Gradient Spin Echo and the Ocsilating Gradient Spin Echo.    

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
