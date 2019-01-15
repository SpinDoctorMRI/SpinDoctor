function [M_boundary, neumann2area]=flux_matrixP1_3D(neumann,coordinates,coeffs)
% Copyright (c) 2018, Jan Valdman

%coeffs can be only P0 (elementwise constant) function 
%represented by a collumn vector with size(elements,1) entries
%if coeffs is not provided then coeffs=1 is assumed globally
%Note: P1 coeffs needs a higher integration rule (not implemented yet)

% [I,J,V] = find(X,...) also returns a vector V containing the values
%  that correspond to the row and column indices I and J.

%computes areas on Neumann faces
neumann2area=evaluate_area(neumann,coordinates);
% nargin: Number of function input arguments.
% this will create the Q but just for the boundary nodes
if nargin==2 
    M=mass_matrixP1_2D(neumann,neumann2area);
elseif nargin==3
    M=mass_matrixP1_2D(neumann,neumann2area,coeffs);
end
% this step will create the Q for all nodes
[X,Y,Z]=find(M);
M_boundary=sparse(X,Y,Z,size(coordinates,1),size(coordinates,1));

% M1=M;M2=M_boundary;ind=find(M1>10^-13);ind1=find(M2>10^-13);
