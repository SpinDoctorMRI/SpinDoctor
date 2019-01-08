Compute mesh normals
====================
 
Adam H. Aitkenhead
adam.aitkenhead@physics.cr.man.ac.uk
The Christie NHS Foundation Trust
1st Decmeber 2010
 
Calculate the normal vectors for each facet of a triangular mesh.  The ordering of the vertices (clockwise/anticlockwise) is also checked for all facets if this is requested as one of the outputs. 
 
 
USAGE:
======

[coordNORMALS] = COMPUTE_mesh_normals(meshdataIN);
  ..or..
[coordNORMALS,meshdataOUT] = COMPUTE_mesh_normals(meshdataIN);
 
 
INPUT PARAMETERS
================
 
% meshdataIN   - structure   - Structure containing the faces and
%                              vertices of the mesh, in the same format
%                              as that produced by the isosurface
%                              command.
%      ..or..  - Nx3x3 array - The vertex coordinates for each facet, with:
%                                1 row for each facet
%                                3 columns for the x,y,z coordinates
%                                3 pages for the three vertices
 
 
OUTPUT PARAMETERS
=================
 
% coordNORMALS - Nx3 array   - The normal vectors for each facet, with:
%                                1 row for each facet
%                                 3 columns for the x,y,z components
%
% meshdataOUT  - (optional)  - The mesh data with the ordering of the vertices
%                              (clockwise/anticlockwise) checked.  Uses
%                              the same format as <meshdataIN>.
 
 
EXAMPLES
========

To run an example of the code:

>>  EXAMPLE_mesh_normals
 
 
NOTES
=====

- Computing <meshdataOUT> to check the ordering of the vertices in each facet may be slow for large meshes.
- Also, it may not be possible to compute <meshdataOUT> for non-manifold meshes.
 

