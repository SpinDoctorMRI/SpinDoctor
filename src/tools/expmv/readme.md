`expmv` - Matrix exponential times a vector.
==========

About
-----

`expmv` contains two MATLAB functions for computing expm(t\*A)\*b without
explicitly forming expm(t\*A), where A is an n-by-n matrix and b is an
n-by-1 vector. This is the problem of computing the action of the matrix
exponential on a vector.

expmv(t,A,B) computes expm(t\*A)\*B, while expmv_tspan(A,b,t0,tmax,q)
computes expm(t\*A)*b for q+1 >= 2 equally spaced values of t between t0 and
tmax.

The functions work for any matrix A, and use just matrix-vector products
with A and A^*.

Function test.m runs a simple test of the codes.

Details on the underlying algorithms can be found in

A. H. Al-Mohy and N. J. Higham, "[Computing the action of the matrix
exponential, with an application to exponential
integrators](https://doi.org/10.1137/100788860)" SIAM
J. Sci. Comput., 33(2):488--511, 2011.
