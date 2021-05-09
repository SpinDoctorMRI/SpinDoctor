function nu = compute_nu(n, d)
%COMPUTE_NU Compute the angular eigenvalue nu as a function of the angular index n.

switch d
  case 2
    nu = n^2;
  case 3
    nu = n * (n + 1); 
end





