function [in, on, outer, out_near] = create_mask(mask,dist,r_min,eps)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin == 3
    eps = 1e-8;
end

in = mask & (dist < -eps);
on = mask & (-eps <=dist ) & (dist <= eps);
outer = mask & (dist > eps);
out_near = mask & ( dist > eps) & ( dist< .01*r_min);

end

