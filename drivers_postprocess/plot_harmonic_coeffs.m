function plot_harmonic_coeffs(a,L_max,title_str)


m = -L_max:L_max;
l = 0:L_max;
imagesc(m,l,abs(a));
colorbar;
if nargin == 3
title(sprintf('Harmonic coefficients for %s',title_str));
else
title('Harmonic coefficients');
end
