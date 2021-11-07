function fullstring = gradient_fieldstring(ug, b)
%GRADIENT_FIELDSTRING Get fieldname string representation of gradient.
%
%   ug: [3 x 1] gradient direction
%   b: double   gradient amplitude
%
%   fullstring: string


if nargin == 2
    if abs(b) < 1e-14
        fullstring = sprintf("b0");
    else
        fullstring = sprintf("b%g_ug%.4f_%.4f_%.4f", b, ug);
    end
elseif nargin == 1
    fullstring = sprintf("ug%.4f_%.4f_%.4f", ug);
else
    error('Illegal input.')
end

% replace illegal character
% replace decimal point '.' with 'o'
fullstring = replace(fullstring, '.', 'o');
% replace minus sign '-' with 'n'
fullstring = replace(fullstring, '-', 'n');
end
