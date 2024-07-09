function nu = evolve_laplace_coef_direction_varying(nu,seq, q,moments, LT2,dtype)
%EVOLVE_LAPLACE_COEF_direction_varying Compute Laplace coefficients of
%final magnetization for a sequences with varying direction vectors.
%   nu: double(neig, 1)
%   seq: Sequence
%   moments: complex(neig, neig,3) moments in eigenfunction basis
%   LT2: double(neig, neig) or double(1,neig)
%   ninterval: int
%
%   nu: double(neig, 1)
LT2_is_matrix = size(LT2, 1) ~= 1;
if not(LT2_is_matrix)
    L = zeros([size(LT2,2),size(LT2,2)] ,func2str(dtype));
    s = size(L);
    diag_ind = 1:s(1)+1:s(1)*s(2);
    L(diag_ind) = L(diag_ind) + LT2;
    LT2 = L;
end

[time, ~, ~] = intervals(seq);
ninterval = length(time) - 1;
for i = 1:ninterval
    % Time step and time profile on given interval
    dt = time(i + 1) - time(i);
    time_midpoint = mean([time(i+1),time(i)]);
    ug = seq.vec_call(time_midpoint);
    % Obtain Bloch-Torrey operator.
    A = sum(moments .* shiftdim(ug, -2), 3);
    A = dtype(A); 
    K = LT2 + 1i*q*A;
    % Laplace coefficients of magnetization at end of interval
    nu = expm(-dt * K) * nu;
end

end