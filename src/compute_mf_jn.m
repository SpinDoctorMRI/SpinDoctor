function mf_jn = compute_mf_jn(lap_eig, setup)
%COMPUTE_MF_JN Compute the quantity J(lambda_n, f).
%
%   lap_eig: struct
%   setup: struct
%
%   mf_jn


disp("Computing the quantity J(lambda_n, f)");

% Extract experiment parameters
sequences = setup.gradient.sequences;
ninterval = setup.mf.ninterval;

% Sizes
nsequence = length(sequences);
n_lap_eig = length(lap_eig);

% Initialize output argument
mf_jn = cell(1, n_lap_eig);

for ilapeig = 1:n_lap_eig
    eigvals = lap_eig(ilapeig).values;
    neig = length(eigvals);
    mf_jn_ilapeig = zeros(nsequence, neig);

    % Compute Jn
    for iseq = 1:nsequence
        if n_lap_eig > 1
            fprintf(" Compartment %d: experiment %d of %d\n", ilapeig, iseq, nsequence);
        else
            fprintf(" Experiment %d of %d\n", iseq, nsequence);
        end

        % Gradient sequences
        seq = sequences{iseq};

        for ieig = 1:neig
            lambda = eigvals(ieig);
            if abs(lambda) < 1e-16
                mf_jn_ilapeig(iseq, ieig) = 0;
            elseif ismember(class(seq), ["PGSE", "DoublePGSE", "CosOGSE", "SinOGSE"])
                mf_jn_ilapeig(iseq, ieig) = seq.J(lambda);
            else
                mf_jn_ilapeig(iseq, ieig) = seq.J(lambda, ninterval);
            end
        end
    end
    mf_jn{ilapeig} = mf_jn_ilapeig;
end
