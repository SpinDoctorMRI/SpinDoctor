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
ncompartment = length(lap_eig);

% Initialize output argument
mf_jn = cell(1, ncompartment);

parfor icmpt = 1:ncompartment
    eigvals = lap_eig(icmpt).values;
    neig = length(eigvals);
    mf_jn_icmpt = zeros(nsequence, neig);

    % Compute Jn
    for iseq = 1:nsequence
        if ncompartment > 1
            fprintf(" Compartment %d: experiment %d of %d\n", icmpt, iseq, nsequence);
        else
            fprintf(" Experiment %d of %d\n", iseq, nsequence);
        end

        % Gradient sequences
        seq = sequences{iseq};

        for ieig = 1:neig
            lambda = eigvals(ieig);
            if abs(lambda) < 1e-16
                mf_jn_icmpt(iseq, ieig) = 0;
            elseif ismember(class(seq), ["PGSE", "DoublePGSE", "CosOGSE", "SinOGSE"])
                mf_jn_icmpt(iseq, ieig) = seq.J(lambda);
            else
                mf_jn_icmpt(iseq, ieig) = seq.J(lambda, ninterval);
            end
        end
    end
    mf_jn{icmpt} = mf_jn_icmpt;
end
