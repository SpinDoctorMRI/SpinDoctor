function magnetization_avg = average_magnetization(magnetization)
%AVERAGE_MAGNETIZATION Compute direction averaged magnetization
%
%   magnetization: cell(ncompartment, namplitude, nsequence, ndirection)[nnode]
%
%   magnetization_avg: cell(ncompartment, namplitude, nsequence)[nnode]


% Sizes
[ncompartment, namplitude, nsequence, ~] = size(magnetization);

% Average over directions
magnetization_avg = cell(ncompartment, namplitude, nsequence);
for iseq = 1:nsequence
    for iamp = 1:namplitude
        for icmpt = 1:ncompartment
            magnetization_avg{icmpt, iamp, iseq} = mean([magnetization{icmpt, iamp, iseq, :}], 2);
        end
    end
end
