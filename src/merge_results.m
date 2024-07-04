function results = merge_results(camino,const,nsequence_camino,nsequence_const,totaltime,save_magnetization)
%MERGE_RESULTS Takes output of experiments from direction varying sequences
%and constant direction vector sequences and creates a singles results
%struct
% 
%   camino: struct
%   const:struct
%   nsequence_camino: int
%   nsequecne_const: int
% 
%   results: struct



% Create output structure
results.totaltime = totaltime;
if nsequence_camino == 0
    results.signal = const.signal;
    results.signal_allcmpts = const.signal_allcmpts;
    results.itertimes = const.itertimes;
    if save_magnetization
        results.magnetization = const.magnetization;
        results.magnetization_avg = average_magnetization(const.magnetization);
    end
    % results.signal_weighted = const.signal_weighted;
    % results.signal_allcmpts_weighted = const.signal_allcmpts_weighted;

elseif nsequence_const == 0
    results.signal = camino.signal;
    results.signal_allcmpts = camino.signal_allcmpts;
    results.itertimes = camino.itertimes;
    if save_magnetization
        results.magnetization = camino.magnetization;
        results.magnetization_avg = average_magnetization(camino.magnetization);
    end
    % results.signal_weighted = camino.signal_weighted;
    % results.signal_allcmpts_weighted = camino.signal_allcmpts_weighted;
else
    results.camino.signal = camino.signal;
    results.camino.signal_allcmpts = camino.signal_allcmpts;
    results.camino.itertimes = camino.itertimes;
    % results.camino.signal_weighted = camino.signal;
    % results.camino.signal_allcmpts_weighted = camino.signal_allcmpts;

    results.const.signal = const.signal;
    results.const.signal_allcmpts = const.signal_allcmpts;
    results.const.itertimes = const.itertimes;
    if save_magnetization
        results.camino.magnetization = camino.magnetization;
        results.camino.magnetization_avg = average_magnetization(camino.magnetization);
        results.const.magnetization = const.magnetization;
        results.const.magnetization_avg = average_magnetization(const.magnetization);
    end
    % results.const.signal_weighted = const.signal_weighted;
    % results.const.signal_allcmpts_weighted = const.signal_allcmpts_weighted;
end
end