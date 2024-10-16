function [process_signals,process_signals_um] = merge_process_signals(results,results_um,femesh_neurites,femesh_neurites_um,inds,inds_um)
process_signals_copy = results.mf_neurites;
process_signals_copy_um = results_um.mf_neurites;
ngroups = length(inds);

sz = size(process_signals_copy{1}.signal);
process_signals = zeros([ngroups,sz]);
process_signals_um = zeros([ngroups,sz]);


for ii = 1:ngroups
    signals1 = zeros([length(inds{ii}),sz]);
    group1 = inds{ii};
    V = 0;
    for j = 1:length(inds{ii})
        signals1(j,:) = process_signals_copy{group1(j)}.signal;
        V = V + femesh_neurites{group1(j)}.total_volume;
    end
    process_signals(ii,:) = real(sum(signals1,1)/V);
    
    signals2 = zeros([length(inds_um{ii}),sz]);
    group2 = inds_um{ii};
    V = 0;
    for j = 1:length(inds_um{ii})
        signals2(j,:) = process_signals_copy_um{group2(j)}.signal;
        V = V + femesh_neurites_um{group2(j)}.total_volume;
    end
    process_signals_um(ii,:)= real(sum(signals2,1))/V;

end
process_signals = squeeze(process_signals); process_signals_um = squeeze(process_signals_um);