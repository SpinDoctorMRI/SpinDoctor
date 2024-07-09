function p = align_dendrites(femesh_dendrites,femesh_dendrites_r)
%%ALIGN_DENDRITES Obtains a permuation p such that femesh_dendrites_r(p) 
% corresponds to femesh_dendrites. Used to align meshes for the same cell
% at different refinement levels.
%   
%   femesh_dendrites: (N,1) cell array of finite element meshes
%   femesh_dendrites_r: (N,1) cell array of finite element meshes
% 
%   p: (N,1) array of unique intergers in range 1 to N

ndendrites = length(femesh_dendrites);
ndendrites_r = length(femesh_dendrites_r);
if ndendrites_r ~= ndendrites
    error('Refinement has split/merged dendrites');
end
centers = zeros(ndendrites,3);
for i =1:ndendrites
    [total_volume, volumes, c] = get_volume_mesh(femesh_dendrites{i}.points{1}, femesh_dendrites{i}.elements{1});
    centers(i,:) = (volumes*c')/total_volume;
end
centers_r = zeros(ndendrites,3);
for i =1:ndendrites
    [total_volume, volumes, c] = get_volume_mesh(femesh_dendrites_r{i}.points{1}, femesh_dendrites_r{i}.elements{1});
    centers_r(i,:) = (volumes*c')/total_volume;
end

P = perms(1:ndendrites);
ind = repmat(1:ndendrites,size(P,1),1);
centers_base = centers(ind,:);
centers_base = reshape(centers_base,size(P,1),ndendrites,3);
centers_new =centers_r(P,:);
centers_new = reshape(centers_new,size(P,1),ndendrites,3);
dist = vecnorm(centers_base-centers_new,2,3);
dist =  vecnorm(dist,2,2);
[M,I] = min(dist);
fprintf('Minimum distance is %f at permutation %d\n',M,I);
next_min = min(dist(1:size(P,1) ~= I));
fprintf('Next smallest = %f\n',next_min);
if next_min < 5*M
    error('Centroid method inconclusive to align dendrites')
end
p = P(I,:);

