function [points_gdir,graddir_index,negii] = HARDI_PTS(ngdir_total)
% compute the diffusion-encoding directions uniformly distributed on the sphere.
[points_gdir,C,v] = spheresurface_regularpoints(1,ngdir_total);
ngdir_total = size(points_gdir,1);
% ii = find(points_gdir(:,3) >= 0);
% % negii
% for j = 1:length(ii)
%     for k = 1:ngdir_total
%         if (norm(points_gdir(ii(j),1:2)+points_gdir(k,1:2)) < 1e-10  ...
%                 && points_gdir(ii(j),3)+points_gdir(k,3) < 1e-10)
%             negii{ii(j)} = k;
%         end
%     end
% end
% graddir_index = ii;
graddir_index = 1:ngdir_total;
ndir = length(graddir_index);
clear negii;
for idir = 1:ndir
    negii{idir} = [];
end
