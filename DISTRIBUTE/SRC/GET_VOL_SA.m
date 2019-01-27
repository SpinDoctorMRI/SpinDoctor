function [VOL_cmpts,SA_cmpts,SAu_cmpts,VOL_allcmpts,VF_cmpts,SoV_cmpts] ...
    = GET_VOL_SA(mymesh,gdir)

% compute surface to volume ratio
% 
% Input:
%     1. mymesh is a structure with 10 elements:
%         Nnode
%         Nele
%         Nface
%         Pts_cmpt_reorder
%         Ele_cmpt_reorder
%         Pts_ind
%         Pts_boundary_reorder
%         Fac_boundary_reorder
%         Nboundary
%         Ncmpt
%     2. gdir
%     
% Output: 
%     1. VOL_cmpts
%     2. SA_cmpts
%     3. SAu_cmpts
%     4. VOL_allcmpts
%     5. VF_cmpts
%     6. SoV_cmpts

UG = gdir';
UG = UG/norm(UG);

for icmpt = 1:mymesh.Ncmpt
    Fac = [];
    for iboundary = 1:mymesh.Nboundary
        Fac = [Fac,mymesh.Fac_boundary_reorder{icmpt}{iboundary}];
    end
    [VOL_cmpts(icmpt)] ...
        = get_volume_mesh(mymesh.Pts_cmpt_reorder{icmpt},mymesh.Ele_cmpt_reorder{icmpt});
    [SA_cmpts(icmpt),SAu_cmpts(icmpt,:)] ...
        = get_surface_mesh(mymesh.Pts_cmpt_reorder{icmpt},Fac,UG);
end

VOL_allcmpts = 0;

for icmpt = 1:mymesh.Ncmpt
    VOL_allcmpts  = VOL_allcmpts + VOL_cmpts(icmpt);
end

for icmpt = 1:mymesh.Ncmpt
    VF_cmpts(icmpt) = VOL_cmpts(icmpt)/VOL_allcmpts;
end

for icmpt = 1:mymesh.Ncmpt
    SoV_cmpts(icmpt,:) = SAu_cmpts(icmpt,:)/VOL_cmpts(icmpt);
end