function [VOL_cmpts,SA_cmpts,SAu_cmpts,VOL_allcmpts,VF_cmpts,SoV_cmpts] ...
    = GET_VOL_SA(mymesh,gdir)

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