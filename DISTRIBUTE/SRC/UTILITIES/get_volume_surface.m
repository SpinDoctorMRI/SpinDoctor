function [VOL,SA,SAu,VOL_allcmpts,VOL_frac,SoV] ...
    = get_volume_surface(mymesh,gdir)

UG = gdir';
UG = UG/norm(UG);

for icmpt = 1:mymesh.Ncmpt
    Fac = [];
    for iboundary = 1:mymesh.Nboundary
        Fac = [Fac,mymesh.Fac_boundary_reorder{icmpt}{iboundary}];
    end
    [VOL(icmpt)] ...
        = get_volume_mesh(mymesh.Pts_cmpt_reorder{icmpt},mymesh.Ele_cmpt_reorder{icmpt});
    [SA(icmpt),SAu(icmpt,:)] ...
        = get_surface_mesh(mymesh.Pts_cmpt_reorder{icmpt},Fac,UG);
end

VOL_allcmpts = 0;

for icmpt = 1:mymesh.Ncmpt
    VOL_allcmpts  = VOL_allcmpts + VOL(icmpt);
end

for icmpt = 1:mymesh.Ncmpt
    VOL_frac(icmpt) = VOL(icmpt)/VOL_allcmpts;
end

for icmpt = 1:mymesh.Ncmpt
    SoV(icmpt,:) = SAu(icmpt,:)/VOL(icmpt);
end