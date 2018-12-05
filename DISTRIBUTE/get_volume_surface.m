function [VOL,SA,SAu,VOL_allcmpts,VOL_frac,SoV] ...
    = get_volume_surface(Nboundary,Ncmpt,...
    Pts_cmpt_reorder,Ele_cmpt_reorder,Fac_boundary_reorder,gdir)

UG = gdir';
UG = UG/norm(UG);

for icmpt = 1:Ncmpt
    Fac = [];
    for iboundary = 1:Nboundary
        Fac = [Fac,Fac_boundary_reorder{icmpt}{iboundary}];
    end
    [VOL(icmpt)] ...
        = get_volume_mesh(Pts_cmpt_reorder{icmpt},Ele_cmpt_reorder{icmpt});
    [SA(icmpt),SAu(icmpt,:)] ...
        = get_surface_mesh(Pts_cmpt_reorder{icmpt},Fac,UG);
end

VOL_allcmpts = 0;

for icmpt = 1:Ncmpt
    VOL_allcmpts  = VOL_allcmpts + VOL(icmpt);
end

for icmpt = 1:Ncmpt
    VOL_frac(icmpt) = VOL(icmpt)/VOL_allcmpts;
end

for icmpt = 1:Ncmpt
    SoV(icmpt,:) = SAu(icmpt,:)/VOL(icmpt);
end