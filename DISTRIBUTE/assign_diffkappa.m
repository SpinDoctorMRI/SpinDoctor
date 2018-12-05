function [DIFF_cmpts,kappa_vec,IC_cmpts,Cell_cmpt,Box_cmpt,Nucleus_cmpt] ...
    = assign_diffkappa(Ncmpt,Nboundary,ncell,Rratio_nucleus,...
    dcoeff_nucleus,dcoeff_cytoplasm,dcoeff_exterior,...
    ic_nucleus,ic_cytoplasm,ic_exterior,kappa_nc,kappa_ce,...
    include_box,cell_shape_name)

DIFF_cmpts = zeros(1,Ncmpt);
kappa_vec = zeros(1,Nboundary);

if (include_box ~= 0)
    Box_cmpt = Ncmpt;
    Box_boundary = Nboundary;
else
    Box_cmpt = [];
    Box_boundary = [];
end
DIFF_cmpts(1,Box_cmpt) = dcoeff_exterior;
kappa_vec(1,Box_boundary) = 0;
IC_cmpts(1,Box_cmpt) = ic_exterior;

if (strcmp(cell_shape_name,'cylinders'))
    Cell_cmpt = 1:ncell;
elseif (strcmp(cell_shape_name,'ellipses'))
    Cell_cmpt = 1:ncell;
end

IC_cmpts(1,Cell_cmpt) = ic_cytoplasm;
DIFF_cmpts(1,Cell_cmpt) = dcoeff_cytoplasm;

if (Rratio_nucleus <= 0)
    CE_boundary = [1:ncell];
    kappa_vec(1,CE_boundary) = kappa_ce;
    Nucleus_cmpt = [];
else
    if (strcmp(cell_shape_name,'cylinders'))
        Myelin_cmpt = ncell+1:2*ncell;
        DIFF_cmpts(1,Myelin_cmpt) = dcoeff_nucleus;
        CM_boundary = [1:4:4*ncell];
        kappa_vec(1,CM_boundary) = kappa_nc;
        ME_boundary = [3:4:4*ncell];
        kappa_vec(1,ME_boundary) = kappa_ce;
        IC_cmpts(1,Myelin_cmpt) = ic_nucleus;
        Nucleus_cmpt = Myelin_cmpt;
    elseif (strcmp(cell_shape_name,'ellipses'))
        Nucleus_cmpt = ncell+1:2*ncell;
        DIFF_cmpts(1,Nucleus_cmpt) = dcoeff_nucleus;
        NC_boundary = 2:2:2*ncell;
        kappa_vec(1,NC_boundary) = kappa_nc;
        CE_boundary = 1:2:2*ncell;
        kappa_vec(1,CE_boundary) = kappa_ce;
        IC_cmpts(1,Nucleus_cmpt) = ic_nucleus;
    else
        disp('wrong');
        stop
    end
end