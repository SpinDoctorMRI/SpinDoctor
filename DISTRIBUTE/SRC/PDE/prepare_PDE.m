function [DIFF_cmpts,kappa_vec,IC_cmpts,Cell_cmpt,Box_cmpt,Nucleus_cmpt,Ncmpt,Nboundary] ...
    = prepare_PDE(ncell,cell_shape,domain_geom,domain_pde)
	
Rratio_nucleus = domain_geom.Rratio_nucleus;
include_box = domain_geom.include_box;	

if (Rratio_nucleus <= 0 | Rratio_nucleus >= 1)
	Ncmpt = ncell;
else
	Ncmpt = 2*ncell;
end

if (include_box ~= 0)
	Ncmpt = Ncmpt + 1;
end

if (Rratio_nucleus <= 0)
	if (cell_shape == 2)
		Nboundary = 2*ncell;
	elseif (cell_shape == 1)
		Nboundary = ncell;
	else
		disp('wrong');
		stop
	end
else
    if (cell_shape == 2)
		Nboundary = 4*ncell;
    elseif (cell_shape == 1)
		Nboundary = 2*ncell;
    else
        disp('wrong');
        stop
    end
end
if (include_box ~= 0)
    Nboundary = Nboundary+1;
end 
	  
DIFF_cmpts = zeros(1,Ncmpt);
kappa_vec = zeros(1,Nboundary);

if (include_box ~= 0)
    Box_cmpt = Ncmpt;
    Box_boundary = Nboundary;
else
    Box_cmpt = [];
    Box_boundary = [];
end

DIFF_cmpts(1,Box_cmpt) = domain_pde.dcoeff_exterior;
kappa_vec(1,Box_boundary) = 0;
IC_cmpts(1,Box_cmpt) = domain_pde.ic_exterior;

if (cell_shape == 2)
    Cell_cmpt = 1:ncell;
elseif (cell_shape == 1)
    Cell_cmpt = 1:ncell;
end

IC_cmpts(1,Cell_cmpt) = domain_pde.ic_cytoplasm;
DIFF_cmpts(1,Cell_cmpt) = domain_pde.dcoeff_cytoplasm;

if (Rratio_nucleus <= 0)
	if (cell_shape == 2)
		CE_boundary = 1:2:2*ncell;
		kappa_vec(1,CE_boundary) = domain_pde.kappa_ce;
		Nucleus_cmpt = [];
	elseif (cell_shape == 1)
		CE_boundary = [1:ncell];
		kappa_vec(1,CE_boundary) = domain_pde.kappa_ce;
		Nucleus_cmpt = [];
	else
        disp('wrong');
        stop
    end
else
    if (cell_shape == 2)
        Myelin_cmpt = ncell+1:2*ncell;
        DIFF_cmpts(1,Myelin_cmpt) = domain_pde.dcoeff_nucleus;
        CM_boundary = [1:4:4*ncell];
        kappa_vec(1,CM_boundary) = domain_pde.kappa_nc;
        ME_boundary = [3:4:4*ncell];
        kappa_vec(1,ME_boundary) = domain_pde.kappa_ce;
        IC_cmpts(1,Myelin_cmpt) = domain_pde.ic_nucleus;
        Nucleus_cmpt = Myelin_cmpt;
    elseif (cell_shape == 1)
        Nucleus_cmpt = ncell+1:2*ncell;
        DIFF_cmpts(1,Nucleus_cmpt) = domain_pde.dcoeff_nucleus;
        NC_boundary = 2:2:2*ncell;
        kappa_vec(1,NC_boundary) = domain_pde.kappa_nc;
        CE_boundary = 1:2:2*ncell;
        kappa_vec(1,CE_boundary) = domain_pde.kappa_ce;
        IC_cmpts(1,Nucleus_cmpt) = domain_pde.ic_nucleus;
    else
        disp('wrong');
        stop
    end
end