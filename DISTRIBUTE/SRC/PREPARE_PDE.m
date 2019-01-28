function [DIFF_cmpts,kappa_bdys,IC_cmpts,OUT_cmpts_index,ECS_cmpts_index,IN_cmpts_index,Ncmpt,Nboundary] ...
    = PREPARE_PDE(ncell,cell_shape,params_domain_geom,params_domain_pde)

% create initial data for solving PDE
% 
% Input:
%     1. ncell
%     2. cell_shape
%     3. params_domain_geom is a structure with 3 elements:
%         Rratio_IN
%         include_ECS
%         ECS_gap    
%     4. params_domain_pde is a structure with 8 elements:
%         dcoeff_IN
%         dcoeff_OUT
%         dcoeff_ECS
%         ic_IN
%         ic_OUT
%         ic_ECS
%         kappa_IN_OUT
%         kappa_OUT_ECS
% 
% Output:
%     1. DIFF_cmpts
%     2. kappa_bdys
%     3. IC_cmpts
%     4. OUT_cmpts_index
%     5. ECS_cmpts_index
%     6. IN_cmpts_index
%     7. Ncmpt
%     8. Nboundary

Rratio_IN = params_domain_geom.Rratio_IN;
include_ECS = params_domain_geom.include_ECS;

if (Rratio_IN <= 0 || Rratio_IN >= 1)
    Ncmpt = ncell;
else
    Ncmpt = 2*ncell;
end

if (include_ECS ~= 0)
    Ncmpt = Ncmpt + 1;
end

if (Rratio_IN <= 0 || Rratio_IN >= 1)
    if (cell_shape == 2)
        % for a cylinder, the side is one boundary, the top and bottom another boundary
        Nboundary = 2*ncell;
    elseif (cell_shape == 1)
        % for a sphere, there is one boundary
        Nboundary = ncell;
    else
        disp('wrong');
        stop
    end
else
    if (cell_shape == 2)
        % for a cylinder, there is the axon side wall, the axon top and bottom,
        % the myelin side wall, the myelin top and bottom.
        Nboundary = 4*ncell;
    elseif (cell_shape == 1)
        % for a sphere, there is the outer sphere and the inner sphere
        Nboundary = 2*ncell;
    else
        disp('wrong');
        stop
    end
end

if (include_ECS ~= 0)
    Nboundary = Nboundary+1;
end

DIFF_cmpts = zeros(1,Ncmpt);
kappa_bdys = zeros(1,Nboundary);

if (include_ECS ~= 0)
    % the ECS is the last cmpt
    % the ECS outer boundary is the last boundary
    ECS_cmpts_index = Ncmpt;
    Box_boundary = Nboundary;
else
    ECS_cmpts_index = [];
    Box_boundary = [];
end

DIFF_cmpts(1,ECS_cmpts_index) = params_domain_pde.dcoeff_ECS;
% the outer boundary of ECS is always impermeable
kappa_bdys(1,Box_boundary) = 0;
IC_cmpts(1,ECS_cmpts_index) = params_domain_pde.ic_ECS;

if (cell_shape == 2)
    if (Rratio_IN <= 0 || Rratio_IN >= 1)
        OUT_cmpts_index = 1:ncell;
        IC_cmpts(1,OUT_cmpts_index) = params_domain_pde.ic_OUT;
        DIFF_cmpts(1,OUT_cmpts_index) = params_domain_pde.dcoeff_OUT;
        OUT_ECS_boundary = 1:2:2*ncell;
		if (include_ECS ~= 0)
            kappa_bdys(1,OUT_ECS_boundary) = params_domain_pde.kappa_OUT_ECS;
		end 
		IN_cmpts_index = [];        
    else
        IN_cmpts_index = 1:ncell;
        IC_cmpts(1,IN_cmpts_index) = params_domain_pde.ic_IN;
        DIFF_cmpts(1,IN_cmpts_index) = params_domain_pde.dcoeff_IN;
        OUT_cmpts_index = ncell+1:2*ncell;
        DIFF_cmpts(1,OUT_cmpts_index) = params_domain_pde.dcoeff_OUT;
        IC_cmpts(1,OUT_cmpts_index) = params_domain_pde.ic_OUT;        
        IN_OUT_boundary = [1:4:4*ncell];
        kappa_bdys(1,IN_OUT_boundary) = params_domain_pde.kappa_IN_OUT;
        OUT_ECS_boundary = [3:4:4*ncell];
        if (include_ECS ~= 0)
            kappa_bdys(1,OUT_ECS_boundary) = params_domain_pde.kappa_OUT_ECS;    
        end
    end
    
elseif (cell_shape == 1)
    OUT_cmpts_index = 1:ncell;
    IC_cmpts(1,OUT_cmpts_index) = params_domain_pde.ic_OUT;
    DIFF_cmpts(1,OUT_cmpts_index) = params_domain_pde.dcoeff_OUT;
    if (Rratio_IN <= 0 || Rratio_IN >= 1)        
        OUT_ECS_boundary = [1:ncell];
        if (include_ECS ~= 0)
            kappa_bdys(1,OUT_ECS_boundary) = params_domain_pde.kappa_OUT_ECS;
        end
        IN_cmpts_index = [];        
    else        
        IN_cmpts_index = ncell+1:2*ncell;
        DIFF_cmpts(1,IN_cmpts_index) = params_domain_pde.dcoeff_IN;
        IC_cmpts(1,IN_cmpts_index) = params_domain_pde.ic_IN;        
        IN_OUT_boundary = 2:2:2*ncell;
        kappa_bdys(1,IN_OUT_boundary) = params_domain_pde.kappa_IN_OUT;
        OUT_ECS_boundary = 1:2:2*ncell;
        if (include_ECS ~= 0)
            kappa_bdys(1,OUT_ECS_boundary) = params_domain_pde.kappa_OUT_ECS;  
        end
    end
end

% if (Rratio_IN <= 0)
%     if (cell_shape == 2)
%         % this is the axon
%         IN_ECS_boundary = 1:2:2*ncell;
%         kappa_bdys(1,IN_ECS_boundary) = params_domain_pde.kappa_IN_ECS;
%         OUT_cmpts_index = [];
%     elseif (cell_shape == 1)
%         OUT_ECS_boundary = [1:ncell];
%         kappa_bdys(1,OUT_ECS_boundary) = params_domain_pde.kappa_OUT_ECS;
%         IN_cmpts_index = [];
%     else
%         disp('wrong');
%         stop
%     end
% else
%     if (cell_shape == 2)
%         OUT_cmpt_index = ncell+1:2*ncell;
%         DIFF_cmpts(1,OUT_cmpt_index) = params_domain_pde.dcoeff_OUT;
%         IC_cmpts(1,OUT_cmpt_index) = params_domain_pde.ic_OUT;
%         
%         IN_OUT_boundary = [1:4:4*ncell];
%         kappa_bdys(1,IN_OUT_boundary) = params_domain_pde.kappa_IN_OUT;
%         OUT_ECS_boundary = [3:4:4*ncell];
%         kappa_bdys(1,OUT_ECS_boundary) = params_domain_pde.kappa_OUT_ECS;
%         
%         
%     elseif (cell_shape == 1)
%         IN_cmpts_index = ncell+1:2*ncell;
%         DIFF_cmpts(1,IN_cmpts_index) = params_domain_pde.dcoeff_IN;
%         IC_cmpts(1,IN_cmpts_index) = params_domain_pde.ic_IN;
%         
%         IN_OUT_boundary = 2:2:2*ncell;
%         kappa_bdys(1,IN_OUT_boundary) = params_domain_pde.kappa_IN_OUT;
%         OUT_ECS_boundary = 1:2:2*ncell;
%         kappa_bdys(1,OUT_ECS_boundary) = params_domain_pde.kappa_OUT_ECS;
%     else
%         disp('wrong');
%         stop
%     end
% end