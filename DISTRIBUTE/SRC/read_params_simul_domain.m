function [params_domain_geom,params_domain_pde,params_domain_femesh] ...
    = read_params_simul_domain(fname_domain)


% read simulation domain parameters
% 
% Input:
%         fname_domain
%         
% Output:
%     1. params_domain_geom is a structure with 3 elements:
%         Rratio_IN
%         include_ECS
%         ECS_gap
%         
%     2. params_domain_pde is a structure with 8 elements:
%         dcoeff_IN
%         dcoeff_OUT
%         dcoeff_ECS
%         ic_IN
%         ic_OUT
%         ic_ECS
%         kappa_IN_OUT
%         kappa_OUT_ECS
%         
%     3. params_domain_femesh is a structure with 2 elements:
%         Htetgen
%         tetgen_cmd

ndim = 3;
fid=fopen(fname_domain);

tline = fgetl(fid);
params_domain_geom.Rratio_IN = sscanf(tline,'%f',1);
  
if (params_domain_geom.Rratio_IN < 0 | params_domain_geom.Rratio_IN > 0.99)
    params_domain_geom.Rratio_IN = 0;
end

tline = fgetl(fid);
params_domain_geom.include_ECS = sscanf(tline,'%f',1);

tline = fgetl(fid);
params_domain_geom.ECS_gap = sscanf(tline,'%f',1); 

tline = fgetl(fid);
params_domain_pde.dcoeff_IN = sscanf(tline,'%f',1);

tline = fgetl(fid);
params_domain_pde.dcoeff_OUT = sscanf(tline,'%f',1);

tline = fgetl(fid);
params_domain_pde.dcoeff_ECS = sscanf(tline,'%f',1);

% initial conditions
tline = fgetl(fid);
params_domain_pde.ic_IN = sscanf(tline,'%f',1);

tline = fgetl(fid);
params_domain_pde.ic_OUT = sscanf(tline,'%f',1);

tline = fgetl(fid);
params_domain_pde.ic_ECS = sscanf(tline,'%f',1);

tline = fgetl(fid);
params_domain_pde.kappa_IN_OUT = sscanf(tline,'%f',1);
tline = fgetl(fid);
params_domain_pde.kappa_OUT_ECS = sscanf(tline,'%f',1);

tline = fgetl(fid);
params_domain_femesh.Htetgen = sscanf(tline,'%f',1);


tline = fgetl(fid);
[strpos] = regexp(tline,"'");
params_domain_femesh.tetgen_cmd = tline(strpos(1)+1:strpos(2)-1);
fclose(fid); 