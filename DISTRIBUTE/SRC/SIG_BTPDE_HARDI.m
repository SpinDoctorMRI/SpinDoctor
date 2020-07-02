function [SIG_BTPDE_cmpts_alldir,SIG_BTPDE_allcmpts_alldir,ctime_alldir] ...
    = SIG_BTPDE_HARDI(experi_btpde,mymesh,DIFF_cmpts,kappa_bdys,IC_cmpts,...
    points_gdir)

% compute the signals from Bloch-Torrey equation for ngdir directions.
% 
% Input:
%     1. experiment_btpde is a structure with 10 elements:
%         ngdir_total 
%         gdir        
%         sdeltavec   
%         bdeltavec   
%         seqvec       
%         npervec     
%         rtol        
%         atol        
%         qvalues     
%         bvalues        
%     2. mymesh is a structure with 10 elements:
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
%     3. DIFF_cmpts
%     4. kappa_bdys
%     5. IC_cmpts
%     6. points_gdir
%     7. graddir_index
%     8. negii
% 
% Output:
%     1. SIG_BTPDE_cmpts_alldir
%     2. SIG_BTPDE_allcmpts_alldir
%     3. ctime_alldir

Ncmpt = length(DIFF_cmpts);

nexperi = length(experi_btpde.sdeltavec);
nb = 1;
ngdir_total = size(points_gdir,1);

SIG_BTPDE_cmpts_alldir = nan*ones(ngdir_total,Ncmpt,nexperi,nb);
SIG_BTPDE_allcmpts_alldir = nan*ones(ngdir_total,nexperi,nb);
ctime_alldir = nan*ones(ngdir_total,nexperi,nb);

OUTPUT_MAGNETIZATION = false;
for idir = 1:ngdir_total
    experi_btpde.gdir = points_gdir(idir,:); 
    [~,SIG_cmpts,SIG_allcmpts,~,ctime] ...
        = SIG_BTPDE(experi_btpde,mymesh,DIFF_cmpts,kappa_bdys,IC_cmpts,OUTPUT_MAGNETIZATION);      
    ctime_alldir(idir,:,:) = ctime;    
    SIG_BTPDE_cmpts_alldir(idir,:,:,:) = SIG_cmpts;
    SIG_BTPDE_allcmpts_alldir(idir,:,:) = SIG_allcmpts(:,:);
end
