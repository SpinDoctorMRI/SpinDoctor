function [points,SIG_EIG_cmpts_alldir,SIG_EIG_allcmpts_alldir] ...
    = MF_HARDI(experi_btpde,mymesh,IC_cmpts,EIG_value_cmpts,EIG_proj_cmpts)

% compute the Matrix Formalism signal 
% for ngdir directions and interpolate to 900 directions uniformly distributed on the sphere.
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
% 
% Output:
%     1. points (ngdir directions)
%     2. SIG_EIG_cmpts_alldir
%     3. SIG_EIG_allcmpts_alldir

Ncmpt = mymesh.Ncmpt;
[points,C,v] = spheresurface_regularpoints(1,experi_btpde.ngdir_total);
ngdir_total = size(points,1);
ii = find(points(:,3) >= 0);
% negii
for j = 1:size(ii,1)
    for k = 1:ngdir_total
        if (norm(points(j,1:2)-points(k,1:2)) < 1e-10  && points(j,3)+points(k,3) < 1e-10)
            negii(ii(j)) = k;
        end
    end
end
jc = 0;
for ic = 1:size(C,1)
    jj = find(C(ic,1) == ii);
    kk = find(C(ic,2) == ii);
    ll = find(C(ic,3) == ii);
    if (~isempty(jj) & ~isempty(kk) & ~isempty(ll))
        Ckeep(jc+1,1:3) = C(ic,1:3);
        jc = jc+1;
    end
end
graddir_index = ii;
ndir = length(graddir_index);
nexperi = length(experi_btpde.sdeltavec);
nb = size(experi_btpde.bvalues,2);

SIG_EIG_cmpts_alldir = nan*ones(ngdir_total,Ncmpt,nexperi,nb);
SIG_EIG_allcmpts_alldir = nan*ones(ngdir_total,nexperi,nb);
for idir = 1:ndir
    experi_btpde.gdir = points(graddir_index(idir),:)';    
    [SIG_EIG_cmpts,SIG_EIG_allcmpts] = BTPDE_EIG(experi_btpde,mymesh,...
    IC_cmpts,EIG_value_cmpts,EIG_proj_cmpts);   
    SIG_EIG_cmpts_alldir(graddir_index(idir),:,:,:) = SIG_EIG_cmpts;
    SIG_EIG_cmpts_alldir(negii(graddir_index(idir)),:,:,:) = SIG_EIG_cmpts;
    SIG_EIG_allcmpts_alldir(graddir_index(idir),:,:) = SIG_EIG_allcmpts(:,:);
    SIG_EIG_allcmpts_alldir(negii(graddir_index(idir)),:,:) = SIG_EIG_allcmpts(:,:);
end
SIG_EIG_cmpts_alldir(find(SIG_EIG_cmpts_alldir==0)) = nan;
SIG_EIG_allcmpts_alldir(find(SIG_EIG_allcmpts_alldir==0)) = nan;