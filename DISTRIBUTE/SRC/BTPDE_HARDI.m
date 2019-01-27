function [points,ADC_BT_cmpts_alldir,ADC_BT_allcmpts_alldir] ...
    = BTPDE_HARDI(experi_btpde,mymesh,DIFF_cmpts,kappa_bdys,IC_cmpts)

% compute the ADC from Bloch-Torrey equation for ngdir directions and interpolate to 900 directions uniformly distributed on the sphere.
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
%     2. ADC_BT_cmpts_alldir
%     3. ADC_BT_allcmpts_alldir

Ncmpt = length(DIFF_cmpts);
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

ADC_BT_cmpts_alldir = nan*ones(ngdir_total,Ncmpt,nexperi);
ADC_BT_allcmpts_alldir = nan*ones(ngdir_total,nexperi);
for idir = 1:ndir
    experi_btpde.gdir = points(graddir_index(idir),:)';
    [TOUT,YOUT,MF_cmpts,MF_allcmpts,difftime,BTPDE_elapsed_time] ...
        = BTPDE(experi_btpde,mymesh,DIFF_cmpts,kappa_bdys,IC_cmpts);
    [ADC_BT_cmpts,ADC_BT_allcmpts,ADC_BT_allcmpts_S0] = FIT_SIGNAL(MF_cmpts,MF_allcmpts,experi_btpde.bvalues);
    ADC_BT_cmpts_alldir(graddir_index(idir),:,:) = ADC_BT_cmpts;
    ADC_BT_cmpts_alldir(negii(graddir_index(idir)),:,:) = ADC_BT_cmpts;
    ADC_BT_allcmpts_alldir(graddir_index(idir),:) = ADC_BT_allcmpts(:,1)';
    ADC_BT_allcmpts_alldir(negii(graddir_index(idir)),:) = ADC_BT_allcmpts(:,1)';
end
ADC_BT_cmpts_alldir(find(ADC_BT_cmpts_alldir==0)) = nan;
ADC_BT_allcmpts_alldir(find(ADC_BT_allcmpts_alldir==0)) = nan;