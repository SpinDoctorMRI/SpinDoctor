function [points,ADC_HADC_cmpts_alldir,ADC_HADC_allcmpts_alldir] ...
    = HADC_HARDI(experi_hadc,mymesh,DIFF_cmpts,IC_cmpts)

% compute the ADC from HADC model for ngdir directions and interpolate to 900 directions uniformly distributed on the sphere.
% 
% Input:
%     1. experiment_hadc is a structure with 8 elements:
%         ngdir_total 
%         gdir       
%         sdeltavec   
%         bdeltavec    
%         seqvec       
%         npervec     
%         rtol        
%         atol               
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
%     4. IC_cmpts
% 
% Output:
%     1. points (ngdir directions)
%     2. ADC_HADC_cmpts_alldir
%     3. ADC_HADC_allcmpts_alldir

Ncmpt = length(DIFF_cmpts);
[points,C,v] = spheresurface_regularpoints(1,experi_hadc.ngdir_total);
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
nexperi = length(experi_hadc.sdeltavec);

ADC_HADC_cmpts_alldir = nan*ones(ngdir_total,Ncmpt,nexperi);
ADC_HADC_allcmpts_alldir = nan*ones(ngdir_total,nexperi);
for idir = 1:ndir
    experi_hadc.gdir = points(graddir_index(idir),:)';    
    [ADC_HADC_cmpts,ADC_HADC_allcmpts,HADC_elapsed_time] ...
        = HADC(experi_hadc,mymesh,DIFF_cmpts,IC_cmpts);
    ADC_HADC_cmpts_alldir(graddir_index(idir),:,:) = ADC_HADC_cmpts;
    ADC_HADC_cmpts_alldir(negii(graddir_index(idir)),:,:) = ADC_HADC_cmpts;
    ADC_HADC_allcmpts_alldir(graddir_index(idir),:) = ADC_HADC_allcmpts(:,1)';
    ADC_HADC_allcmpts_alldir(negii(graddir_index(idir)),:) = ADC_HADC_allcmpts(:,1)';
        
end
ADC_HADC_cmpts_alldir(find(ADC_HADC_cmpts_alldir==0)) = nan;
ADC_HADC_allcmpts_alldir(find(ADC_HADC_allcmpts_alldir==0)) = nan;