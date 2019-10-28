function [ADC_HADC_cmpts_alldir,ADC_HADC_allcmpts_alldir,ctime_alldir] ...
    = HADC_HARDI(experi_hadc,mymesh,DIFF_cmpts,IC_cmpts,...
    points_gdir,graddir_index,negii)

% compute the ADC from HADC model for ngdir directions.
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
%     5. points_gdir
%     6. graddir_index
%     7. negii
%
% Output:
%     1. ADC_HADC_cmpts_alldir
%     2. ADC_HADC_allcmpts_alldir
%     3. ctime_alldir

Ncmpt = length(DIFF_cmpts);
nexperi = length(experi_hadc.sdeltavec);

ngdir_total = size(points_gdir,1);
ndir = length(graddir_index);


ADC_HADC_cmpts_alldir = nan*ones(ngdir_total,Ncmpt,nexperi);
ADC_HADC_allcmpts_alldir = nan*ones(ngdir_total,nexperi);
ctime_alldir = nan*ones(ngdir_total,nexperi);

for idir = 1:ndir
    experi_hadc.gdir = points_gdir(graddir_index(idir),:)';
    experi_hadc.gdir = experi_hadc.gdir/norm(experi_hadc.gdir);      
    [ADC_HADC_cmpts,ADC_HADC_allcmpts,ctime] ...
        = HADC(experi_hadc,mymesh,DIFF_cmpts,IC_cmpts);
    ctime_alldir(graddir_index(idir),:) = ctime;    
    ADC_HADC_cmpts_alldir(graddir_index(idir),:,:) = ADC_HADC_cmpts;
    
    ADC_HADC_allcmpts_alldir(graddir_index(idir),:) = ADC_HADC_allcmpts(:,1)';
     
    if (~isempty(negii{idir}))
        ADC_HADC_cmpts_alldir(negii(graddir_index(idir)),:,:) = ADC_HADC_cmpts;
        ADC_HADC_allcmpts_alldir(negii(graddir_index(idir)),:) = ADC_HADC_allcmpts(:,1)';
    end   
    
end
