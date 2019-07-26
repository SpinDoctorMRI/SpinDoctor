function [SIG_EIG_cmpts_alldir,SIG_EIG_allcmpts_alldir,ctime_alldir] ...
    = SIG_MFGA_HARDI(experi_btpde,VOL_cmpts,IC_cmpts,DTENSOR_cmpts,...
    points_gdir,graddir_index,negii)

% compute the Matrix Formalism GAUSSIAN APPROXIMATION signal 
% for ngdir_total directions uniformly distributed on the sphere.
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
%     2. VOL_cmpts
%     3. IC_cmpts
%     4. DTENSOR_cmpts
%     5. points_gdir
%     6. graddir_index
%     7. negii
%
% Output:
%     1. SIG_EIG_cmpts_alldir
%     2. SIG_EIG_allcmpts_alldir
%     3. ctime_alldir

Ncmpt = length(VOL_cmpts);

nexperi = length(experi_btpde.sdeltavec);
nb = size(experi_btpde.bvalues,2);

ngdir_total = size(points_gdir,1);
ndir = length(graddir_index);

SIG_EIG_cmpts_alldir = nan*ones(ngdir_total,Ncmpt,nexperi,nb);
SIG_EIG_allcmpts_alldir = nan*ones(ngdir_total,nexperi,nb);

ctime_alldir = nan*ones(Ncmpt,ngdir_total,nexperi,nb);
tic
for idir = 1:ndir
    experi_btpde.gdir = points_gdir(graddir_index(idir),:)';    
    experi_btpde.gdir = experi_btpde.gdir/norm(experi_btpde.gdir);
    
    [SIG_EIG_cmpts,SIG_EIG_allcmpts,ctime] = SIG_MFGA(experi_btpde,VOL_cmpts,...
        IC_cmpts,DTENSOR_cmpts);
    ctime_alldir(:,graddir_index(idir),:,:) = ctime;
    
    SIG_EIG_cmpts_alldir(graddir_index(idir),:,:,:) = SIG_EIG_cmpts;
    SIG_EIG_allcmpts_alldir(graddir_index(idir),:,:) = SIG_EIG_allcmpts(:,:);
    if (~isempty(negii{idir}))
        SIG_EIG_cmpts_alldir(negii{idir},:,:,:) = SIG_EIG_cmpts;
        SIG_EIG_allcmpts_alldir(negii{idir},:,:) = SIG_EIG_allcmpts(:,:);
    end
end
toc
SIG_EIG_cmpts_alldir(find(SIG_EIG_cmpts_alldir==0)) = nan;
SIG_EIG_allcmpts_alldir(find(SIG_EIG_allcmpts_alldir==0)) = nan;