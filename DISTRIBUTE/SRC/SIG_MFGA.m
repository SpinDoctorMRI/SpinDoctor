function [SIG_MFGA_cmpts,SIG_MFGA_allcmpts,ctime] = SIG_MFGA(experiment,VOL_cmpts,...
    IC_cmpts,DTENSOR_cmpts)

% compute the Matrix Formalism GAUSSIAN APPROXIMATION signal 
%
% Input:
%     1. experiment is a structure with 10 elements:
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
%
% Output:
%     1. SIG_MFGA_cmpts
%     2. SIG_MFGA_allcmpts
%     3. ctime

gdir = experiment.gdir;
sdeltavec = experiment.sdeltavec;
bdeltavec = experiment.bdeltavec;
seqvec = experiment.seqvec;
npervec = experiment.npervec;


Ncmpt = length(VOL_cmpts);

VOL_allcmpts = sum(VOL_cmpts);
VOL_frac = VOL_cmpts/VOL_allcmpts;

gdir = gdir/norm(gdir);

nexperi = length(sdeltavec);

nb = size(experiment.bvalues,2);

SIG_MFGA_cmpts = zeros(Ncmpt,nexperi,nb);
MF_EIG_allcmpts = zeros(nexperi,nb);

ctime = nan*ones(Ncmpt,nexperi,nb);
% MF_GA: Matrix Formalism Gaussian Approximation
nexperi = length(experiment.sdeltavec);
for iexperi = 1:nexperi
    bvec = experiment.bvalues(iexperi,:);
    nb = length(bvec);
    for ib = 1:nb
        for icmpt = 1:Ncmpt
            
            b_start_time = clock;
            ADC_TMP_cmpts(icmpt,iexperi) = experiment.gdir'*squeeze(DTENSOR_cmpts(icmpt,iexperi,:,:))*experiment.gdir;
            SIG_MFGA_cmpts(icmpt,iexperi,ib) = VOL_cmpts(icmpt)*IC_cmpts(icmpt)*exp(-ADC_TMP_cmpts(icmpt,iexperi)*bvec(ib));
            ctime(icmpt,iexperi,ib) = etime(clock, b_start_time);
            
        end
        SIG_MFGA_allcmpts(iexperi,ib) = 0;
        for icmpt = 1:Ncmpt
            SIG_MFGA_allcmpts(iexperi,ib) = SIG_MFGA_allcmpts(iexperi,ib) + SIG_MFGA_cmpts(icmpt,iexperi,ib);
        end
    end
end
