function [SIG_MF_cmpts,SIG_MF_allcmpts,ctime] = SIG_MF(experiment,VOL_cmpts,...
    IC_cmpts,EIG_value_cmpts,EIG_proj_cmpts)

% compute the Matrix Formalism signal 
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
%     4. EIG_value_cmpts
%
% Output:
%     1. SIG_MF_cmpts
%     2. SIG_MF_allcmpts
%     3. ctime

global QVAL UG 
global BDELTA SDELTA SEQ OGSEPER 

gdir = experiment.gdir;
sdeltavec = experiment.sdeltavec;
bdeltavec = experiment.bdeltavec;
seqvec = experiment.seqvec;
npervec = experiment.npervec;
UG = gdir';
UG = UG/norm(UG);

Ncmpt = length(VOL_cmpts);

VOL_allcmpts = sum(VOL_cmpts);
VOL_frac = VOL_cmpts/VOL_allcmpts;


gdir = gdir/norm(gdir);
nexperi = length(sdeltavec);

nb = size(experiment.bvalues,2);

SIG_MF_cmpts = zeros(Ncmpt,nexperi,nb);
ctime = nan*ones(Ncmpt,nexperi,nb);

for icmpt = 1:Ncmpt
    nexperi = length(sdeltavec);
    for iexperi = 1:nexperi
               
        bvec = experiment.bvalues(iexperi,:);
        
        nb = length(bvec);
        SDELTA = sdeltavec(iexperi);
        BDELTA = bdeltavec(iexperi);
        TE = SDELTA+BDELTA;
        SEQ = seqvec(iexperi);% for choosing case PGSE, OGSEcos or OGSEsin
        omega = 2*pi*npervec(iexperi)/SDELTA;
        OGSEPER = 1./omega*2*pi;%% set up number for OGSE
        bvaluenoq(iexperi) = seqbvaluenoq;
        
        W_mat = squeeze(EIG_proj_cmpts{icmpt}(1,:,:)*gdir(1)...
            +EIG_proj_cmpts{icmpt}(2,:,:)*gdir(2)...
            +EIG_proj_cmpts{icmpt}(3,:,:)*gdir(3));
        L_mat = diag(EIG_value_cmpts{icmpt});
        
        for ib = 1:nb
            tic           
            b_start_time = clock;
             
            gv = sqrt(bvec(ib)/bvaluenoq(iexperi));            
            V_mat = L_mat+sqrt(-1)*W_mat*gv;            
            [V,D] = eig(V_mat);        
%             if (icmpt == Ncmpt)
%                 diag(D)
%             end
            %V_mat_trans = L_mat-sqrt(-1)*W_mat*gv;
            invV = inv(V);            
            H1 = V(1,:)*diag(exp(-SDELTA*diag(D)))*invV*diag(exp(-(BDELTA-SDELTA)*EIG_value_cmpts{icmpt}))*(invV'*diag(exp(-SDELTA*diag(D)))'*V(1,:)');            
            SIG_MF_cmpts(icmpt,iexperi,ib)= H1(1,1)*VOL_cmpts(icmpt)*IC_cmpts(icmpt);
                        
            ctime(icmpt,iexperi,ib) = etime(clock, b_start_time);
            toc
        end
        
    end
    
end

SIG_MF_allcmpts = zeros(nexperi,nb);
for iexperi = 1:nexperi
    for ib = 1:nb
        for icmpt = 1:Ncmpt
            SIG_MF_allcmpts(iexperi,ib) = SIG_MF_allcmpts(iexperi,ib) + SIG_MF_cmpts(icmpt,iexperi,ib);
        end
    end
end