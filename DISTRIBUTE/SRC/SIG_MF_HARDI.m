function [SIG_EIG_cmpts_alldir,SIG_EIG_allcmpts_alldir,ctime_alldir] ...
    = SIG_MF_HARDI(experiment,VOL_cmpts,...
    IC_cmpts,EIG_value_cmpts,EIG_proj_cmpts,...
    points_gdir,graddir_index,negii)

% compute the Matrix Formalism signal 
% for ngdir_total directions uniformly distributed on the sphere.
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
%     5. EIG_proj_cmpts
%     6. points_gdir
%     7. graddir_index
%     8. negii
%
% Output:
%     1. SIG_EIG_cmpts_alldir
%     2. SIG_EIG_allcmpts_alldir
%     3. ctime_alldir

global BDELTA SDELTA SEQ OGSEPER

sdeltavec = experiment.sdeltavec;
bdeltavec = experiment.bdeltavec;
seqvec = experiment.seqvec;
npervec = experiment.npervec;

Ncmpt = length(VOL_cmpts);

VOL_allcmpts = sum(VOL_cmpts);
VOL_frac = VOL_cmpts/VOL_allcmpts;

nexperi = length(sdeltavec);
nb = size(experiment.bvalues,2);

ngdir_total = size(points_gdir,1);
ndir = length(graddir_index);

SIG_EIG_cmpts_alldir = nan*ones(ngdir_total,Ncmpt,nexperi,nb);


ctime_alldir = nan*ones(Ncmpt,ngdir_total,nexperi,nb);

for icmpt = 1:Ncmpt
    L_mat = diag(EIG_value_cmpts{icmpt});
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
        
        for ib = 1:nb
            
            gv = sqrt(bvec(ib)/bvaluenoq(iexperi));
            tic
            for idir = 1:ndir
                
                b_start_time = clock;
                
                gdir = points_gdir(graddir_index(idir),:)';
                gdir = gdir/norm(gdir);
                W_mat = squeeze(EIG_proj_cmpts{icmpt}(1,:,:)*gdir(1)...
                    +EIG_proj_cmpts{icmpt}(2,:,:)*gdir(2)...
                    +EIG_proj_cmpts{icmpt}(3,:,:)*gdir(3));                                
                V_mat = L_mat+sqrt(-1)*W_mat*gv;                
                [V,D] = eig(V_mat);                
                % V_mat_trans = L_mat-sqrt(-1)*W_mat*gv;                
                invV = inv(V);
                H1 = V(1,:)*diag(exp(-SDELTA*diag(D)))...
                    *invV*diag(exp(-(BDELTA-SDELTA)*EIG_value_cmpts{icmpt}))...
                    *(invV'*diag(exp(-SDELTA*diag(D)))'*V(1,:)');              
                SIG_EIG_cmpts_alldir(graddir_index(idir),icmpt,iexperi,ib)...
                    = H1(1,1)*VOL_cmpts(icmpt)*IC_cmpts(icmpt);                
                if (~isempty(negii{idir}))                
                    SIG_EIG_cmpts_alldir(negii{idir},icmpt,iexperi,ib) ...
                    = SIG_EIG_cmpts_alldir(graddir_index(idir),icmpt,iexperi,ib);
                end
                 
                ctime_alldir(icmpt,graddir_index(idir),iexperi,ib) = etime(clock, b_start_time);
                
            end
            toc
        end
    end
end

SIG_EIG_allcmpts_alldir = zeros(ngdir_total,nexperi,nb);
for iexperi = 1:nexperi
    for ib = 1:nb
        for icmpt = 1:Ncmpt
            SIG_EIG_allcmpts_alldir(:,iexperi,ib) ...
                = SIG_EIG_allcmpts_alldir(:,iexperi,ib) ...
                + squeeze(SIG_EIG_cmpts_alldir(:,icmpt,iexperi,ib));
        end
    end
end
