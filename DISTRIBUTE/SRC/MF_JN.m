function [MF_JN_cmpts] = MF_JN(EIG_value_cmpts,DIFF_cmpts,experiment)

% compute the quantity J(lambda_n,f)
% 
% Input:
%     1. EIG_value_cmpts
%     2. DIFF_cmpts
%     3. experiment is a structure with 10 elements:
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
% 
% Output:
%     MF_JN_cmpts


Ncmpt = length(DIFF_cmpts);

global BDELTA SDELTA SEQ OGSEPER

gdir = experiment.gdir;
sdeltavec = experiment.sdeltavec;
bdeltavec = experiment.bdeltavec;
seqvec = experiment.seqvec;
npervec = experiment.npervec;

for icmpt = 1:Ncmpt
    neig = length(EIG_value_cmpts{icmpt});
    nexperi = length(sdeltavec);
    for iexperi = 1:nexperi
        SDELTA = sdeltavec(iexperi);
        BDELTA = bdeltavec(iexperi);
        TE = SDELTA+BDELTA;
        SEQ = seqvec(iexperi);% for choosing case PGSE, OGSEcos or OGSEsin
        omega = 2*pi*npervec(iexperi)/SDELTA;
        OGSEPER = 1./omega*2*pi;%% set up number for OGSE
        bvaluenoq(iexperi) = seqbvaluenoq;
        dtime = 10;
        TLIST = linspace(0,SDELTA+BDELTA,round((SDELTA+BDELTA)/dtime)+1);
        for ieig = 1:neig
            const = EIG_value_cmpts{icmpt}(ieig);
            if (seqvec(iexperi) == 1)
                if (abs(const)<1e-16)
                    dtmp = 0;
                else
                    ff = -1/(const^2)*(exp(-const*(BDELTA+SDELTA))+...
                        exp(-const*(BDELTA-SDELTA))-...
                        2*(const*SDELTA+exp(-const*SDELTA)+exp(-const*BDELTA)-1));
                    dtmp = 1/bvaluenoq(iexperi)*ff;
                end
            else
                fint = zeros(1,length(TLIST));
                if (const*TLIST(end)<=1e-3)
                    for idt=1:length(TLIST)
                        fint(idt) = fexp_taylor(const,TLIST(idt));
                    end
                else
                    for idt=1:length(TLIST)
                        fint(idt) = fexp(const,TLIST(idt),TE);
                    end
                end
                functime_eig{iexperi}(ieig,:) = 1/bvaluenoq(iexperi)*(fint);
                fint = functime_eig{iexperi}(ieig,:);
                dtmp = dtime*(sum(fint(2:end-1))+fint(1)/2+fint(end)/2);
            end
            MF_JN_cmpts{icmpt}(iexperi,ieig) = dtmp/DIFF_cmpts(icmpt);
        end
    end
    
end

