function [ADC_cmpts,ADC_allcmpts,ADC_allcmpts_S0] = FIT_SIGNAL(MF_cmpts,MF_allcmpts,bvalues)

% fit the ADC from the dMRI signal
% 
% Input:
%     1. MF_cmpts
%     2. MF_allcmpts
%     3. bvalues
%     
% Output:
%     1. ADC_cmpts
%     2. ADC_allcmpts
%     3. ADC_allcmpts_S0

nexperi = size(MF_cmpts,2);
nb = size(MF_cmpts,3);
Ncmpt = size(MF_cmpts,1);
% for iexperi = 1:nexperi
    % bvec = bvalues(iexperi,:);  
    % nb = length(bvec);
    % for ib = 1:nb
        % for icmpt = 1:Ncmpt
            % MF_cmpts(icmpt,iexperi,ib) = MT{iexperi}{ib}{icmpt}(end);
            % M0(icmpt,iexperi,ib) = MT{iexperi}{ib}{icmpt}(1);
        % end
        % MF_allcmpts(iexperi,ib) = 0;
        % for icmpt = 1:Ncmpt
            % MF_allcmpts(iexperi,ib) = MF_allcmpts(iexperi,ib) + MF_cmpts(icmpt,iexperi,ib);
        % end
        % M0_allcmpts(iexperi,ib) = 0;
        % for icmpt = 1:Ncmpt
            % M0_allcmpts(iexperi,ib) = M0_allcmpts(iexperi,ib) + M0(icmpt,iexperi,ib);
        % end
    % end   
    % ib0 = find(abs(bvec)<=1e-16);
    % ibn0 = find(abs(bvec)>1e-16);
    % if (length(ib0) >= 1)
        % for icmpt = 1:Ncmpt
            % S0(icmpt,iexperi) = MF_cmpts(icmpt,iexperi,ib0(1));
        % end
        % S0_allcmpts(iexperi) = MF_allcmpts(iexperi,ib0(1));
    % else
        % S0(1:Ncmpt,iexperi) = nan;
        % S0_allcmpts(iexperi) = nan;
    % end
% end
ADC_cmpts = nan*ones(Ncmpt,nexperi);
ADC_allcmpts = nan*ones(nexperi,1);
ADC_polydeg = nan*ones(Ncmpt,nexperi);
ADC_allcmpts_polydeg = nan*ones(nexperi,1);
for iexperi = 1:nexperi
    bvec = bvalues(iexperi,:);
    if (length(bvec) >= 2)       
        bmin = bvec(1);
        bmax = bvec(end);        
        for icmpt = 1:Ncmpt
            data1d = real(squeeze(MF_cmpts(icmpt,iexperi,:)))';
            [fit_poly,ADC01d,KUR1d,KUR01d,S01d,Cfit1d,errfit,ndeg,ADC0_err1d,KUR_err1d] ...
                = process_signal_POLY(data1d,bvec,bmin,bmax);
            ADC_cmpts(icmpt,iexperi) = ADC01d;
            ADC_polydeg(icmpt,iexperi) = ndeg;
            ADC_S0(icmpt,iexperi) = S01d;
        end
        data1d = real(MF_allcmpts(iexperi,:));
        [fit_poly,ADC01d,KUR1d,KUR01d,S01d,Cfit1d,errfit,ndeg,ADC0_err1d,KUR_err1d] ...
            = process_signal_POLY(data1d,bvec,bmin,bmax);
        ADC_allcmpts(iexperi,1) = ADC01d;
        ADC_allcmpts_polydeg(iexperi,1) = ndeg;
        ADC_allcmpts_S0(iexperi,1) = S01d;
    end
end