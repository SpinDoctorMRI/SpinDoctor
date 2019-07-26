function PLOT_DTENSOR(DTENSOR_cmpts,DTENSOR_allcmpts,DIFF_cmpts)

% plot the diffusion tensor
%
% Input:
%     1. DTENSOR_cmpts
%     2. DTENSOR_allcmpts
%     3. DIFF_cmpts
%
% Output:
%     1 figure for TENSOR of nexperi experiments

nexperi = size(DTENSOR_allcmpts,1);
Ncmpt = length(DIFF_cmpts);
figure; hold on;
for iexperi = 1:nexperi   
    for icmpt = 1:Ncmpt
        subplot(Ncmpt,nexperi,icmpt+(iexperi-1)*Ncmpt); hold on;
        if icmpt == 1
            title(['Experi ',num2str(iexperi)]);
        end
        yvec = squeeze(DTENSOR_cmpts(icmpt,iexperi,:,:))/DIFF_cmpts(icmpt);
        yvec = yvec(:);
        yc = zeros(6,1);
        yo = zeros(6,1);
        yc([1:3]) = yvec([1,5,9]);
        yo([4:6]) = yvec([2,3,6]);
        bar(1:6,yc,'r');
        bar(1:6,yo,'b');
        %title(['Cmpt ',num2str(icmpt)]);
        set(gca,'ylim',[min(0,min(yvec)),1]);
        set(gca,'Ytick',linspace(0,1,3));
        set(gca,'Xtick',[1:6]);
        ylabel('');
        xlabel('matrix index');
        grid on;
    end
end



