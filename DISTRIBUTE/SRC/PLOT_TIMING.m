function PLOT_TIMING(elapsed_time,mymesh)

% plot elapsed time
% 
% Input:
%     1. elapsed_time
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
% 
% Output:
%     1 figure of computational times of nexperi experiments 

nexperi = size(elapsed_time,2);
nb = size(elapsed_time,1);
figure;
for iexperi = 1:nexperi
    subplot(nexperi,1,iexperi); hold on;
    bar(1:nb,[elapsed_time(:,iexperi)],'b');
	bar(nb+1,[sum(elapsed_time(:,iexperi))],'r');
    title(['Experi ',num2str(iexperi),': ',...
        num2str(mymesh.Nnode(1)),' nodes, ',...
        num2str(mymesh.Nele(1)), ' elements']);
    ylabel('Computational time (s)');
    xlabel('Index (last: total)');
end