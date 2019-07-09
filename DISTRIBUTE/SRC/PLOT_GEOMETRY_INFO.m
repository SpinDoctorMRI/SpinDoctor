function PLOT_GEOMETRY_INFO(cmpts_bdys_mat,OUT_cmpts_index,IN_cmpts_index,ECS_cmpts_index,...
	VOL_cmpts,SA_cmpts,SAu_cmpts)

% plot information of the geometry
% 
% Input:
%     1. cmpts_bdys_mat
%     2. OUT_cmpts_index
%     3. IN_cmpts_index
%     4. ECS_cmpts_index
%     5. VOL_cmpts
%     6. SA_cmpts
%     7. SAu_cmpts
%     
% Output:
%     1. 1 figure with title of "Connections boundary-compartment"
%     2. 1 figure with title of "VOL_cmpts"
%     3. 1 figure with title of "Surface Area"
%     4. 1 figure with title of "SAu_cmpts/VOL_cmpts"	

Ncmpt = length(VOL_cmpts);
Nboundary = size(cmpts_bdys_mat,2);

figure; 
hold all;
cmpts_bdys_mat_plot=zeros(size(cmpts_bdys_mat));
cmpts_bdys_mat_plot(OUT_cmpts_index,:)=cmpts_bdys_mat(OUT_cmpts_index,:);
spy(cmpts_bdys_mat_plot,'ko',10);


cmpts_bdys_mat_plot(:,:)=0;
cmpts_bdys_mat_plot(ECS_cmpts_index,:)=cmpts_bdys_mat(ECS_cmpts_index,:);
spy(cmpts_bdys_mat_plot,'ks',10);

cmpts_bdys_mat_plot(:,:)=0;
cmpts_bdys_mat_plot(IN_cmpts_index,:)=cmpts_bdys_mat(IN_cmpts_index,:);
spy(cmpts_bdys_mat_plot,'kd',10); 


xlabel('iboundary'); ylabel('icmpt');
title('Connections boundary-compartment');
set(gca,'Ytick',[1:Ncmpt]);
set(gca,'Xtick',[1:Nboundary]);
grid on;

figure; 
hold on;
bar(1:Ncmpt,VOL_cmpts,'b');
bar(Ncmpt+1,sum(VOL_cmpts),'r');
title('VOL\_cmpts');
set(gca,'Xtick',[1:Ncmpt+1]);
xlabel('icmpt (last: all cmpts) ');
grid on;

figure; 
hold on;
bar(1:Ncmpt,SA_cmpts,'b');
bar(Ncmpt+1,sum(SA_cmpts),'r');
title('Surface Area');
set(gca,'Xtick',[1:Ncmpt+1]);
xlabel('icmpt (last: all cmpts)');
grid on;

figure; 
hold on;
bar(1:Ncmpt,SAu_cmpts./VOL_cmpts','b');
bar(Ncmpt+1,sum(SAu_cmpts)./sum(VOL_cmpts),'r');
title('SAu\_cmpts/VOL\_cmpts');
set(gca,'Xtick',[1:Ncmpt+1]);
xlabel('icmpt (last: all cmpts)');
grid on;