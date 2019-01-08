function PLOT_GEOMETRY_INFO(boundary_mat,Cell_cmpt,Nucleus_cmpt,Box_cmpt,...
	VOL,SA,SAu)
	
Ncmpt = length(VOL);
Nboundary = size(boundary_mat,2);

figure; 
hold all;
boundary_mat_plot=zeros(size(boundary_mat));
boundary_mat_plot(Cell_cmpt,:)=boundary_mat(Cell_cmpt,:);
spy(boundary_mat_plot,'ko',10);


boundary_mat_plot(:,:)=0;
boundary_mat_plot(Box_cmpt,:)=boundary_mat(Box_cmpt,:);
spy(boundary_mat_plot,'ks',10);

boundary_mat_plot(:,:)=0;
boundary_mat_plot(Nucleus_cmpt,:)=boundary_mat(Nucleus_cmpt,:);
spy(boundary_mat_plot,'kd',10); 


xlabel('iboundary'); ylabel('icmpt');
title('Connections boundary-compartment');
set(gca,'Ytick',[1:Ncmpt]);
set(gca,'Xtick',[1:Nboundary]);
grid on;

figure; 
hold on;
bar(1:Ncmpt,VOL,'b');
bar(Ncmpt+1,sum(VOL),'r');
title('VOL');
set(gca,'Xtick',[1:Ncmpt+1]);
xlabel('icmpt (last: all cmpts) ');
grid on;

figure; 
hold on;
bar(1:Ncmpt,SA,'b');
bar(Ncmpt+1,sum(SA),'r');
title('Surface Area');
set(gca,'Xtick',[1:Ncmpt+1]);
xlabel('icmpt (last: all cmpts)');
grid on;

figure; 
hold on;
bar(1:Ncmpt,SAu./VOL','b');
bar(Ncmpt+1,sum(SAu)./sum(VOL),'r');
title('SAu/VOL');
set(gca,'Xtick',[1:Ncmpt+1]);
xlabel('icmpt (last: all cmpts)');
grid on;