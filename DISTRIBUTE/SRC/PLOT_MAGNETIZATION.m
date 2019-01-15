function PLOT_MAGNETIZATION(mymesh,YOUT,OUT_cmpts_index,ECS_cmpts_index,IN_cmpts_index)

xmin=0;
xmax=0;
ymin=0;
ymax=0;
zmin=0;
zmax=0;
for icmpt = 1:mymesh.Ncmpt
    xx=max(mymesh.Pts_cmpt_reorder{icmpt}(1,:));
    yx=max(mymesh.Pts_cmpt_reorder{icmpt}(2,:));
    zx=max(mymesh.Pts_cmpt_reorder{icmpt}(3,:));
    xn=min(mymesh.Pts_cmpt_reorder{icmpt}(1,:));
    yn=min(mymesh.Pts_cmpt_reorder{icmpt}(2,:));
    zn=min(mymesh.Pts_cmpt_reorder{icmpt}(3,:));
    xmin = min(xmin,xn);
    xmax = max(xmax,xx);
    ymin = min(ymin,yn);
    ymax = max(ymax,yx);
    zmin = min(zmin,zn);
    zmax = max(zmax,zx);
end
if (~isempty(OUT_cmpts_index))
figure; 
hold on;
cmptvec = OUT_cmpts_index;
for ict = 1:length(cmptvec)
    icmpt = cmptvec(ict);
    Fac = [];
    for iboundary = 1:mymesh.Nboundary
        Fac = [Fac,mymesh.Fac_boundary_reorder{icmpt}{iboundary}];
    end
    h = trisurf(Fac',mymesh.Pts_cmpt_reorder{icmpt}(1,:),mymesh.Pts_cmpt_reorder{icmpt}(2,:),...
        mymesh.Pts_cmpt_reorder{icmpt}(3,:),real(YOUT{end}{end}{icmpt}(:,end)));
    set(h,'facealpha',0.6);
    axis equal;
    axis([xmin,xmax,ymin,ymax,zmin,zmax]); colorbar('southoutside');
    view(3);
    title(['Inner cmpts: ',num2str(OUT_cmpts_index)]);
end
end
if (~isempty(IN_cmpts_index))
figure; 
hold on;
cmptvec = IN_cmpts_index;
for ict = 1:length(cmptvec)
    icmpt = cmptvec(ict);
    Fac = [];
    for iboundary = 1:mymesh.Nboundary
        Fac = [Fac,mymesh.Fac_boundary_reorder{icmpt}{iboundary}];
    end
    h = trisurf(Fac',mymesh.Pts_cmpt_reorder{icmpt}(1,:),mymesh.Pts_cmpt_reorder{icmpt}(2,:),...
        mymesh.Pts_cmpt_reorder{icmpt}(3,:),real(YOUT{end}{end}{icmpt}(:,end)));
    set(h,'facealpha',0.6);
    axis equal;
    axis([xmin,xmax,ymin,ymax,zmin,zmax]); colorbar('southoutside');
    view(3);
    title(['Outer cmpts: ',num2str([IN_cmpts_index])]);
end
end
if (~isempty(ECS_cmpts_index))
figure; 
hold on;
cmptvec = ECS_cmpts_index;
for ict = 1:length(cmptvec)
    icmpt = cmptvec(ict);
    Fac = [];
    for iboundary = 1:mymesh.Nboundary
        Fac = [Fac,mymesh.Fac_boundary_reorder{icmpt}{iboundary}];
    end
    h = trisurf(Fac',mymesh.Pts_cmpt_reorder{icmpt}(1,:),mymesh.Pts_cmpt_reorder{icmpt}(2,:),...
        mymesh.Pts_cmpt_reorder{icmpt}(3,:),real(YOUT{end}{end}{icmpt}(:,end)));
    set(h,'facealpha',0.6);
    axis equal;
    axis([xmin,xmax,ymin,ymax,zmin,zmax]); colorbar('southoutside');
    view(3);
    title(['Box cmpt: ',num2str([ECS_cmpts_index])]);
end
end
