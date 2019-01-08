function PLOT_FEMESH(mymesh,Cell_cmpt,Box_cmpt,Nucleus_cmpt)

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

if (~isempty(Cell_cmpt))
figure; 
%subplot(1,3,1); 
hold on;
cmptvec = Cell_cmpt;
for ict = 1:length(cmptvec)
    icmpt = cmptvec(ict);
    Fac = [];
    for iboundary = 1:mymesh.Nboundary
        Fac = [Fac,mymesh.Fac_boundary_reorder{icmpt}{iboundary}];
    end
    h = trisurf(Fac',mymesh.Pts_cmpt_reorder{icmpt}(1,:),mymesh.Pts_cmpt_reorder{icmpt}(2,:),...
        mymesh.Pts_cmpt_reorder{icmpt}(3,:));
    set(h,'facealpha',0.9);
    axis equal;
    axis([xmin,xmax,ymin,ymax,zmin,zmax]); 
    view(3);
    title(['FE Mesh Inner cmpts: ',num2str(Cell_cmpt)]);
end
end
if (~isempty(Nucleus_cmpt))
figure;
%subplot(1,3,2); 
hold on;
cmptvec = Nucleus_cmpt;
for ict = 1:length(cmptvec)
    icmpt = cmptvec(ict);
    Fac = [];
    for iboundary = 1:mymesh.Nboundary
        Fac = [Fac,mymesh.Fac_boundary_reorder{icmpt}{iboundary}];
    end
    h = trisurf(Fac',mymesh.Pts_cmpt_reorder{icmpt}(1,:),mymesh.Pts_cmpt_reorder{icmpt}(2,:),...
        mymesh.Pts_cmpt_reorder{icmpt}(3,:));
    set(h,'facealpha',0.9);
    axis equal;
    axis([xmin,xmax,ymin,ymax,zmin,zmax]);
    view(3);
    title(['FE Mesh Outer cmpts: ',num2str([Nucleus_cmpt])]);
end
end
if (~isempty(Box_cmpt))
figure; 
%subplot(1,3,3); 
hold on;
cmptvec = Box_cmpt;
for ict = 1:length(cmptvec)
    icmpt = cmptvec(ict);
    Fac = [];
    for iboundary = 1:mymesh.Nboundary
        Fac = [Fac,mymesh.Fac_boundary_reorder{icmpt}{iboundary}];
    end
    h = trisurf(Fac',mymesh.Pts_cmpt_reorder{icmpt}(1,:),mymesh.Pts_cmpt_reorder{icmpt}(2,:),...
        mymesh.Pts_cmpt_reorder{icmpt}(3,:));
    set(h,'facealpha',0.9);
    axis equal;
    axis([xmin,xmax,ymin,ymax,zmin,zmax]); 
    view(3);
    title(['FE Mesh ECS cmpt: ',num2str([Box_cmpt])]);
end
end
