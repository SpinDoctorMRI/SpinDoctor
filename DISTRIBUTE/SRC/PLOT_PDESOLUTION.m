function PLOT_PDESOLUTION(mymesh,SOL,OUT_cmpts_index,ECS_cmpts_index,IN_cmpts_index,nindex,title_str)

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
            mymesh.Pts_cmpt_reorder{icmpt}(3,:),real(SOL{icmpt}(:,nindex)));
        set(h,'facealpha',0.6);set(h,'EdgeColor','none');
        axis equal;
        axis([xmin,xmax,ymin,ymax,zmin,zmax]); colorbar('eastoutside');
        view(3);
        title([title_str,' Inner cmpts: ',num2str(OUT_cmpts_index)]);
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
            mymesh.Pts_cmpt_reorder{icmpt}(3,:),real(SOL{icmpt}(:,nindex)));
        set(h,'facealpha',0.6);set(h,'EdgeColor','none');
        axis equal;
        axis([xmin,xmax,ymin,ymax,zmin,zmax]); colorbar('eastoutside');
        view(3);
        title([title_str,' Outer cmpts: ',num2str([IN_cmpts_index])]);
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
            mymesh.Pts_cmpt_reorder{icmpt}(3,:),real(SOL{icmpt}(:,nindex)));
        set(h,'facealpha',0.6);set(h,'EdgeColor','none');
        axis equal;
        axis([xmin,xmax,ymin,ymax,zmin,zmax]); c = colorbar('eastoutside');
        %set(c,'Position',[0.9251 0.1095 0.0381 0.8167]);
        view(3);
        title([title_str,' ECS cmpt: ',num2str([ECS_cmpts_index])]);
        %xlabel('x');
        %ylabel('y');
        view(-40,60);
        %view(90,90);
    end
end
