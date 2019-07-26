function [Pts_new,Fac_new,meshlimits] = FEMESH2PLY(mymesh,fname_tetgen)

for icmpt = 1:mymesh.Ncmpt
                 
    Fac = [];
    Pts_index = [];
    for iboundary = 1:mymesh.Nboundary
        Fac = [Fac,mymesh.Fac_boundary_reorder{icmpt}{iboundary}]; 
        Pts_index = [Pts_index;mymesh.Pts_boundary_reorder{icmpt}{iboundary}]; 
    end
    
    Pts_index = unique(Pts_index);
    npt = length(Pts_index)
    
    Pts_new{icmpt} = zeros(3,npt);
    for jj = 1:3
        Pts_new{icmpt}(jj,:) = mymesh.Pts_cmpt_reorder{icmpt}(jj,Pts_index);
    end
    
    nnodes = size(mymesh.Pts_cmpt_reorder{icmpt},2);    
    nodesindex = zeros(nnodes,1); 
    for ipt = 1:npt
        nodesindex(Pts_index(ipt),1) = ipt; 
    end
    
    nfac = size(Fac,2);
    Fac_new{icmpt} = zeros(size(Fac));
    
    for ifac = 1:nfac
        for jj = 1:3
            index_new = nodesindex(Fac(jj,ifac),1);
            Fac_new{icmpt}(jj,ifac) = index_new;
        end
    end
    
    figure; 
    h = trisurf(Fac_new{icmpt}(1:3,1:nfac)',Pts_new{icmpt}(1,1:npt),Pts_new{icmpt}(2,1:npt),...
        Pts_new{icmpt}(3,1:npt),'facealpha',0.1);
    axis equal;
    xlabel('x');
    ylabel('y');
    zlabel('z');
        
    filename = [fname_tetgen,'_cmpt',num2str(icmpt),'.ply']
    fid = fopen(filename,'w');
    
    ply_header = {'ply',...
        'comment closed surface',...
        'format ascii 1.0',...
        'comment created by SpinDoctor',...
        ['element vertex ',num2str(npt)],...
        'property float x',...
        'property float y',...
        'property float z',...
        ['element face ',num2str(nfac)],...
        'property list uchar short vertex_indices',...
        'end_header'};
    
    if (fid ~= -1)
        for iii = 1:11
            fprintf(fid, '%s\n', ply_header{iii});
        end  
    end
    
    for ipt = 1:npt
        fprintf(fid, '%26.16e %26.16e %26.16e\n', (1e-6)*Pts_new{icmpt}(1:3,ipt));
        %fprintf(fid, '%26.16e %26.16e %26.16e\n', Pts_new{icmpt}(1:3,ipt));
    end
    
    meshlimits{icmpt} = [min(Pts_new{icmpt}(1,:)),max(Pts_new{icmpt}(1,:)),...
        min(Pts_new{icmpt}(2,:)),max(Pts_new{icmpt}(2,:)),...
        min(Pts_new{icmpt}(3,:)),max(Pts_new{icmpt}(3,:))];
    
    for ii = 1:nfac
        fprintf(fid, '%d ', 3);
        for jj = 1:3
            fprintf(fid, '%d ', Fac_new{icmpt}(jj,ii)-1);
        end
        fprintf(fid, '\n');
    end
        
    fclose(fid);
end