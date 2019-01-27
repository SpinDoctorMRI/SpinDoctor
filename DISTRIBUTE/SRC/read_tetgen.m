function [mymesh,cmpts_bdys_mat] = read_tetgen(fname,para_deform,Ncmpt,Nboundary)

% create FE mesh on canonical configuration; bend and twist the FE mesh nodes by analytical transformation
% 
% Input:
%     1. fname
%     2. para_deform
%     3. Ncmpt
%     4. Nboundary
%     
% Output:
%     1. mymesh is a structure with 10 elements:
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
%     2. cmpts_bdys_mat


% [Pts_cmpt_reorder,Ele_cmpt_reorder,Pts_ind,Pts_boundary_reorder,Fac_boundary_reorder, Nboundary,Ncmpt]

disp(['Reading from Tetgen FE mesh from ',fname]);

filename = [fname,'.node'];
fid = fopen(filename,'r');
if (fid ~= -1)
    Nnode = fscanf(fid,'%d',[4,1]);
    Rnode = fscanf(fid, '%f', [4,inf]);
end
fclose(fid);
mymesh.Nnode = Nnode;
Pts_all = Rnode(2:end,:);

filename = [fname,'.ele'];
fid = fopen(filename,'r');
if (fid ~= -1)
    Nele = fscanf(fid,'%d',[3,1]);
    Rele = fscanf(fid, '%f', [6,inf]);
end
fclose(fid);
mymesh.Nele = Nele;

Ele_all = Rele(2:5,:);
Ele_attrib = Rele(6,:);

filename = [fname,'.face'];
fid = fopen(filename,'r');
if (fid ~= -1)
    Nface = fscanf(fid,'%d',[2,1]);
    Rface = fscanf(fid, '%f', [5,inf]);
end
fclose(fid);
mymesh.Nface = Nface;

Fac_all = Rface(2:4,:);
Fac_attrib = Rface(5,:);

temp = [];
Cmpt_attrib = unique(Ele_attrib);

Ncmpt_tmp = length(Cmpt_attrib);

disp(['Separating FE mesh into ',num2str(Ncmpt_tmp), ' compartments ']);

if (Ncmpt_tmp ~= Ncmpt)
    disp(['FE mesh not good, use smaller hmax or change surf triangulation']);
    mymesh = [];
    cmpts_bdys_mat = [];
else

for icmpt = 1:Ncmpt
    jj = find(Ele_attrib == Cmpt_attrib(icmpt));
    Ele_cmpt{icmpt} = Ele_all(:,jj); 
    Pts_ind{icmpt} = unique(Ele_cmpt{icmpt}(:));
    Pts_cmpt_reorder{icmpt} = Pts_all(:,Pts_ind{icmpt});
end

if (~isempty(find(para_deform ~= 0)))
    for icmpt = 1:Ncmpt
        Pts_cmpt_reorder{icmpt} = deform_domain(Pts_cmpt_reorder{icmpt},para_deform);
    end
end

Boundary_attrib = unique(Fac_attrib);
Nboundary_tmp = length(Boundary_attrib);

disp(['Separating FE mesh with ',num2str(Nboundary_tmp), ' boundaries ']);

if (Nboundary_tmp ~= Nboundary)
    disp(['FE mesh not good, use smaller hmax or change surf triangulation']);
    mymesh = [];
    cmpts_bdys_mat = [];
else

for iboundary = 1:Nboundary
    jj = find(Fac_attrib == Boundary_attrib(iboundary));
    Fac_boundary{iboundary} = Fac_all(:,jj);
    Pts_boundary{iboundary} = unique(Fac_all(:,jj));
end

for icmpt = 1:Ncmpt
    Ele_cmpt_reorder{icmpt} = Ele_cmpt{icmpt};
    for ii = 1:length(Pts_ind{icmpt})
        jj = find(Ele_cmpt{icmpt}==Pts_ind{icmpt}(ii));
        Ele_cmpt_reorder{icmpt}(jj) = ii;
    end
    for iboundary = 1:Nboundary
        onbd = 1;
        for ib = 1:size(Fac_boundary{iboundary},2)
            ind = find(Fac_boundary{iboundary}(ib)==Pts_ind{icmpt});
            if (isempty(ind))
                onbd = 0;
            end
        end
        if (onbd == 1)
            Fac_boundary_reorder{icmpt}{iboundary} = Fac_boundary{iboundary};
            for ii = 1:length(Pts_ind{icmpt})
                jj = find(Fac_boundary{iboundary}==Pts_ind{icmpt}(ii));
                Fac_boundary_reorder{icmpt}{iboundary}(jj) = ii;
            end
            Pts_boundary_reorder{icmpt}{iboundary} = Pts_boundary{iboundary};
            for ii = 1:length(Pts_ind{icmpt})                
                jj = find(Pts_boundary{iboundary}==Pts_ind{icmpt}(ii));
                Pts_boundary_reorder{icmpt}{iboundary}(jj) = ii;                
            end
        else
            Fac_boundary_reorder{icmpt}{iboundary} = [];
            Pts_boundary_reorder{icmpt}{iboundary} = [];
        end
    end
end

mymesh.Pts_cmpt_reorder = Pts_cmpt_reorder;
mymesh.Ele_cmpt_reorder = Ele_cmpt_reorder;
mymesh.Pts_ind = Pts_ind;
mymesh.Pts_boundary_reorder = Pts_boundary_reorder;
mymesh.Fac_boundary_reorder = Fac_boundary_reorder;
mymesh.Nboundary = Nboundary;
mymesh.Ncmpt = Ncmpt;

cmpts_bdys_mat = zeros(mymesh.Ncmpt,mymesh.Nboundary);
for icmpt = 1:mymesh.Ncmpt
    for ibd = 1:mymesh.Nboundary
        cmpts_bdys_mat(icmpt,ibd)=~isempty(mymesh.Pts_boundary_reorder{icmpt}{ibd});
    end
end
% %%%%%%%%%%%%%%%%%%% Start checking the mesh quality %%%%%%%%%%%%%%%%%%%%%%%
% aspect_ratio_lim = 0.05; % the worst -> [0, 1] <- the best
% Qmesh=cell(1,mymesh.Ncmpt);
% hmax=0;
% for icmpt = 1:mymesh.Ncmpt
%     Qmesh{icmpt} = mesh_quality(mymesh.Pts_cmpt_reorder{icmpt},mymesh.Ele_cmpt_reorder{icmpt});
%     hmax = max(max(Qmesh{icmpt}.hout),hmax);
% %    if Qmesh{icmpt}.quality1(1)<aspect_ratio_lim
%         disp(['  Compartment ',num2str(icmpt),' - FE mesh with minimum aspect ratio of ',num2str(Qmesh{icmpt}.quality1(1),'%.1e')]);
% %    end;
% end;
% %%%%%%%%%%%%%%%%%%% End of checking the mesh quality %%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Start checking the mesh quality %%%%%%%%%%%%%%%%%%%%%%%
aspect_ratio_lim = 0.05; % the worst -> [0, 1] <- the best
Qmesh = cell(1,mymesh.Ncmpt);
hmax = 0;
for icmpt = 1:mymesh.Ncmpt
    Qmesh{icmpt}=mesh_quality(mymesh.Pts_cmpt_reorder{icmpt},mymesh.Ele_cmpt_reorder{icmpt}) ;
    hmax = max(max(Qmesh{icmpt}.hout),hmax);
    disp(['  Compartment ',num2str(icmpt),' - FE mesh with minimum aspect ratio of ',num2str(Qmesh{icmpt}.quality(1),'%.1e')]);
end
%%%%%%%%%%%%%%%%%%% End of checking the mesh quality %%%%%%%%%%%%%%%%%%%%%%

end
end