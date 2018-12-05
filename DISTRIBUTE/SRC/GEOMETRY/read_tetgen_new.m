function [Pts_cmpt_reorder,Ele_cmpt_reorder,Pts_ind,Pts_boundary_reorder,Fac_boundary_reorder,...
    Nboundary,Ncmpt] = read_tetgen_new(fname,para_deform)



disp(['Reading from file ',fname]);

filename = [fname,'.node'];
fid = fopen(filename,'r');
if (fid ~= -1)
    Nnode = fscanf(fid,'%d',[4,1]);
    Rnode = fscanf(fid, '%f', [4,inf]);
end
fclose(fid);

Pts_all = Rnode(2:end,:);


filename = [fname,'.ele'];
fid = fopen(filename,'r');
if (fid ~= -1)
    Nele = fscanf(fid,'%d',[3,1]);
    Rele = fscanf(fid, '%f', [6,inf]);
end
fclose(fid);

Ele_all = Rele(2:5,:);
Ele_attrib = Rele(6,:);

filename = [fname,'.face'];
fid = fopen(filename,'r');
if (fid ~= -1)
    Nface = fscanf(fid,'%d',[2,1]);
    Rface = fscanf(fid, '%f', [5,inf]);
end
fclose(fid);

Fac_all = Rface(2:4,:);
Fac_attrib = Rface(5,:);

temp=[];
Cmpt_attrib = unique(Ele_attrib);
Ncmpt = length(Cmpt_attrib);

disp(['Separating into ',num2str(Ncmpt), ' compartments ']);

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
Nboundary = length(Boundary_attrib);

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


