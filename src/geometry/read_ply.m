function [f,v,c] = read_ply(file)

fid = fopen(file);
tline = string(fgetl(fid));
vertex_properties= [];
face_properties = [];
while tline ~= "end_header"
    tline = string(fgetl(fid));
    if contains(tline, "element vertex")
        nvertex = str2num(erase(tline,"element vertex"));
        property_vertex_count = 0;
    elseif contains(tline, "element face")
        nface = str2num(erase(tline,"element face"));
        property_face_count = 0;
    elseif contains(tline, "property ")
        if ~exist('property_face_count')
            property_vertex_count = property_vertex_count +1;
            vertex_properties = [vertex_properties; erase(tline,"property ")];
        else 
            property_face_count = property_face_count +1;
            face_properties = [face_properties; erase(tline,"property ")];
        end
    end

end
V = fscanf(fid,"%f",[length(vertex_properties),nvertex]);
try
F = fscanf(fid,"%d %d %d %d "+repmat("%f",1,length(face_properties) - 1),[3+(length(face_properties)),nface]);
catch
F = fscanf(fid,"%d %d %d %d",[4,nface]);
end

fclose(fid);
v = V(1:3,:)';
f = F(2:4,:)' + 1;
c = V(4:end,:)';
