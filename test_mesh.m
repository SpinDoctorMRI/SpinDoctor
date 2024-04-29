function flag = test_mesh(mesh_name)
%TEST_MESH This function tests if tetgen will produce a FE mesh from a
%surface input.
mesh_name = string(mesh_name);
if isfile(mesh_name)
    call_tetgen(mesh_name);
    [folder,name,~] = fileparts(mesh_name);
    fe_mesh = folder+"/"+name+".1";
    if isfile(fe_mesh+".ele")
        flag = true;
    else
        flag = false;
    end
    delete(fe_mesh+".ele");
    delete(fe_mesh+".smesh");
    delete(fe_mesh+".node");
    delete(fe_mesh+".face");
    delete(fe_mesh+".edge");
else
    warning('%s not found',mesh_name);
    flag = false;
end


end

