function  [name] = savefile_name(SD,BD,bval,ng_total,i)
% create the filename string for the saving data
    name = ['d',num2str(SD),'_','D',num2str(BD),'_','b',num2str(bval),'_gdir',num2str(ng_total),'_cmpt',num2str(i),'gdir2D']
end