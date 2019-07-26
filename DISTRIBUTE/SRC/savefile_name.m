function [fname] = savefile_name(SDELTA,BDELTA,bvalue,gdir,icmpt)

seq_name = ['d',num2str(SDELTA),'_D',num2str(BDELTA)];
bvalue_name = ['b',num2str(bvalue)];
if (length(gdir) == 3)
    gdir_name = ['gdir',num2str(gdir(1)),'_',num2str(gdir(2)),'_',num2str(gdir(3))];
else
    gdir_name = ['gdir',num2str(gdir(1))];
end
cmpt_name = ['cmpt',num2str(icmpt)];
fname = [seq_name,'_',bvalue_name,'_',gdir_name,'_',cmpt_name];