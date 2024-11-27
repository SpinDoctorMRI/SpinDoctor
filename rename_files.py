import os

if __name__=='__main__':
    for root,subdirs,files in os.walk(r'C:\Users\amcsween\SpinDoctor_saved_simul'):
        for file in files:
            if "FEM_" in file:
                os.rename(os.path.join(root,file),os.path.join(root,'FEM_lap_eig.mat'))
                # print(os.path.join(root,file))