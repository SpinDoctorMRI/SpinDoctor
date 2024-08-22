meshes=$(<neuron_meshing_paper/cells.txt); tet="-pq1.2a0.5VCn";conda activate defaut
for mesh in $meshes; do python alphaSwc/get_mesh_quality.py $mesh $tet; done

mesh=mesh_files/selected/Homo-sapiens_-4203_BA46_1-3_a.CNG.ply; tet="-pq1.2a0.5VCn";conda activate defaut
python alphaSwc/get_mesh_quality.py $mesh $tet;