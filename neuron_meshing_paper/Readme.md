The neuron meshing paper present simulations from the following three cells


```bash
mesh_files/ultraliser/1-2-2-watertight.ply
mesh_files/ultraliser_modified/1-2-2.CNG_um.ply
mesh_files/selected/1-2-2.CNG.ply
```

which were meshed from the Ultraliser tool, modified Ultraliser and our method respectively. 

For each mesh individual simulations were run by running:

```matlab
driver_mf(mesh,"setup_comparison","-pq1.2a0.5O9VCn","3");
driver_btpde(mesh,"setup_comparison","-pq1.2a0.1O9VCn","1");
```

Convergence of solutions were established by saving the errors to  a .mat file with:

```
results = driver_get_PGSE_errors(mesh,"1","setup_comparison","neuron_meshing_paper/figures","3","0.5") ;
save("neuron_meshing_paper/"+cellname+".mat","results");
```

Directly compare volume weighted signals with the script 
```
cd neuron_meshing_paper; compare_ultraliser_signals;
```

Meshing and comparing convergence of individual compartments of the cell can be done by running:


```matlab
driver_mf_segmented(mesh,"setup_comparison","-pq1.2a0.5O9VCn","3");
driver_btpde(mesh,"setup_comparison","-pq1.2a0.1O9VCn","0");
results = driver_get_PGSE_errors(mesh,"1","setup_comparison","neuron_meshing_paper/figures","3","0.5");
```

Once convergence has been established, individual expermeints can be run by running:
```matlab
[results,femesh]= load_mf_segmented(mesh,"setup_comparison","-pq1.2a0.1O9VCn",'3');
```
or 
```matlab
[results,femesh]= load_mf(mesh,"setup_comparison","-pq1.2a0.1O9VCn",'3');
```