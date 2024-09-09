To run a simulation on a single cell run:

```
driver_mf(mesh_path,setup_file,tetgen_options,ls);
```
The parameters already computed are tetgen_options = "-pq1.2a1.0VCn" and ls="1". Choose the diffusion encoding gradient by changing the setup.gradient.sequences in the setup_file where the parameters for the eigenfunction decomposition can also be set.

For the current data set, we do not support runinng driver_mf_segmented, instead change the mesh_path to the appropriate component.

Signals can be loaded with:

```matlab
[results,femesh] = load_mf(mesh_path,setup_file,tetgen_options,ls);
```
and then saved for analysis as a .mat file.

To check convergence it is neccessary to first run:

```matlab
driver_btpde(mesh_path,"setup_morez_ref_sol", "-pq1.2a0.1VCn","1");
```
before comparing results with
```matlab
driver_get_morez_errors(mesh_path,"1","1',"1","figures_error");  
```
to save the resulting figures to figures_error.
