list="
driver_Ex4_T2_20e3_20e3_20e3
driver_Ex4_T2_40e3_20e3_80e3
driver_Ex4_T2_40e3_40e3_40e3
driver_Ex4_T2_80e3_80e3_80e3
driver_Ex4_T2_Inf_Inf_Inf
"

for f in ${list} 
do
    /Applications/MATLAB_R2018a.app/bin/matlab -nodesktop -nosplash -nodisplay -r $f &
done


