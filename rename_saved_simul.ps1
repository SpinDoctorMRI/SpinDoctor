$add = "mf_D_out0.002_permea0_af280b\relax_oInf_density_o1_36acdf";
$source = pwd;
Get-ChildItem -Path C:\Users\amcsween\SpinDoctor_saved_simul -Directory | ForEach-Object {
    $old_folder = $_.FullName;
    # $folder_name = $_.BaseName;
    # $parts = $folder_name -Split('_tet');
    # $cell = $parts[0]; $tet = $parts[1];
    # $new_folder_name = $cell + "_ply_tet" + $tet + "_no_ecs";
    # echo $new_folder_name;
    $new_folder_name = $old_folder + "\" + $add;
    # echo $new_folder_name;
    # Rename-Item -Path $old_folder -NewName $new_folder_name
    # $new_folder = $old_folder + $add;
    # echo $new_folder
    # echo $new_folder_name;
    # if (Test-Path $new_folder_name){
        echo $new_folder_name;
    #     rm $new_folder_name;
    # }
    cd $old_folder;
    # mkdir $new_folder_name;
    Get-ChildItem -Path $old_folder -Recurse | ForEach-Object {
        $filename = Resolve-Path $_.FullName -Relative;
        $check_filename = $filename.Replace('mf_D_out0.002_permea0_af280b\relax_oInf_density_o1_36acdf','')
        if (Test-Path $check_filename){
        echo $filename;
        $new_filename =  $add + "\"+ $filename;
        $new_filename = $new_filename.Replace('\.','');
        echo $new_filename;
        # mv $filename $new_filename
        Rename-Item -Path $filename -NewName $new_filename
        } 
        else{
            echo "Already moved";
        }

        # $new_filename =  $add + "\"+ $filename;
        # $new_filename = $new_filename.Replace('\.','');
        # echo $new_filename;
        # mv $filename $new_filename
        # Rename-Item -Path $filename -NewName $new_filename

    }
}
Set-Location $source;

