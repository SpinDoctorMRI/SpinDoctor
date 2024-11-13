$list = Get-Content .\cells_human_um.txt -Raw;

$list=$list -Split "\n";

$output_dir = "C:\Users\amcsween\SpinDoctor_saved_simul";

foreach ($cell in $list){

$path = $cell.Split('/');
$name = $path[$path.Length - 1];

$parts = $name.Split("."); 

$L = $parts.Length - 2;
$cellname = $parts[0..$L] -Join ".";


$sourceDirec = -join('amcsween@margaret.saclay.inria.fr:/scratch/amcsween/SpinDoctor_Neuron_Meshing_Paper/saved_simul/',$cellname,'_tet-pq1.2a0.05O9VCn');
$targetDirec = -join($output_dir,'/',$cellname,'_tet-pq1.2a0.05O9VCn');

if (Test-Path -Path $targetDirec){
    echo 'is path';
    rm $targetDirec
}



echo $sourceDirec;

scp -r $sourceDirec $output_dir;

}