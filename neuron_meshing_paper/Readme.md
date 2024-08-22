To download data from Margaret.

$list = Get-Content .\neuron_meshing_paper\cells.txt -Raw;

$list=$list -Split "\n";

foreach ($cell in $list){

$sourceDirec = -join('amcsween@margaret.saclay.inria.fr:/scratch/amcsween/SpinDoctor/',$cell);
scp -r $sourceDirec  $cell;


$name = Split-Path -Path $cell -Leaf -Resolve;

$parts = $name.Split("."); $L = $parts.Length - 2;

$cellname = $parts[0..$L] -Join ".";

echo $cellname;
$Direc = Split-Path -Path $cell;
$sourceDirec = -join('amcsween@margaret.saclay.inria.fr:/scratch/amcsween/SpinDoctor/',$Direc,'/',$cellname,'_ply_dir');
$targetDirec = -join($Direc,'/',$cellname,'_ply_dir');

scp -r $sourceDirec $targetDirec


}

foreach ($cell in $list){
$name = Split-Path -Path $cell -Leaf -Resolve;
$parts = $name.Split("."); $L = $parts.Length - 2;
$cellname = $parts[0..$L] -Join ".";
echo $cellname;

$sourceDirec = -join('amcsween@margaret.saclay.inria.fr:/scratch/amcsween/SpinDoctor/saved_simul/',$cellname,'_tet*');

echo $sourceDirec;
scp -r $sourceDirec ./saved_simul;
}



