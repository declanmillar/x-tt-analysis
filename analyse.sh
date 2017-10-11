p="qq-tt-bbllvv_SM_13TeV_CT14LL.txt"
r="NuW"
b="2"
for f in $(ls -1 /scratch/dam1g09/zprime/qq-tt-bbmuevv_SM_13TeV_CT14LL_???.lhef.gz)
do
    c=$(echo $f | cut -d '/' -f 5)
    k=$(echo $c | cut -d '.' -f 1)
    s="${k}_pythia_delphes_${r}_b${b}"
    ./analyse.py "${k}_pythia_delphes_${r}_b${b}" -i "${k}_pythia_delphes.root" -p "$p" -r "$r" -b "$b"
done
