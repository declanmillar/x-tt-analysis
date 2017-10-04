p="gg-tt-bbllvv_SM_13TeV_CT14LL.txt"
r="NuW"
b="0"
for f in $(ls -1 /scratch/dam1g09/zprime/gg-tt-bbmuevv_SM_13TeV_CT14LL_???.lhef.gz)
do
    c=$(echo $f | cut -d '/' -f 5)
    k=$(echo $c | cut -d '.' -f 1)
    s="${k}_pythia_delphes_${r}_b${b}"

    ./analyse.py "${k}" -i "${k}_pythia_delphes.root" -p "$p" -r "$r" -b "$b"
done

./analyse.py "gg-tt-bbmuevv_SM_13TeV_CT14LL_001" -i "gg-tt-bbmuevv_SM_13TeV_CT14LL_001_pythia_delphes.root" -p "gg-tt-bbllvv_SM_13TeV_CT14LL.txt" -r "NuW"
