p="gg-tt-bbllvv_SM_13TeV_CT14LL.txt"
for f in $(ls -1 /scratch/dam1g09/zprime/gg-tt-bbmuevv_SM_13TeV_CT14LL_???.lhef.gz)
do
    c=$(echo $f | cut -d '/' -f 5)
    k=$(echo $c | cut -d '.' -f 1)

    ./analyse.py "${k}" -i "${k}_pythia_delphes.root" -p "$p" -r "NuW"
done

./analyse.py "gg-tt-bbmuevv_SM_13TeV_CT14LL_001" -i "/scratch/dam1g09/zprime/gg-tt-bbmuevv_SM_13TeV_CT14LL_001.lhef.gz" -p "gg-tt-bbllvv_SM_13TeV_CT14LL.txt" -r "NuW"
