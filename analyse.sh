p="uu-AZX-tt-bbllvv_GLR-R-3_13TeV_CT14LL.txt"
r="NuW"
b="2"
for f in $(ls -1 /scratch/dam1g09/zprime/uu-AZX-tt-bbemuvv_GLR-R-2.5_13TeV_CT14LL_???.lhef.gz)
do
    c=$(echo $f | cut -d '/' -f 5)
    k=$(echo $c | cut -d '.' -f 1)
    ./analyse.py "${k}_pythia_delphes_${r}_b${b}" -i "${k}_pythia_delphes.root" -p "$p" -r "$r" -b "$b"
done


p="uu-AZX-tt-bbllvv_GLR-R-2.5-20pc_13TeV_CT14LL.txt"
r="NuW"
b="2"
for f in $(ls -1 /scratch/dam1g09/zprime/uu-AZX-tt-bb*vv_GLR-R-2.5-20pc_13TeV_CT14LL_???.lhef.gz)
do
    c=$(echo $f | cut -d '/' -f 5)
    k=$(echo $c | cut -d '.' -f 1)
    l=$(echo $c | cut -d '.' -f 2)
    m=$(echo $k.$l)
    ./analyse.py "${m}_pythia_delphes_${r}_b${b}" -i "${m}_pythia_delphes.root" -p "$p" -r "$r" -b "$b"
done

./analysis -i "dd-AZX-tt-bbmumuvv_GLR-R-2.5_13TeV_CT14LL_026_pythia_delphes.root" -p "dd-AZX-tt-bbllvv_GLR-R-2.5_13TeV_CT14LL.txt" -r "KIN" -b 2
