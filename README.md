# Acolyte

[Prophet](https://gitlab.cern.ch/demillar/prophet) -> [Delphes](https://gitlab.cern.ch/demillar/delphes) -> Acolyte -> [Disciple](https://gitlab.cern.ch/demillar/disciple)

Analyse Characteristics Of Leptonically decaYing Ttbar Events. Code to analyse top pair production events following processing by [Delphes](https://cp3.irmp.ucl.ac.be/projects/delphes).

## Installation

Dependencies: `boost`, `delphes`

## Important files

* `analysis`: complied `C++` executable. Use directly when running locally.
* `analyse.py`: run file to submit analysis as a batch job on `lxplus` or `iridis`.

## Running the Program

* Run locally: `./analysis [options]`
* Submit job: `./analyse.py handler_name [options]`

## Reminder 

The input is stored in a `std::vector` of `std::tuple` with number of entries equal to the number of separate event files: `{string event_file, int proc_id}`.
While an additional `std::vector` of `std::tuple` is created for each individual subprocess: `{string proc_file, int proc_id, int nfiles, double cross_section, double uncertainty, double weight}`, where `nfiles` is the number of input files for the process with `proc_id`.

## Examples commands

For example
```bash
./analysis -i "uu-X-tt-bbeevv_GLR-R-2.5-20pc_13TeV_CT14LL_000_pythia_delphes.root" -p "uu-X-tt-bbllvv_GLR-R-2.5-20pc_13TeV_CT14LL.txt" -r "NuW" -b 2
```

For example
```bash
./analyse.py "gg-tt-bbmumuvv_SM_13TeV_CT14LL_148_delphes_NuW_b1" -i "gg-tt-bbmumuvv_SM_13TeV_CT14LL_148_delphes.root" -r "NuW" -b 1
```

## Batch submission commands

Submit jobs using a wildcard.
```bash
r="NuW"
b="1"
for f in $(ls -1 /scratch/dam1g09/zprime/gg-tt-bbmumuvv_SM_13TeV_CT14LL_???.lhef.gz)
do
    c=$(echo $f | cut -d '/' -f 5)
    k=$(echo $c | cut -d '.' -f 1)
    ./analyse.py "${k}_delphes_${r}_b${b}" -i "${k}_delphes.root" -r "$r" -b "$b"
done
```

Submit jobs using a wildcard and a process file.
```bash
p="gg-tt-bbllvv_SM_13TeV_CT14LL.txt"
r="NuW"
b="1"
for f in $(ls -1 /scratch/dam1g09/zprime/gg-tt-bbmumuvv_SM_13TeV_CT14LL_???.lhef.gz)
do
    c=$(echo $f | cut -d '/' -f 5)
    k=$(echo $c | cut -d '.' -f 1)
    ./analyse.py "${k}_pythia_delphes_${r}_b${b}" -i "${k}_pythia_delphes.root" -p "$p" -r "$r" -b "$b"
done
```

Submit jobs using a wildcard with a period in the filename.
```bash
r="NuW"
b="2"
for f in $(ls -1 /scratch/dam1g09/zprime/??-AZX-tt-bbmumuvv_GLR-R-2.5_13TeV_CT14LL_???.lhef.gz)
do
    c=$(echo $f | cut -d '/' -f 5)
    k=$(echo $c | cut -d '.' -f 1)
    l=$(echo $c | cut -d '.' -f 2)
    m=$(echo $k.$l)
    ./analyse.py "${m}_pythia_delphes_${r}_b${b}" -i "${m}_pythia_delphes.root" -r "$r" -b "$b"
done
```

## Post-processing

Merge all files for each process - Z' only, 1 b-tag
```bash
hadd -f uu-X-tt-bbeevv_GLR-R-2.5-20pc_13TeV_CT14LL_pythia_delphes_NuW_b1.root uu-X-tt-bbeevv_GLR-R-2.5-20pc_13TeV_CT14LL_???_pythia_delphes_NuW_b1.root
hadd -f uu-X-tt-bbemuvv_GLR-R-2.5-20pc_13TeV_CT14LL_pythia_delphes_NuW_b1.root uu-X-tt-bbemuvv_GLR-R-2.5-20pc_13TeV_CT14LL_???_pythia_delphes_NuW_b1.root
hadd -f uu-X-tt-bbmuevv_GLR-R-2.5-20pc_13TeV_CT14LL_pythia_delphes_NuW_b1.root uu-X-tt-bbmuevv_GLR-R-2.5-20pc_13TeV_CT14LL_???_pythia_delphes_NuW_b1.root
hadd -f uu-X-tt-bbmumuvv_GLR-R-2.5-20pc_13TeV_CT14LL_pythia_delphes_NuW_b1.root uu-X-tt-bbmumuvv_GLR-R-2.5-20pc_13TeV_CT14LL_???_pythia_delphes_NuW_b1.root
```

Merge all files for each process - Z' only
```bash
hadd -f uu-X-tt-bbeevv_GLR-R-2.5-20pc_13TeV_CT14LL_pythia_delphes_NuW.root uu-X-tt-bbeevv_GLR-R-2.5-20pc_13TeV_CT14LL_???_pythia_delphes_NuW.root
hadd -f uu-X-tt-bbemuvv_GLR-R-2.5-20pc_13TeV_CT14LL_pythia_delphes_NuW.root uu-X-tt-bbemuvv_GLR-R-2.5-20pc_13TeV_CT14LL_???_pythia_delphes_NuW.root
hadd -f uu-X-tt-bbmuevv_GLR-R-2.5-20pc_13TeV_CT14LL_pythia_delphes_NuW.root uu-X-tt-bbmuevv_GLR-R-2.5-20pc_13TeV_CT14LL_???_pythia_delphes_NuW.root
hadd -f uu-X-tt-bbmumuvv_GLR-R-2.5-20pc_13TeV_CT14LL_pythia_delphes_NuW.root uu-X-tt-bbmumuvv_GLR-R-2.5-20pc_13TeV_CT14LL_???_pythia_delphes_NuW.root
```

Merge all files for each process - AZZ' only
```bash
hadd -f uu-AZX-tt-bbeevv_GLR-R-2.5-20pc_13TeV_CT14LL_pythia_delphes_NuW.root uu-AZX-tt-bbeevv_GLR-R-2.5-20pc_13TeV_CT14LL_???_pythia_delphes_NuW.root
hadd -f uu-AZX-tt-bbemuvv_GLR-R-2.5-20pc_13TeV_CT14LL_pythia_delphes_NuW.root uu-AZX-tt-bbemuvv_GLR-R-2.5-20pc_13TeV_CT14LL_???_pythia_delphes_NuW.root
hadd -f uu-AZX-tt-bbmuevv_GLR-R-2.5-20pc_13TeV_CT14LL_pythia_delphes_NuW.root uu-AZX-tt-bbmuevv_GLR-R-2.5-20pc_13TeV_CT14LL_???_pythia_delphes_NuW.root
hadd -f uu-AZX-tt-bbmumuvv_GLR-R-2.5-20pc_13TeV_CT14LL_pythia_delphes_NuW.root uu-AZX-tt-bbmumuvv_GLR-R-2.5-20pc_13TeV_CT14LL_???_pythia_delphes_NuW.root
```

Merge all subprocesses - AZZ' only
```bash
hadd -f uu-AZX-tt-bbllvv_GLR-R-2.5-20pc_13TeV_CT14LL_pythia_delphes_NuW.root uu-AZX-tt-bbeevv_GLR-R-2.5-20pc_13TeV_CT14LL_pythia_delphes_NuW.root uu-AZX-tt-bbemuvv_GLR-R-2.5-20pc_13TeV_CT14LL_pythia_delphes_NuW.root uu-AZX-tt-bbmuevv_GLR-R-2.5-20pc_13TeV_CT14LL_pythia_delphes_NuW.root uu-AZX-tt-bbmumuvv_GLR-R-2.5-20pc_13TeV_CT14LL_pythia_delphes_NuW.root
```

## Rename trees

```bash
root
TFile *file = new TFile("dd-AZX-tt-bbllvv.GLR-R-3.13TeV.CT14LL.root", "update")
RootTuple->SetName("events")
RootTuple->SetTitle("events")
gDirectory->Delete("RootTuple;1")
file->Write()
file->Close()
TTree* process = new TTree("process", "process")
double cross_section = 1.1911627439683250E-003
process->Branch("cross_section", &cross_section)
process->Fill()
file->Write()
file->Close()
```