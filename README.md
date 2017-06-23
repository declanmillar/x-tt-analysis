# Branchus Readme

Previously know as `zprime-top-delphes-analysis`. Code to analyse top quark events in models containing *Z'* bosons. Reads `.root` output file from [Delphes](https://cp3.irmp.ucl.ac.be/projects/delphes).

The input is stored in a `std::tuple` with number of entries equal to the number of separate event files: `{string event_file, int proc_id}`

While an additional `std::tuple` is created for each individual process: `{string grid_file, int proc_id, int nfiles}`, where `nfiles` is the number of input files for the process with `proc_id`.

# Running the Program

* Run locally: `./analysis [options]`
* Submit job: `./analyse.py [options]`
