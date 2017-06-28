# Apsis

Previously know as `zprime-top-delphes-analysis`.

Code to analyse top quark events in models containing *Z'* bosons. Reads `.root` output file from [Delphes](https://cp3.irmp.ucl.ac.be/projects/delphes).

The input is stored in a `std::vector` of `std::tuple` with number of entries equal to the number of separate event files: `{string event_file, int proc_id}`.

While an additional `std::vector` of `std::tuple` is created for each individual subprocess: `{string proc_file, int proc_id, int nfiles, double cross_section, double uncertainty, double weight}`, where `nfiles` is the number of input files for the process with `proc_id`.

# Running the Program

* Run locally: `./analysis [options]`
* Submit job: `./analyse.py [options]`

# Important files

* `analysis`: complied `C++` executable. Use directly when running locally.
* `analyse.py`: run file to submit analysis as a batch job on `lxplus` or `iridis`.
