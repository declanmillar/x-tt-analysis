#!/usr/bin/env python

import os, StringIO, re, optparse, subprocess, sys, glob, socket, string

handler_name = sys.argv[2] + ".sh" 
executable = "analysis"
queue = "8nh"
walltime = "03:00:00"
run_directory = "."
data_directory = "."


# set directories
hostname = socket.gethostname()
run_directory = "."
data_directory = "."
if "lxplus" in hostname:
    run_directory = "/afs/cern.ch/user/d/demillar/zprime-top-analysis"
    data_directory = "/afs/cern.ch/work/d/demillar/zprime"
elif "cyan" in hostname:
    run_directory = "/home/dam1g09/zprime-top-analysis"
    data_directory = "/scratch/dam1g09/zprime"
else:
    sys.exit("Unrecognised hostname.")
if os.path.isdir(data_directory) is False:
    sys.exit("The target data directory '%s' does not exist" % data_directory)

argstring = ""
iterarg = iter(sys.argv)
next(iterarg)
for arg in iterarg:
    argstring = argstring + " " + arg

# print handler
handler = StringIO.StringIO()
print >> handler, "#!/bin/bash"
# if "lxplus" in hostname:
#     print >> handler, "source /afs/cern.ch/sw/lcg/external/gcc/4.8/x86_64-slc6/setup.sh"
if "cyan" in hostname:
    print >> handler, "module load gcc/4.9.1; source /local/software/cern/root_v6.06.06/bin/thisroot.sh"    
print >> handler, "cd %s" % run_directory
print >> handler, "%s/%s %s" % (run_directory, executable, argstring)

# write handler
try:
    with open('%s' % handler_name, 'w') as handler_file:
        handler_file.write(handler.getvalue())
    print "Handler file written to %s." % handler_name
except IOERROR:
    sys.exit("ERROR! Cannot write handler file.")

# run command
subprocess.call("chmod a+x %s" % handler_name, shell = True)
print "Submitting batch job."
if "lxplus" in hostname: subprocess.call('bsub -q %s %s/%s' % (queue, run_directory, handler_name), shell = True)
elif "cyan03" in hostname: subprocess.call('qsub -l walltime=%s %s/%s' % (walltime, run_directory, handler_name), shell = True)
