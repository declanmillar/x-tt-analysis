#!/usr/bin/env python

import os, StringIO, subprocess, sys, socket

executable = "analysis"
walltime = "00:30:00"
queue = "1nw"
nodes = "1"
ppn = "1"

# set directories
run_directory = "."
data_directory = "."
hostname = socket.gethostname()
if "lxplus" in hostname:
    run_directory = "/afs/cern.ch/user/d/demillar/acolyte/build"
    data_directory = "/afs/cern.ch/work/d/demillar/zprime"
elif "cyan" in hostname:
    run_directory = "/home/dam1g09/acolyte/build"
    data_directory = "/scratch/dam1g09/zprime"
else:
    sys.exit("Error: unrecognised hostname")

# check directories exist
if not os.path.isdir(run_directory):
    sys.exit("Error: No directory %s" % run_directory)
if not os.path.isdir(data_directory):
    sys.exit("Error: No directory %s" % data_directory)

# collect command line arguments
argstring = ""
iterarg = iter(sys.argv)
next(iterarg)
file_name = next(iterarg)
for arg in iterarg:
    argstring = argstring + " " + arg

# print handler
handler = StringIO.StringIO()
if "lxplus" in hostname:
    print >> handler, "#!/bin/bash"
    print >> handler, "cd %s" % run_directory
    print >> handler, "source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6/setup.sh"
    print >> handler, "%s/%s %s" % (run_directory, executable, argstring)
if "cyan" or "blue" in hostname:
    print >> handler, "#PBS -N {}".format(file_name)
    print >> handler, "#PBS -l walltime={}".format(walltime)
    print >> handler, "#PBS -l nodes={}".format(nodes)
    print >> handler, "#!/bin/bash"
    print >> handler, "cd %s" % run_directory
    print >> handler, "source /home/dam1g09/.bash_profile"
    print >> handler, "%s/%s %s" % (run_directory, executable, argstring)
else:
    sys.exit("ERROR: Unrecognised hostname")

# write handler
handler_name = file_name + ".sh"
try:
    with open('%s/%s' % (run_directory, handler_name), 'w') as handler_file:
        handler_file.write(handler.getvalue())
    print "run script: %s" % handler_name
except IOERROR:
    sys.exit("ERROR: Cannot write run script.")

# run command
subprocess.call("chmod a+x %s/%s" % (run_directory,handler_name), shell = True)
if "lxplus" in hostname:
    subprocess.call('bsub -q %s %s/%s' % (queue, run_directory, handler_name), shell = True)
elif "cyan" in hostname:
    subprocess.call('qsub %s/%s' % (run_directory, handler_name), shell = True)
else:
    sys.exit("Error: unrecognised hostname")
