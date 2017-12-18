#!/usr/bin/env python

import os, StringIO, subprocess, sys, socket

# handler_name = sys.argv[1] + ".sh"
executable = "analysis"
walltime = "02:00:00"
queue = "1nw"

# set directories
hostname = socket.gethostname()
run_directory = "."
data_directory = "."
if "lxplus" in hostname:
    run_directory = "/afs/cern.ch/user/d/demillar/apsis"
    data_directory = "/afs/cern.ch/work/d/demillar/zprime"
elif "cyan" in hostname:
    run_directory = "/home/dam1g09/apsis"
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
filename = next(iterarg)
for arg in iterarg:
    argstring = argstring + " " + arg

handler_name = filename + ".sh"
# logfile = filename + ".log"
logfile = "/dev/null"

# print handler
handler = StringIO.StringIO()
print >> handler, "#!/bin/bash"
if "lxplus" in hostname:
    print >> handler, "source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6/setup.sh"
    print >> handler, "cd %s" % run_directory
    print >> handler, "%s/%s %s" % (run_directory, executable, argstring)
if "cyan" or "blue" in hostname:
    print >> handler, "source /home/dam1g09/.bash_profile"
    print >> handler, "echo 'changing to run directory ...'"
    print >> handler, "cd %s" % run_directory
    print >> handler, "echo 'running code ...'"
    print >> handler, "%s/%s %s > %s" % (run_directory, executable, argstring, logfile)
    print >> handler, 'rm -- "$0"'
else:
    sys.exit("ERROR: Unrecognised hostname")

# write handler
try:
    with open('%s' % handler_name, 'w') as handler_file:
        handler_file.write(handler.getvalue())
    print "handler file: %s" % handler_name
    print "log file:     %s" % logfile
except IOERROR:
    sys.exit("ERROR: Cannot write handler file.")

# run command
subprocess.call("chmod a+x %s" % handler_name, shell = True)
print "submitting batch job ..."
if "lxplus" in hostname: subprocess.call('bsub -q %s %s/%s' % (queue, run_directory, handler_name), shell = True)
elif "cyan03" in hostname: subprocess.call('qsub -l walltime=%s %s/%s' % (walltime, run_directory, handler_name), shell = True)
