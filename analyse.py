#!/usr/bin/env python

import os, StringIO, subprocess, sys, socket

# handler_name = sys.argv[1] + ".sh"
executable = "analysis"
# walltime = "02:00:00" # NuW
walltime = "00:00:30" # KIN
queue = "1nw"
nodes = "1"
ppn = "1"

# set directories
hostname = socket.gethostname()
run_directory = "."
data_directory = "."
if "lxplus" in hostname:
    run_directory = "/afs/cern.ch/user/d/demillar/acolyte"
    data_directory = "/afs/cern.ch/work/d/demillar/zprime"
elif "cyan" in hostname:
    run_directory = "/home/dam1g09/acolyte"
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
filen_ame = next(iterarg)
for arg in iterarg:
    argstring = argstring + " " + arg

handler_name = filen_ame + ".sh"
# logfile = filen_ame + ".log"
logfile = "/dev/null"

# print handler
handler = StringIO.StringIO()
if "lxplus" in hostname:
    print >> handler, "#!/bin/bash"
    print >> handler, "cd %s" % run_directory
    print >> handler, "source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6/setup.sh"
    print >> handler, "%s/%s %s" % (run_directory, executable, argstring)
if "cyan" or "blue" in hostname:
    print >> handler, "#PBS -N {}".format(filen_ame)
    print >> handler, "#PBS -l walltime={}".format(walltime)
    print >> handler, "#PBS -l nodes={}".format(nodes)
    print >> handler, "#!/bin/bash"
    print >> handler, "cd %s" % run_directory
    print >> handler, "source /home/dam1g09/.bash_profile"
    print >> handler, "%s/%s %s > %s" % (run_directory, executable, argstring, logfile)
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
elif "cyan" in hostname: subprocess.call('qsub %s/%s' % (run_directory, handler_name), shell = True)
