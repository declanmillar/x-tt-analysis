#!/usr/bin/env python

import os, StringIO, subprocess, sys, socket, argparse

# handler_name = sys.argv[1] + ".sh"
executable = "analysis"
# walltime = "60:00:00"
# queue = "1nw"

parser = argparse.ArgumentParser(description="generate ttbar events")
parser.add_argument("-a", "--arguments", help = "argument list to be passed to C++ program", default = "")
parser.add_argument("-w", "--walltime", help = "cut off time for job completion", default = "60:00:00")
parser.add_argument("-Q", "--queue", help = "queue to join", default = "1nw")
parser.add_argument("-P", "--parallel", help = "submit each input file as an individual job", default = True,  action = "store_false")
parser.add_argument("-m", "--model",          default = "SM")
parser.add_argument("-g", "--gg",             default = False, action = "store_true")
parser.add_argument("-q", "--qq",             default = False, action = "store_true")
parser.add_argument("-u", "--uu",             default = False, action = "store_true")
parser.add_argument("-d", "--dd",             default = False, action = "store_true")
parser.add_argument("-f", "--final_state",    default = "tt-bbeevv")
parser.add_argument("-E", "--energy",         default = 13)
parser.add_argument("-L", "--luminosity",     default = -1)
parser.add_argument("-b", "--minimumBtags",   default = 2)
parser.add_argument("-o", "--options",        default = "")
parser.add_argument("-r", "--reconstruction", default = "TRN")
parser.add_argument("-t", "--tag",            default = "")
parser.add_argument("-s", "--slice",          default = False, action = "store_true")
parser.add_argument("-i", "--inputfilename",  default = "")
parser.add_argument("-p", "--processfilename",default = "")
args = parser.parse_args()

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

handler_name = args.inputfilename + ".sh"
logfile = args.inputfilename + ".log"

# print handler
handler = StringIO.StringIO()
print >> handler, "#!/bin/bash"
if "lxplus" in hostname:
    print >> handler, "source /afs/cern.ch/sw/lcg/external/gcc/4.9/x86_64-slc6/setup.sh"
    print >> handler, "cd %s" % run_directory
    print >> handler, "%s/%s %s" % (run_directory, executable, args.arguments)
if "cyan" or "blue" in hostname:
    print >> handler, "source /home/dam1g09/.bash_profile"
    print >> handler, "echo 'changing to run directory ...'"
    print >> handler, "cd %s" % run_directory
    print >> handler, "echo 'running code ...'"
    print >> handler, "%s/%s %s > %s/%s" % (run_directory, executable, args.arguments, run_directory, logfile)
else:
    sys.exit("Error: Unrecognised hostname")

# write handler
try:
    with open('%s' % handler_name, 'w') as handler_file:
        handler_file.write(handler.getvalue())
    print "Handler file written to %s." % handler_name
    print "Log written to %s." % logfile
except IOERROR:
    sys.exit("ERROR! Cannot write handler file.")

# run command
subprocess.call("chmod a+x %s" % handler_name, shell = True)
print "Submitting batch job."
if "lxplus" in hostname: subprocess.call('bsub -q %s %s/%s' % (queue, run_directory, handler_name), shell = True)
elif "cyan03" in hostname: subprocess.call('qsub -l walltime=%s %s/%s' % (walltime, run_directory, handler_name), shell = True)
