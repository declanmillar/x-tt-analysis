#!/usr/bin/env python

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

fname = "SM_gg-tt-bbllvv_2-4_5x10M.lhe"

print "%s has %i lines." % (fname, file_len(fname))

