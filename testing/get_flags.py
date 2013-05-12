#!/usr/bin/env python

"""Get the on switches for a merger tree flag integer.

Usage: get_flags.py <flag>
"""

from docopt import docopt

__author__ = "Simon Mutch"
__date__   = "25 Mar 2013"

args = docopt(__doc__)

switches = {
    "TREE_CASE_SIMPLE": 1,
    "TREE_CASE_MAIN_PROGENITOR": 2,
    "TREE_CASE_MERGER": 4,
    "TREE_CASE_DROPPED": 8,
    "TREE_CASE_STRAYED": 16,
    "TREE_CASE_SPUTTERED": 32,
    "TREE_CASE_BRIDGED": 64,
    "TREE_CASE_EMERGED_CANDIDATE": 128,
    "TREE_CASE_FOUND": 256,
    "TREE_CASE_NO_PROGENITORS": 512,
    "TREE_CASE_FRAGMENTED_LOST": 1024,
    "TREE_CASE_FRAGMENTED_RETURNED": 2048,
    "TREE_CASE_FRAGMENTED_EXCHANGED": 4096,
    "TREE_CASE_MATCHED_TO_BRIDGE": 8192,
    "TREE_CASE_BRIDGE_DEFAULT": 16384,
    "TREE_CASE_GHOST": 32768,
}

flag = int(args['<flag>'])
print "flag=%d matches:" % flag

for s,v in switches.iteritems():
    if (flag & v)==v:
        print s
