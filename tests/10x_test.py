import bcrdist
import sys

if (len(sys.argv) < 1):
    print ("Give the filename as the first argument")

infile = sys.argv[1]

arr = bcrdist.load10x(infile);
arr.savePCs()
arr.umapplot(colorby = "clonotype")
#arr.tsneplot(colorby = "clonotype")
