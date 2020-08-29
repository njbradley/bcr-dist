import bcrdist
import sys

if (len(sys.argv) < 1):
    print ("Give the filename as the first argument")

infile = sys.argv[1]

arr = bcrdist.cell.load10x(infile);
arr.savePCs()
arr.tsneplot()
arr.umapplot()
