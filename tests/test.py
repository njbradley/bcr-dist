import bcrdist
import numpy as np
import sys
import os

#array = bcrdist.load_bd_data("../data/data-test-heavy.csv", "../data/data-test-light.csv")
arr = bcrdist.cell.array()
arr.load_bd_data("../data/data-test-heavy.csv", "../data/data-test-light.csv")
print (arr)
arr.save("../data/save.tsv")

arr2 = bcrdist.cell.array()
arr2.load("../data/save.tsv")

print (arr2)
