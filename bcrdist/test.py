import cell
from cmodule.build import bcrdist
import numpy as np

#array = bcrdist.load_bd_data("../data/data-test-heavy.csv", "../data/data-test-light.csv")
arr = cell.array()

arr.load_bd_data("../data/data-test-heavy.csv", "../data/data-test-light.csv")



print (arr)
