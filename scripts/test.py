from cmodule.build import bcrdist
import numpy as np

#array = bcrdist.load_bd_data("../data/data-test-heavy.csv", "../data/data-test-light.csv")
array = bcrdist.dsbcellarray()
print (array)
array.load_bd_data("../data/Hutchinson-BCRigh_DominantCDR3.csv", "../data/Hutchinson-BCRigl_DominantCDR3.csv")
print (array)
array.generate_dist_matrix();
dist, ids = array.dist_matrix();
print (dist)
print(ids)
