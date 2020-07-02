from cmodule.build import bcrdist

bcrdist.load_bd_data("../data/Hutchinson-BCRigh_DominantCDR3.csv")

bcrdist.save_dist_matrix("../data/hutch-dist.csv")
