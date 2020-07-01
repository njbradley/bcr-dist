from cmodule.build import bcrdist

bcrdist.load_bd_data("../data/data-test.csv")

bcrdist.save_dist_matrix("../data/dist.csv")
