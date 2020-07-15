import bcrdist

arr = bcrdist.cell.array()

arr.load_bd_data("../data/Hutchinson-BCRigh_DominantCDR3.csv", "../data/Hutchinson-BCRigl_DominantCDR3.csv")

arr.generate_kpca_data()
