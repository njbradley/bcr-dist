import bcrdist

arr = bcrdist.cell.array("../dekosky/donor3-limited-2")
arr.load_dekosky_data("../dekosky/donor3.csv");
#arr.save("../dekosky/save.tsv");
arr.generate_kpca_data()
