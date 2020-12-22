import bcrdist

arr = bcrdist.cellarray("test")
print(arr)
arr.append(("cell", "clono", "AA", "AAA", "AAAA", "GG", "GGG", "GGGG"))
print(arr)
arr.append(("hi","clono","IGLV3-21*01", "CQVWDTSGDHHVF", "IGLV3-21*01", "CCCQQ"))
print(arr)
arr.save()

