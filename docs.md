# Bcrdist docs

## bcrdist.cell.array
This is the main class that everything in bcrdist is run on. It is an array of bcells, and it can hold either single or paired bcell sequences.

## bcrdist.cell.load10x(string path)
This method creates a bcrdist.cell.array object and fills it with the cells in the input path. The input path is expected to be a 10x contig annotations file.
This method is shorthand for
```python
arr = bcrdist.cell.array()
arr.load10X(path)
```

## bcrdist.cell.loadBD(string heavypath, [string lightpath]) -> None
This method, like load10X, creates a cell array and populates it with the bd data files. BD data has on eheavy and one light file, so if the data you are using is paired, you have to put in both filenames
This method is shorthand for
```python
arr = bcrdist.cell.array()
arr.loadBD(path)
```

## bcrdist.cell.array([string name]) -> None
The constructor for bcrdist.cell.array can take a name parameter. This name will be stored and used in other functions. If a name is not specified, it remains uuninitialized until data is loaded.
When the array is created, it looks to see if a save file exists for the specified name. If it does, that file is loaded.

## bcrdist.cell.array.load10x(string path) -> None
This method loads bcells from a 10x type data file, and appends the cells onto this cell array.

## bcrdist.cell.array.loadBD(string heavypath, [string lightpath]) -> None
This method loads bcells from one or two bd data files and appends the cells to the cell array. An error is thrown if the loaded cells are not the same type (paired/not paired) as the cells in the array.

## bcrdist.cell.array.generate_distmatrix() -> None
This method generates a distance matrix with the current cells, and writes it to a file. The file path is the name + ".bcrdistmatrix.tsv"
Unlike bcrdist.cell.array.distmatrix(), this method will generate a matrix even if there is already one saved to file, so it is useful to regenerate the distance matrix. This can also be accomplished by deleting the distance matrix file.

## bcrdist.cell.array.distmatrix() -> ( list id_list, np.array distmatrix )
This method returns the distance matrix in the form of a numpy array. The id_list is a python list of cell ids, specifying the order of the cells in the distance matrix.
When this method is run, it checks if a distance matrix file is found, and if there is not, it calls bcrdist.cell.array.generate_distmatrix(). once the file is generated, it loads in the distances into a numpy array.

## bcrdist.cell.array.name() -> string name
This method returns the name of the array

## bcrdist.cell.array.save([string path]) -> None
This method writes out a save file, that includes all of the cells in the array. This file can later be loaded by calling bcrdist.cell.array.load() or by creating a cell array with the same name as this one. If path is specified, it writes it out to that path, otherwise it writes it out to name + ".bcrsavefile.tsv"

## bcrdist.cell.array.load(string path) -> None
This method loads in data from a bcrdist save file. The cell array this method is called on has to be empty.

## bcrdist.cell.array.tsneplot([string path]) -> None
This method uses the distance matrix to generate a 2d tsne plot and save it to a file. If path is not specified, it defaults to the name + "-tsne-plot.png".

## bcrdist.cell.array.umapplot([string path]) -> None
This method uses the distance matrix to generate a 2d umap plot and save it to a file. If path is not specified, it defaults to the name + "-umap-plot.png".
This method needs the umap library to be installed.

## bcrdist.cell.array.savePCs([n_pcs = 250, string path]) -> None
This method uses the distance matrix to create a high dimentional projection of the cells with kernel pca. The number of pcs is set with n_pcs, and if a path is not specified, the file is saved to the name + "-pcs.csv"
