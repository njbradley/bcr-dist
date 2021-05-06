# Bcrdist docs

## bcrdist.cellarray
This is the main class that everything in bcrdist is run on. It is an array of bcells, and it can hold either single or paired bcell sequences.

## bcrdist.cellarray([string name])
The constructor for bcrdist.cellarray can take a name parameter. This name will be stored and used in other functions. If a name is not specified, it remains uuninitialized until data is loaded.
When the array is created, it looks to see if a save file exists for the specified name. If it does, that file is loaded.

## bcrdist.cellarray(list data)
This constructor populates the cell array with a list of data. The data is in the same format as data passed to append.


## bcrdist.load10x(string path) -> new cellarray
This method creates a bcrdist.cellarray object and fills it with the cells in the input path. The input path is expected to be a 10x contig annotations file.
This method is shorthand for
```python
arr = bcrdist.cellarray()
arr.load10X(path)
```

## bcrdist.loadBD(string heavypath, [string lightpath]) -> new cellarray
This method, like load10X, creates a cell array and populates it with the bd data files. BD data has one heavy and one light file, so if the data you are using is paired, you have to put in both filenames
This method is shorthand for
```python
arr = bcrdist.cellarray()
arr.loadBD(path)
```

## bcrdist.cellarray.append(tuple data)
This method appends a new cell into the dataset. The tuple of data can have two different formats:
```python
(cell-name, cell-clonotype, heavy-v-gene, heavy-cdr3, light-v-gene, light-cdr3)
(cell-name, cell-clonotype, heavy-cdr1, heavy-cdr2, heavy-cdr3, light-cdr1, light-cdr2, light-cdr3)
```
Right now this method only works for paired data.

## bcrdist.cellarray.load10x(string path) -> None
This method loads bcells from a 10x type data file, and appends the cells onto this cell array.

## bcrdist.cellarray.loadBD(string heavypath, [string lightpath]) -> None
This method loads bcells from one or two bd data files and appends the cells to the cell array. An error is thrown if the loaded cells are not the same type (paired/not paired) as the cells in the array.

## bcrdist.cellarray.generate_dist_matrix() -> None
This method generates a distance matrix with the current cells, and writes it to a file. The file path is the name + ".bcrdistmatrix.tsv"
Unlike bcrdist.cellarray.distmatrix(), this method will generate a matrix even if there is already one saved to file, so it is useful to regenerate the distance matrix. This can also be accomplished by deleting the distance matrix file.

## bcrdist.cellarray.distmatrix() -> ( list id_list, np.array distmatrix )
This method returns the distance matrix in the form of a numpy array. The id_list is a python list of cell ids, specifying the order of the cells in the distance matrix.
When this method is run, it checks if a distance matrix file is found, and if there is not, it calls bcrdist.cellarray.generate_distmatrix(). once the file is generated, it loads in the distances into a numpy array.

## bcrdist.cellarray.name() -> string name
This method returns the name of the array

## bcrdist.cellarray.save([string path]) -> None
This method writes out a save file, that includes all of the cells in the array. This file can later be loaded by calling bcrdist.cellarray.load() or by creating a cell array with the same name as this one. If path is specified, it writes it out to that path, otherwise it writes it out to name + ".bcrsavefile.tsv"

## bcrdist.cellarray.load(string path) -> None
This method loads in data from a bcrdist save file. The cell array this method is called on has to be empty.

## bcrdist.cellarray.tsneplot([string path, colorby]) -> None
## bcrdist.cellarray.umapplot([string path, colorby]) -> None
These methods use the distance matrix to generate a 2d umap/tsne plot and save it to a file. If path is not specified, it defaults to the name + "-umap-plot.png"
or "-tsne-plot.png". This method needs the umap library to be installed.
The colorby parameter is used to color points. A function that returns a valid matplotlib color can be passed in. The function is passed a python list
in the format that is returned from bcrdist.cellarray.tolist(), for every cell.

## bcrdist.cellarray.savePCs([n_pcs = 250, string path]) -> None
This method uses the distance matrix to create a high dimentional projection of the cells with kernel pca. The number of pcs is set with n_pcs, and if a path is not specified, the file is saved to the name + "-pcs.csv"

## bcrdist.cellarray.tolist() -> list
This method returns a python representation of this cellarray. It returns a list of tuples, each tuple represents one cell. The tuples are
in these two different formats, depending on whether it is a paired or unpaired array
```python
(cell-barcode, clonotype, cdr1, cdr2, cdr3)
(cell-barcode, clonotype, heavy-cdr1, heavy-cdr2, heavy-cdr3,light-cdr1, light-cdr2, light-cdr3)
```


