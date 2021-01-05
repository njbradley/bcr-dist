# bcr-dist
## the bcr-dist pipeline

It is inspired by tcr-dist, (https://github.com/phbradley/tcr-dist).

This is bcrdist, a python library that computes the relative distances between bcr chains. These distances can then be used to construct a high dimentional mapping of the cells through kernel pca, or create tsne/umap plots to analyze the cells.

# Installation

## Installation requirements

In order to install bcrdist, you need to have a python 3 environment, and a c++ compiler. This should be automatically installed for most linux and mac computers. If you are on windows, you may need to install visual c++.

After cloning or downloading the repository, there are two different ways to install bcrdist:

## Local installation

A local installation creates a library folder that can only be used by python scripts in the same directory. The advantages of this method is that the module is not added to your python environment.

Run the command: ```python setup.py build```

This will create a new build directory, with a lib folder inside. Inside that lib folder is the bcrdist package. You can move this to wherever you need to import bcrdist.

## Global installation

Installing bcrdist globally means you will be able to import it from any python script, regardless of the location. You will need to have root access to your computer

Run the command: ```python setup.py install```

This will install bcrdist globally. Test to make sure it worked by trying to import bcrdist from the python interpreter.

# Usage

Bcrdist is a python library, meaning that you can import it into other python scripts. If you want to get up and running quickly and don't want to customize the process, there are many ready to use scripts in the tests/ folder. To run these, call the python program with the input files as parameters. For example:

```python 10x_test.py filtered_contig_annotations.csv```

This script will load in the specified cell dataset and produce tsne and umap plots, as well as a matrix of principal components for each cell. This example uses the 10x format of input files, but there is also a bd_test.py for bd data.

## Custom scripts

If you want to customize the processing pipeline, you will need to write your own scripts using the bcrdist library. The api documentiation is in the docs.md file, which gives an overview of all of the classes / methods of bcrdist. However, here is a general overview of how a script would work

```python
import bcrdist

cellarray = bcrdist.cell.load10X("filtered_contig_annotations.csv")
cellarray.savePCS()
cellarray.umapplot()
cellarray.tsneplot()
```

In this example, the filtered_contig_annotations.csv file is loaded in, and the PCS are saved, as well as a tsne and umap plot. A different way to do this is to split up the data processing and the plotting and saving.

calculations.py:
```python
import bcrdist

cellarray = bcrdist.cell.load10X("filtered_contig_annotations.csv")
cellarray.distmatrix()
cellarray.save()
```

plots.py:
```python
import bcrdist

cellarray = bcrdist.cell.array("filtered_contig_annotations.csv")
cellarray.savePCS()
cellarray.umapplot()
cellarray.tsneplot()
```

In this example, the calculations are done in calculations.py, and the plotting is done in plots.py. This approach has the benefit that all of the time consuming processes are done first and saved, so if you are tweaking the plot parameters and rerunning plot.py many times, it will not have to rerun all of the analysis.

Hopefully this has given you an idea of what the general structure of a script will be. The docs.md has specific documentation for methods.
