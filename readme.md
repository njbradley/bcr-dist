# bcr-dist
## the bcr-dist pipeline

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

