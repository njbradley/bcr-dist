from . import cellarray

def load10x(path):
    arr = cellarray()
    arr.load10x(path)
    return arr

def loadBD(*args):
    arr = cellarray()
    arr.loadBD(*args)
    return arr

def loaddekosky(path):
    arr = cellarray()
    arr.loaddekosky(path)
    return arr

def load(path):
    return cellarray(path)
