import sys

def help():
    print ("bcrdist command line tools: version 0.1")
    print ("This command can run bcrdist on a 10x or BD dataset")
    print ("Usage: bcrdist [options] input-files")
    print ("Options:")
    print ("  --type: the type of dataset given, either 10x or BD")
    print ("  --pcs: the path to save the PC table")
    print ("  --tsne: the path to save the tsne plot")
    print ("  --umap: the path to save the umap plot")

def main():
    args = []
    kwargs = {}
    
    name = ""
    for arg in sys.argv[1:]:
        if (arg[0] == '-'):
            if (name != ""):
                raise ValueError("No value specified for arg '" + name + "'")
            name = arg.replace('-','')
        elif (name != ""):
            kwargs[name] = arg
            name = ""
        else:
            args.append(arg)
    
    if (name == 'h' or name == 'help'):
        help()
    
    if (len(args) < 1):
        print ("ERR: no input files specified")
        help()
        return
    
    import bcrdist
    
    if (kwargs['type'] == "10x"):
        arr = bcrdist.load10x(args[0])
    elif (kwargs['type'] == "BD"):
        arr = bcrdist.loadBD(*args)
    else:
        print ("ERR: invalid input type")
        return
    
    if ('pcs' in kwargs):
        arr.savePCs(kwargs['pcs'])
    else:
        arr.savePCs()
    
    if ('tsne' in kwargs):
        arr.tsneplot(kwargs['tsne'])
    else:
        arr.tsneplot()
    
    if ('umap' in kwargs):
        arr.umapplot(kwargs['umap'])
    else:
        arr.umapplot()
    
    print ("Completed without errors!")

main()
