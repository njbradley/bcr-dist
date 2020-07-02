import distutils.core as dcore
import os

module1 = dcore.Extension('bcrdist',
                    sources = ['bcrdistmodule.cc', 'scripts/table.cc', 'scripts/cell.cc', 'scripts/data.cc', 'scripts/fileio.cc'])

dcore.setup (name = 'bcr-dist',
       version = '1.0',
       description = 'This is a demo package',
       ext_modules = [module1])

#os.system('cp build/lib.linux-x86_64-3.7/bcrdist.cpython-37m-x86_64-linux-gnu.so ./')
