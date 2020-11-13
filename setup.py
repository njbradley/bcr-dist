import setuptools

cext = setuptools.Extension(
    "bcrdist.cbcrdist",
    sources = [
        'bcrdist/cmodule/bcrdistmodule.cc',
        'bcrdist/cmodule/scripts/table.cc',
        'bcrdist/cmodule/scripts/cell.cc',
        'bcrdist/cmodule/scripts/data.cc',
        'bcrdist/cmodule/scripts/fileio.cc',
        # 'bcrdist/cmodule/scripts/distances.cc',
    ],
    include_dirs = [
        '/app/software/Anaconda3/2020.02/lib/python3.7/site-packages/numpy/core/include',
        'bcrdist/cmodule/scripts/packages'
    ]
)

setuptools.setup(
    name = "bcrdist",
    version = "0.1",
    description = "a scientific library that computes the relative distances of bcr sequences",
    packages = ["bcrdist"],
    ext_modules = [cext],
    package_data = {
        "bcrdist": ["data/*.*"]
    },
    # entry_points = {
    #     "console_scripts": [
    #         'bcrdist:
)
