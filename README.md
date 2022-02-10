![MCvsBC](https://user-images.githubusercontent.com/38473042/152366380-3695f20e-606b-429a-b103-4c125ae17f31.png)


# MC3D – An algorithm for Structured Volume Decomposition

`MC3D` is an implementation of [The 3D Motorcycle Complex for Structured Volume Decomposition \[Brückler et al. 2021\]](https://arxiv.org/abs/2112.05793) (accepted to Eurographics 2022) distributed under GPLv3.

If you make use of `MC3D` in your scientific work, please cite our paper. For your convenience,
you can use the following bibtex snippet:

    @article{DBLP:journals/corr/abs-2112-05793,
        author     = {Hendrik Br{\"{u}}ckler and
                     Ojaswi Gupta and
                     Manish Mandad and
                     Marcel Campen},
        title      = {The 3D Motorcycle Complex for Structured Volume Decomposition},
        journal    = {CoRR},
        volume     = {abs/2112.05793},
        year       = {2021},
        url        = {https://arxiv.org/abs/2112.05793},
        eprinttype = {arXiv},
        eprint     = {2112.05793},
    }

***

![sequence](https://user-images.githubusercontent.com/38473042/152368807-5c51f045-a127-4052-ab04-d3fd60f612fa.png)


## What is the 3D Motorcycle Complex?

The 3D Motorcycle Complex partitions a tetrahedral mesh, equipped with a suitable seamless map, into blocks that are axis-aligned rectangular cuboids under the parametrization. The resulting partition is non-conforming (i.e. one where "T"-junctions exist) and therefore can be made much coarser than similar conforming partitions like the so called base complex (BC).

For a visual impression of what that means, check out our demonstration video:

https://user-images.githubusercontent.com/38473042/152365682-b0052ec2-df41-447a-bec4-c30865ea71b4.mp4


***

### Dependencies
- GMP (NOT included, must be installed on your system)
- [TrulySeamless3D](https://github.com/HendrikBrueckler/TrulySeamless3D) (Included as submodule, together with subdependency libHexEx)
- OpenVolumeMesh (Included as submodule. A PATCHED version is required to fully support selfadjacent blocks! If you use CMake to build, it will automatically apply the patch in ```extern/patches/```, otherwise you have to apply the patch yourself.)
- glog (Included as submodule)
- googletest (Included as submodule)
- CLI11 (Included as submodule)

### Building
In root directory

    mkdir build
    cd build
    cmake [-DMC3D_BUILD_CLI=Off] [-DMC3D_BUILD_TESTS=Off] [-DMC3D_ENABLE_LOGGING=Off] ..
    make

### Usage
An example command-line application is included that reads a tetrahedral mesh including a seamless parametrization from a file in .hexex-format, as used and documented in [libHexEx](https://www.graphics.rwth-aachen.de/software/libHexEx/).
It outputs a file in the same format, containing a refined version of the input mesh, followed by the number and the list of wall triangles (three vertex indices, and a double indicating the wall triangle's distance from the brush fire's origin per line) at the end of the file.

After building the CLI app can be found in ```build/Build/bin/cli``` .
For full information on its usage, execute

    mc3d_cli --help

Example input can be found in folder ```tests/resources```.

### API
For details on the API of the library, check the headers in ```include```, they are thoroughly documented. Apart from that, ```cli/main.cpp``` demonstrates usage of the entire MC3D-Pipeline for both simple and advanced usage.
