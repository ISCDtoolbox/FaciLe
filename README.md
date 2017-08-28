# FaciLe
Reconstruction faciale

## Installation
Clone the repository with the recursive option activated:
```shell
git clone --recursive https://github.com/ISCDtoolbox/FaciLe.git
```
This will fetch all useful dependencies for the whole project:
* [Super4PCS](https://github.com/nmellado/Super4PCS.git) for global registration
* [ISCD Commons Library](https://github.com/ISCDtoolbox/Commons.git)
* [ISCD LinearElasticity Library](https://github.com/ISCDtoolbox/LinearElasticity.git)
* [MshDist](https://github.com/ISCDtoolbox/Mshdist.git) for signed distances
* [tetgen](https://github.com/ufz/tetgen.git) for tetrahedral meshes creation
* [mmgTools](https://github.com/MmgTools/mmg.git) for surface and volume remeshing

To install all the previous utilities, as well as the Facial reconstruction pipeline, run in a terminal:
```shell
mkdir build
cd build
cmake ..
make
```
