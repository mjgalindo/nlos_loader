# C++ NLOS Dataset Loader

This is a simple C++ header-only library for loading synthetic NLOS datasets.
It works with both `hdf5` and `mat` files so long as the following fields can be found in their root.

| Field | Explanation |
|-------|-------------|
| cameraPosition|Position of the SPAD|
| cameraGridNormals | List of capture point normals. Points start at the top left of the grid going to the right and bottom.|
| cameraGridPositions | List of capture point positions. Points start at the top left of the grid going to the right and bottom.|
| cameraGridSize | Width and height of the capture grid.|
| cameraGridPoints |  Number of capture points in the grid in X and Y.|
| laserPosition | Position of the laser that generates the virtual point lights in the grid.|
| laserGridPositions | Positions of each virtual point light. Virtual point lights in the grid start at the top left of the grid going to the right and bottom. Since the grid is in a plane they are the same in most cases.|
| laserGridNormals | List of virtual point light normals in the grid. |
| laserGridSize | Width and height of the virtual point light grid.|
| laserGridPoints |  Number of laser points in the grid in X and Y.|
| data | Monochrome transient image container. Its size is N_LASER_SPOTS x N_SPAD_SPOTS x BOUNCES x TIME_RES. In a non-confocal setting, the size may be N_LASER_SPOTS_X x N_LASER_SPOTS_Y x N_SPAD_SPOTS_X x N_SPAD_SPOTS_X x BOUNCES x TIME_RES. Light bounces are stored separately in order to test how MPI affects algorithms or add filters to later bounces to improve reconstruction quality. |
| deltaT | Time resolution for each pixel on the transient image in distance units. For instance, deltaT=0.001 means each pixels captures light for the time it takes it to travel 0.001m in the scene with c=1m/s.|
| hiddenVolumePosition | Center position of the volume containing the hidden geometry. Useful for backprojection algorithms that only work on a volume.|
| hiddenVolumeRotation | Rotation of a volume containing the hidden geometry with respect to the original geometry. Useful when comparing reconstructions.|
| hiddenVolumeSize | Half-length of the box that tightly bounds the hidden geometry.|
| isConfocal | True if the data was rendered confocally, that is, the spad and laser grids are the same size, have the same number of points and only the positions where both match were rendered/captured.|
| t | Time resolution.|
| t0 | First instant captured in each transient row.|

By default, all data is loaded into structs containing the raw carray, its size and its shape (in case it's an n-dimensional tensor as is the case for the captures).
To work with the data more easily, you can use `xtensor` as described in [dependencies](##Dependencies).

A boolean value `is_row_major` indicates whether the previous data shapes are changed for a row-major scheme with the fastest-changing dimensions in the later axes.

## Dependencies

To use this header you'll need the HDF5 library which you can find here [https://www.hdfgroup.org](https://www.hdfgroup.org/) or in your repositories.
Optionally, you can use [xtensor](http://quantstack.net/xtensor) to load the data directly into `xtensor` objects. In this case, `USE_XTENSOR` must be defined before including `nlos_reader.h` and the appropriate headers for `xtensor` will be included and used.

## Usage

To use it in a C++ project just include `nlos_loader.h` and link the HDF5 library. To simplify builds you can use the cmake commands
```{cmake}
find_package (HDF5 REQUIRED COMPONENTS CXX)
target_include_directories(<your_target> ${HDF5_INCLUDE_DIRS})
target_link_libraries(<your_target> ${HDF5_LIBRARIES})
``` 
as in the CMakeLists.txt file included for building the example.

To load a file, just create an `NLOSData` instance from its filename string and the bounces you want to load (starting from the 3rd and up to the 7th).
For instance, `NLOSData("example.hdf5", {3,4,5})` will load all the metadata from `example.hdf5` including the transient data for the 3rd, 4th and 5th bounces.
The constructor will load all the values from the file, throwing an exception if the file can't be read.