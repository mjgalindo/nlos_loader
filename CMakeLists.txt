project(dsreader)

cmake_minimum_required(VERSION 3.11)

set (CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE Release)
add_executable(loader_example example.cc)

option(USE_XTENSOR OFF)

if (USE_XTENSOR)
    add_definitions(-DUSE_XTENSOR)
    if (USE_XSIMD)
        add_definitions(-DXTENSOR_ENABLE_XSIMD)
    endif(USE_XSIMD)
    message("Using xtensor arrays")
else()
    message("Using raw c-arrays")
endif(USE_XTENSOR)

find_package (HDF5 REQUIRED COMPONENTS CXX)
include_directories(${HDF5_INCLUDE_DIRS})
target_link_libraries(loader_example ${HDF5_LIBRARIES})
