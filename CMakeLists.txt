project(dsreader)

cmake_minimum_required(VERSION 3.11)

add_executable(loader_example example.cc)
target_compile_features(loader_example PRIVATE cxx_std_14)
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
target_include_directories(loader_example PRIVATE ${HDF5_INCLUDE_DIRS} xtl/include xtensor/include)
target_link_libraries(loader_example ${HDF5_LIBRARIES})
