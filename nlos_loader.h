#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <H5Cpp.h>
#include <chrono>

#ifdef USE_XTENSOR
#include <xtensor/xarray.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xindex_view.hpp>
#include <xtensor/xcontainer.hpp>
#include <xtensor/xfunction.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xmath.hpp>
#define array_type xt::xarray
#else
template <typename T>
struct simple_tensor {
    std::vector<uint> shape;
    T *buff;
    uint total_elements;
    
    // Prints the nd-array
    void print() {
        uint depth = 0;
        using std::cout;
        uint dims = this->shape.size();
        std::vector<uint> dim_ctrs(dims);
        for (uint i = 0; i < this->total_elements; i++) {
            for (int dc = dims-1; dc >= 0; dc--) {
                if (dim_ctrs[dc] == this->shape[dc]) {
                    cout << "}";
                    dim_ctrs[dc] = 0;
                    dim_ctrs[dc-1]++;
                    depth--;
                } else {
                    if (dc != dims-1) cout << ",\n";
                    break;
                }
            }
            if (depth < dims) {
                for (int c = 0; c < depth; c++) cout << ' ';
            }
            while (depth < dims) {
                depth++;
                cout << "{";
            }
            cout << buff[i];
            dim_ctrs[dims-1]++;
            if (dim_ctrs[dims-1] != this->shape[dims-1]) cout << ", ";
        }
        for (int dc = dims-1; dc >= 0; dc--) {
            cout << '}';
        }
        cout << '\n';
    }
};
#define array_type simple_tensor
#endif

class nlos_loader {
    private:
    /// Constant field name definitions.
    const std::string DS_CAM_GRID_POSITIONS = "cameraGridPositions";
    const std::string DS_CAM_GRID_NORMALS = "cameraGridNormals";
    const std::string DS_CAM_POSITION = "cameraPosition";
    const std::string DS_CAM_GRID_POINTS = "cameraGridPoints";
    const std::string DS_CAM_GRID_SIZE = "cameraGridSize";
    const std::string DS_LASER_GRID_POSITIONS = "laserGridPositions";
    const std::string DS_LASER_GRID_NORMALS = "laserGridNormals";
    const std::string DS_LASER_POSITION = "laserPosition";
    const std::string DS_LASER_GRID_POINTS = "laserGridPoints";
    const std::string DS_LASER_GRID_SIZE = "laserGridSize";
    const std::string DS_DATA = "data";
    const std::string DS_DELTA_T = "deltaT";
    const std::string DS_T0 = "t0";
    const std::string DS_T = "t";
    const std::string DS_HIDDEN_VOLUME_POSITION = "hiddenVolumePosition";
    const std::string DS_HIDDEN_VOLUME_ROTATION = "hiddenVolumeRotation";
    const std::string DS_HIDDEN_VOLUME_RADIUS = "hiddenVolumeSize";
    const std::string DS_IS_CONFOCAL = "isConfocal";

    template <typename T>
    array_type<T> load_field_array(const H5::DataSet &dataset) {
        H5T_class_t type_class = dataset.getTypeClass(); // Check the data type
        H5::DataSpace dataspace = dataset.getSpace();
        int rank = dataspace.getSimpleExtentNdims();
        std::vector<hsize_t> dimensions(rank);
        dataspace.getSimpleExtentDims(dimensions.data(), nullptr);
        hsize_t num_elements = 1;
        for (int i = 0; i < rank; i++) {
            num_elements *= dimensions[i];
        }
        T *buff = (T*) new uint8_t[num_elements*sizeof(T)];
        auto ptype = H5::PredType::NATIVE_FLOAT;
        switch(type_class) {
            case H5T_INTEGER:
                ptype = H5::PredType::NATIVE_INT32;
                break;
            case H5T_FLOAT:
            default:
                ptype = H5::PredType::NATIVE_FLOAT;
        }
        dataset.read(buff, ptype);
        #ifdef USE_XTENSOR
        array_type<T> retval = xt::adapt(buff, num_elements, xt::acquire_ownership(), dimensions);
        #else
        array_type<T> retval;
        retval.buff = buff;
        retval.total_elements = num_elements;
        retval.shape.resize(rank);
        for (int i = 0; i < rank; i++) {
            retval.shape[i] = dimensions[i];
        }
        #endif
        return retval;
    }

    public:
    nlos_loader(std::string file_path) {
        H5::H5File file(file_path, H5F_ACC_RDONLY);
        data = load_field_array<float>(file.openDataSet(DS_DATA));
        camera_point_positions = load_field_array<float>(file.openDataSet(DS_CAM_GRID_POSITIONS));
        camera_grid_normals = load_field_array<float>(file.openDataSet(DS_CAM_GRID_NORMALS));
        camera_position = load_field_array<float>(file.openDataSet(DS_CAM_POSITION));
        camera_grid_dimensions = load_field_array<float>(file.openDataSet(DS_CAM_GRID_SIZE)); 
        camera_grid_points = load_field_array<float>(file.openDataSet(DS_CAM_GRID_POINTS)); 
        laser_grid_positions = load_field_array<float>(file.openDataSet(DS_LASER_GRID_POSITIONS)); 
        laser_grid_normals = load_field_array<float>(file.openDataSet(DS_LASER_GRID_NORMALS)); 
        laser_position = load_field_array<float>(file.openDataSet(DS_LASER_POSITION));
        laser_grid_dimensions = load_field_array<float>(file.openDataSet(DS_LASER_GRID_SIZE)); 
        laser_grid_points = load_field_array<float>(file.openDataSet(DS_LASER_GRID_POINTS)); 
        hidden_volume_position = load_field_array<float>(file.openDataSet(DS_HIDDEN_VOLUME_POSITION)); 
        hidden_volume_rotation = load_field_array<float>(file.openDataSet(DS_HIDDEN_VOLUME_ROTATION)); 
        hidden_volume_size = load_field_array<float>(file.openDataSet(DS_HIDDEN_VOLUME_RADIUS));
        t0 = load_field_array<int>(file.openDataSet(DS_T0));
        bins = load_field_array<int>(file.openDataSet(DS_T));
        deltat = load_field_array<float>(file.openDataSet(DS_DELTA_T));
        is_confocal = load_field_array<int>(file.openDataSet(DS_IS_CONFOCAL));
    }
    
    // Spad capture volume
    array_type<float> data;

    // Camera/Spad
    array_type<float> camera_point_positions; // Position of every recorded point of the grid
    array_type<float> camera_grid_normals; // Normal of every recorded point of the grid
    array_type<float> camera_position; // Camera origin
    array_type<float> camera_grid_dimensions; // Dimensions of the camera point grid
    array_type<float> camera_grid_points; // Number of capture points in the grid in X and Y.

    // Laser
    array_type<float> laser_grid_positions; // Position of every traced point of the grid
    array_type<float> laser_grid_normals; // Normal of every traced point of the grid
    array_type<float> laser_position; // Laser origin
    array_type<float> laser_grid_dimensions; // Dimensions of the laser point grid
    array_type<float> laser_grid_points; // Number of laser points in the grid in X and Y.

    // Scene info
    array_type<float> hidden_volume_position; // Center of the hidden geometry
    array_type<float> hidden_volume_rotation; // Hidden geometry rotation with respect 
    // to the ground truth
    /// These next are arrays for consistency, but they should be single values ///
    array_type<float> hidden_volume_size; // Half-length of a cube containing the hidden
    // geometry
    array_type<int> t; // Time resolution
    array_type<int> t0; // Time at which the captures start (first data column)
    array_type<int> bins;  // Number of time instants recorded (number of columns in the data)
    array_type<float> deltat;  // Per pixel aperture duration (time resolution)
    array_type<int> is_confocal; // Boolean value. 1 if the dataset is confocal, 0 if all combinations 
    // of laser points and spad points were captured/rendered
};