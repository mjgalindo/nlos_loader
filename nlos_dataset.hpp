#ifndef NLOS_DATASET_HPP
#define NLOS_DATASET_HPP

#include <string>
#include <xtensor/xarray.hpp>

#include "nlos_enums.hpp"

namespace nlos
{
/**
 * @brief Represents an NLOS dataset and metadata
 */
class NLOSDataset
{
public:
    // Spad capture volume
    xt::xarray<float> data;

    // Camera/Spad
    xt::xarray<float> camera_grid_positions; // Position of every recorded point of the grid
    xt::xarray<float> camera_grid_normals; // Normal of every recorded point of the grid
    xt::xarray<float> camera_position; // Camera origin
    xt::xarray<float> camera_grid_dimensions; // Dimensions of the camera point grid
    xt::xarray<float> camera_grid_points; // Number of capture points in the grid in X and Y.

    // Laser
    xt::xarray<float> laser_grid_positions; // Position of every traced point of the grid
    xt::xarray<float> laser_grid_normals; // Normal of every traced point of the grid
    xt::xarray<float> laser_position; // Laser origin
    xt::xarray<float> laser_grid_dimensions; // Dimensions of the laser point grid
    xt::xarray<float> laser_grid_points; // Number of laser points in the grid in X and Y.

    // Scene info
    xt::xarray<float> hidden_volume_position; // Center of the hidden geometry
    xt::xarray<float> hidden_volume_rotation; // Hidden geometry rotation with respect 
    // to the ground truth
    /// These next are arrays for consistency, but they should be single values ///
    xt::xarray<float> hidden_volume_size; // Dimensions of prism containing the hidden geometry
    xt::xarray<int> t; // Time resolution
    xt::xarray<float> t0; // Time at which the captures start 
    xt::xarray<int> bins;  // Number of time instants recorded (number of columns in the data)
    xt::xarray<float> deltat;  // Per pixel aperture duration (time resolution)
    CaptureStrategy capture;

    DataOrder data_order = DataOrder::ColumnMajor;
    std::string engine = "default";
};
}
#endif // NLOS_DATASET_HPP
