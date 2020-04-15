#ifndef NLOS_ENUMS_HPP
#define NLOS_ENUMS_HPP

namespace nlos
{
    
/**
 * @brief Specifies a compute paradigm
 */
enum Compute
{
    CPU, 
    GPU,
    ComputeNone,
};

/**
 * @brief Specifies the strategy to access voxel volume
 */
enum VolumeAccess
{
    Naive,
    Octree,
    VolumeAccessNone,
};

/**
 * @brief Specifies the way a multidimensional array is stored.
 * Specifies the way a multidimensional array is stored. 
 * Either RowMajor or ColumnMajor order.
 * 
 * In the future, other patterns may make sense, like:
 * https://en.wikipedia.org/wiki/Z-order_curve
 */
enum DataOrder
{
    ColumnMajor,
    RowMajor,
    DataOrderNone,
};

/**
 * @brief Specifies the capture strategy for an NLOS scene
 */
enum CaptureStrategy
{
    Exhaustive,
    Confocal,
    CaptureStrategyNone,
};

}; // namespace nlos

#endif // NLOS_ENUMS_HPP