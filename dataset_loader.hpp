#ifndef DATASET_LOADER_HPP
#define DATASET_LOADER_HPP

#include <vector>
#include <H5Cpp.h>
#include <xtensor/xarray.hpp>
#include <xtensor/xfunction.hpp>
#include <xtensor/xadapt.hpp>
#include "nlos_dataset.hpp"

namespace nlos
{

/**
 * A loader for HDF5 NLOS datasets
 * 
 * Note: in HDF5 files each entry or field is called a dataset.
 * This clashes with our notion of dataset here, which refers to the data and
 * metadata from an NLOS scene.
 */ 
class DatasetLoader {
private:
    /// Constant field name definitions.
    static constexpr const char* DS_CAM_GRID_POSITIONS = "cameraGridPositions";
    static constexpr const char* DS_CAM_GRID_NORMALS = "cameraGridNormals";
    static constexpr const char* DS_CAM_POSITION = "cameraPosition";
    static constexpr const char* DS_CAM_GRID_POINTS = "cameraGridPoints";
    static constexpr const char* DS_CAM_GRID_SIZE = "cameraGridSize";
    static constexpr const char* DS_LASER_GRID_POSITIONS = "laserGridPositions";
    static constexpr const char* DS_LASER_GRID_NORMALS = "laserGridNormals";
    static constexpr const char* DS_LASER_POSITION = "laserPosition";
    static constexpr const char* DS_LASER_GRID_POINTS = "laserGridPoints";
    static constexpr const char* DS_LASER_GRID_SIZE = "laserGridSize";
    static constexpr const char* DS_DATA = "data";
    static constexpr const char* DS_DELTA_T = "deltaT";
    static constexpr const char* DS_T0 = "t0";
    static constexpr const char* DS_T = "t";
    static constexpr const char* DS_HIDDEN_VOLUME_POSITION = "hiddenVolumePosition";
    static constexpr const char* DS_HIDDEN_VOLUME_ROTATION = "hiddenVolumeRotation";
    static constexpr const char* DS_HIDDEN_VOLUME_SIZE = "hiddenVolumeSize";
    static constexpr const char* DS_IS_CONFOCAL = "isConfocal";

    static constexpr const char* ATT_THIRD_BOUNCE = "third_bounce";

    /**
     * Loads a full field from the HDF5 file
     * 
     * @param dataset the field to fully load from the file
     * @tparam T field type to load. Admits single precision floats and 32bit integers
     * @return 
     */ 
    template <typename T>
    static xt::xarray<T> load_field_array(const H5::DataSet &dataset) 
    {
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
        xt::xarray<T> retval = xt::adapt(buff, num_elements, xt::acquire_ownership(), dimensions);
        return retval;
    }

    /**
     * Loads the specified field as the data field with the specified parameters.
     * 
     * @warning It is assumed the data is stored as single precision floats.
     * 
     * @param dataset data field in the HDF5 file. This is not enforced however, so it must be used with care.
     * @param bounces the light bounces that will be loaded. 
     * 
     */ 
    template <typename T>
    static xt::xarray<T> load_transient_data_dataset(const H5::DataSet &dataset, 
                                                     const std::vector<uint32_t>& bounces,
                                                     bool sum_bounces=false,
                                                     DataOrder data_order=ColumnMajor) 
    {
        assert(bounces.size() > 0);

        // Assumes the third bounce matches the second element of the bounce 
        int third_bounce = 1;
        if (dataset.attrExists(ATT_THIRD_BOUNCE))
        {
            H5::Attribute tbattr = dataset.openAttribute(ATT_THIRD_BOUNCE);
            tbattr.read(H5::PredType::NATIVE_INT, &third_bounce);
        }

        H5::DataSpace dataspace = dataset.getSpace();
        int rank = dataspace.getSimpleExtentNdims();
        std::vector<hsize_t> dimensions(rank);
        dataspace.getSimpleExtentDims(dimensions.data(), nullptr);
        size_t bounce_axis = data_order == RowMajor ? rank - 3 : 2;
        hsize_t num_elements = 1;
        for (int i = 0; i < rank; i++) {
            // Account only for the chosen bounces
            if (i == bounce_axis) num_elements *= bounces.size();
            else num_elements *= dimensions[i];
        }
        T *buff = (T*) new uint8_t[num_elements*sizeof(T)];
        auto ptype = H5::PredType::NATIVE_FLOAT;
        {
            std::vector<hsize_t> offset(rank, 0);
            std::vector<hsize_t> count(rank, 0);
            for (int i = 0; i < rank; i++)
                count[i] = dimensions[i];
            count[bounce_axis] = 1;
            dataspace.selectNone();
            for (uint32_t b = 0; b < bounces.size(); b++)
            {
                // Bounces in the dataset start at 2, so the 3rd bounce is the 1st element
                offset[bounce_axis] = bounces[b]- (3 - third_bounce);
                dataspace.selectHyperslab(H5S_SELECT_OR, count.data(), offset.data());
            }
        }
        dimensions[bounce_axis] = bounces.size();
        
        H5::DataSpace mspace = H5::DataSpace(rank, dimensions.data());

        dataset.read(buff, ptype, mspace, dataspace);
        xt::xarray<T> retval = xt::adapt(buff, num_elements, xt::no_ownership(), dimensions);
        if (sum_bounces && bounces.size() > 1)
        {
            dimensions[bounce_axis] = 1;
            retval = xt::sum(retval, {bounce_axis});
            retval.reshape(dimensions);
        }
        return retval;
    }

public:
    /**
     * @brief 
     * 
     * @param dataset 
     * @param file_path 
     * @param bounces 
     * @param sum_bounces 
     * @param data_order 
     */
    static void read_NLOS_dataset(NLOSDataset& dataset, 
                                  std::string file_path, 
                                  const std::vector<uint32_t>& bounces, 
                                  bool sum_bounces=false,
                                  DataOrder data_order=DataOrder::DataOrderNone)
    {
        dataset = read_NLOS_dataset(file_path, bounces, sum_bounces, data_order);
    }

    /**
     * @brief 
     * 
     * @param file_path 
     * @param bounces 
     * @param sum_bounces 
     * @param data_order 
     * @return NLOSDataset 
     */
    static NLOSDataset read_NLOS_dataset(std::string file_path, 
                                         const std::vector<uint32_t>& bounces,
                                         bool sum_bounces=false,
                                         DataOrder data_order=DataOrder::DataOrderNone)
    {
        NLOSDataset dataset;
        H5::H5File file(file_path, H5F_ACC_RDONLY);
        if (file.attrExists("data order"))
        {
            H5::Attribute att(file.openAttribute("data order"));
            H5::StrType stype = att.getStrType();
            std::string data_order_string;
            att.read(stype, data_order_string);
            if (data_order_string.compare("row-major") == 0)
            {
                dataset.data_order = RowMajor;
            }
            else if (data_order_string.compare("column-major") == 0)
            {
                dataset.data_order = ColumnMajor;
            }
        }
        if (file.attrExists("engine"))
        {
            H5::Attribute att(file.openAttribute("engine"));
            H5::StrType stype = att.getStrType();
            att.read(stype, dataset.engine);
        }

        dataset.data = load_transient_data_dataset<float>(file.openDataSet(DS_DATA), bounces, sum_bounces, dataset.data_order);
		dataset.camera_grid_positions = load_field_array<float>(file.openDataSet(DS_CAM_GRID_POSITIONS));
        dataset.camera_grid_normals = load_field_array<float>(file.openDataSet(DS_CAM_GRID_NORMALS));
        dataset.camera_position = load_field_array<float>(file.openDataSet(DS_CAM_POSITION));
        dataset.camera_grid_dimensions = load_field_array<float>(file.openDataSet(DS_CAM_GRID_SIZE)); 
        dataset.camera_grid_points = load_field_array<float>(file.openDataSet(DS_CAM_GRID_POINTS)); 
        dataset.laser_grid_positions = load_field_array<float>(file.openDataSet(DS_LASER_GRID_POSITIONS)); 
        dataset.laser_grid_normals = load_field_array<float>(file.openDataSet(DS_LASER_GRID_NORMALS)); 
        dataset.laser_position = load_field_array<float>(file.openDataSet(DS_LASER_POSITION));
        dataset.laser_grid_dimensions = load_field_array<float>(file.openDataSet(DS_LASER_GRID_SIZE)); 
        dataset.laser_grid_points = load_field_array<float>(file.openDataSet(DS_LASER_GRID_POINTS)); 
        dataset.hidden_volume_position = load_field_array<float>(file.openDataSet(DS_HIDDEN_VOLUME_POSITION)); 
        dataset.hidden_volume_rotation = load_field_array<float>(file.openDataSet(DS_HIDDEN_VOLUME_ROTATION)); 
        dataset.hidden_volume_size = load_field_array<float>(file.openDataSet(DS_HIDDEN_VOLUME_SIZE));
        dataset.t0 = load_field_array<float>(file.openDataSet(DS_T0));
        dataset.bins = load_field_array<int>(file.openDataSet(DS_T));
        dataset.deltat = load_field_array<float>(file.openDataSet(DS_DELTA_T));
        dataset.capture = load_field_array<int>(file.openDataSet(DS_IS_CONFOCAL))[0] ? CaptureStrategy::Confocal : CaptureStrategy::Exhaustive;
        
        if (data_order != DataOrderNone && dataset.data_order != data_order)
        {
            std::cout << "Transposing data\n";
            dataset.data = xt::transpose(dataset.data);
            dataset.camera_grid_positions = xt::transpose(dataset.camera_grid_positions);
            dataset.camera_grid_normals = xt::transpose(dataset.camera_grid_normals);
            dataset.laser_grid_positions = xt::transpose(dataset.laser_grid_positions);
            dataset.laser_grid_normals = xt::transpose(dataset.laser_grid_normals);
            dataset.data_order = data_order;
        }
        return std::move(dataset);
    }
};

}; // namespace nlos

#endif // DATASET_LOADER_HPP