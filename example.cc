#include "dataset_loader.hpp"
#include <chrono>

/// Example loading and showing some values from an NLOS dataset.
int main(int argc, char* argv[]) {
    nlos::NLOSDataset data = nlos::DatasetLoader::read_NLOS_dataset(argv[1], {3}, false);
    std::cout << "Bin resolution: " << data.deltat[0] << '\n' <<
                 "Bins: " << data.bins[0] << '\n' << 
                 "Cam grid dimensions: " << data.camera_grid_dimensions[0] << '\n' << 
                 "Data size: ";
    for (auto const i :  data.data.shape()) {
        std::cout << i << ' ';
    }
    std::cout << '\n';
    std:: cout << "Is confocal?: " << (data.capture == nlos::CaptureStrategy::Confocal ? "true" : "false") << '\n' << 
    "Camera positions: " << data.camera_grid_positions[0];
}
