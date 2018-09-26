#include "nlos_loader.h"
#include <chrono>

/// Example loading and showing some values from an NLOS dataset.
int main(int argc, char* argv[]) {
    nlos_loader data(argv[1]);
    #ifdef USE_XTENSOR
    std::cout << "Bin resolution: " << data.deltat << '\n' <<
                 "Bins: " << data.bins << '\n' << 
                 "Cam grid dimensions: " << data.camera_grid_dimensions << '\n' << 
                 "Data size: ";
    for (auto const i :  data.data.shape()) {
        std::cout << i << ' ';
    }
    std::cout << '\n';
    std:: cout << "Is confocal?: " << data.is_confocal << '\n' << 
    "Camera positions: " << data.camera_point_positions;
#else 
    std::cout << "Bin resolution: "; data.deltat.print();
    std::cout << "Bins: "; data.bins.print();
    std::cout << "Cam grid dimensions: "; data.camera_grid_dimensions.print();
    std::cout << "Data size: ";
    for (auto const i : data.data.shape) {
        std::cout << i << ' ';
    }
    std::cout << '\n';
    std::cout << "Is confocal?: "; data.is_confocal.print();
    std::cout << "Camera positions: "; data.camera_point_positions.print();
#endif 
}