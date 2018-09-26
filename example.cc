#include "nlos_reader.h"
#include <chrono>

int main(int argc, char* argv[]) {
    nlos_dataset a(argv[1]);
    #ifdef USE_XTENSOR
    std::cout << "Bin resolution: " << a.deltat << '\n' <<
                 "Bins: " << a.bins << '\n' << 
                 "Cam grid dimensions: " << a.camera_grid_dimensions << '\n' << 
                 "Data size: ";
    for (auto const i :  a.data.shape()) {
        std::cout << i << ' ';
    }
    std::cout << '\n';
    std:: cout << "Is confocal?: " << a.is_confocal << '\n' << 
    "Camera positions: " << a.camera_point_positions;
#else 
    std::cout << "Bin resolution: "; a.deltat.print();
    std::cout << "Bins: "; a.bins.print();
    std::cout << "Cam grid dimensions: "; a.camera_grid_dimensions.print();
    std::cout << "Data size: ";
    for (auto const i : a.data.shape) {
        std::cout << i << ' ';
    }
    std::cout << '\n';
    std::cout << "Is confocal?: "; a.is_confocal.print();
    std::cout << "Camera positions: "; a.camera_point_positions.print();
#endif 
}