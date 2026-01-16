#include "../src/spiral_vm_core.hpp"
#include <iostream>
#include <iomanip>
#include <random>

int main() {
    SpiralVM vm(30, 30);
    vm.initialize_state("neel");
    uint32_t q0 = vm.add_qubit(15, 15);

    auto Z = [&vm](uint32_t id){ return vm.measure_logical_Z(id); };
    auto X = [&vm](uint32_t id){ return vm.measure_logical_X(id); };

    // Prepare |0>_L
    vm.logical_x_pulse(q0, 1);
    vm.logical_x_pulse(q0, 1);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Before any H: Z = " << Z(q0) << "   X = " << X(q0) << "\n\n";

    const int N_steps = 10;
    const double H_base_candidates[] = {0.5, 1.0, 2.0, 5.0, 10.0, 20.0};

    for (double H_base : H_base_candidates) {
        std::cout << "\n=== Testing H_base = " << H_base << " ===\n";
        std::cout << "Step |    φ    |    Z    |    X    | X_norm\n";
        std::cout << "--------------------------------------------\n";

        double Z_start = Z(q0);
        double X_start = X(q0);

        for (int step = 0; step < N_steps; ++step) {
            vm.logical_hadamard_step(q0, H_base);

            double Z_now = Z(q0);
            double X_now = X(q0);
            double X_norm = X_now / (std::hypot(X_now, Z_now) + 1e-6);

            std::cout << std::setw(4) << step
                      << " | " << std::setw(8) << vm.current_orbit_phase(q0)
                      << " | " << std::setw(8) << Z_now
                      << " | " << std::setw(8) << X_now
                      << " | " << std::setw(8) << X_norm
                      << "\n";
        }

        std::cout << "Final with H_base=" << H_base << ": Z = " << Z(q0) << "   X = " << X(q0) << "\n";
    }

    return 0;
}