#include "../src/spiral_vm_core.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

int main() {
    const int rows = 30;
    const int cols = 30;
    const int N_steps = 500;       // incremental Hadamard steps
    const double tol = 1e-3;      // convergence tolerance for X_norm
    const int max_iters = 1;     // bisection iterations

    double H_min = 1;          // min candidate amplitude
    double H_max = 1;         // max candidate amplitude
    double best_H = H_min;

    /*
    // Calibration loop: find H_amp that rotates X to ~1
    for (int iter = 0; iter < max_iters; ++iter) {
        double H_mid = 0.5 * (H_min + H_max);

        // New VM for each trial
        SpiralVM vm(rows, cols);
        vm.initialize_state("neel");
        uint32_t q0 = vm.add_qubit(15, 15);

        auto Z = [&vm](uint32_t id){ return vm.measure_logical_Z(id); };
        auto X = [&vm](uint32_t id){ return vm.measure_logical_X(id); };

        // Initialize |0>_L properly
       // vm.logical_x_pulse(q0, 1);
        //vm.logical_x_pulse(q0, 1);

        // Apply incremental Hadamard
        for (int step = 0; step < N_steps; ++step)
            vm.logical_hadamard_step(q0, H_mid);

        double X_now = X(q0);
        double Z_now = Z(q0);
        double X_norm = X_now / (std::hypot(X_now, Z_now) + 1e-12);

        std::cout << "Iter " << iter
                  << " | H_amp=" << H_mid
                  << " | X_norm=" << X_norm
                  << " | Z=" << Z_now << "\n";

        if (std::abs(X_norm - 1.0) < tol) {
            best_H = H_mid;
            break;
        }

        if (X_norm < 1.0) {
            H_min = H_mid; // too small → increase
        } else {
            H_max = H_mid; // too large → decrease
        }

        best_H = H_mid;
    }*/

    // Step-by-step Hadamard display with calibrated best_H
    SpiralVM vm(rows, cols);
    vm.initialize_state("neel");
    uint32_t q0 = vm.add_qubit(15, 15);

    auto Z = [&vm](uint32_t id){ return vm.measure_logical_Z(id); };
    auto X = [&vm](uint32_t id){ return vm.measure_logical_X(id); };

    vm.logical_x_pulse(q0, 1);
    vm.logical_x_pulse(q0, 1);

    std::cout << "\n=== Step-by-step Hadamard with H_amp=" << best_H << " ===\n";
    std::cout << "Step |    Z    |    X    | X_norm\n";
    std::cout << "--------------------------------\n";

    for (int step = 0; step < N_steps; ++step) {
        vm.logical_hadamard_step(q0, best_H);

        double X_now = X(q0);
        double Z_now = Z(q0);
        double X_norm = X_now / (std::hypot(X_now, Z_now) + 1e-12);

        std::cout << std::setw(4) << step
                  << " | " << std::setw(8) << Z_now
                  << " | " << std::setw(8) << X_now
                  << " | " << std::setw(8) << X_norm
                  << "\n";
    }

    std::cout << "\nFinal: Z=" << Z(q0)
              << "  X=" << X(q0)
              << "  X_norm=" << X(q0)/std::hypot(X(q0), Z(q0)) << "\n";

    return 0;
}
