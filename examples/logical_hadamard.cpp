#include "../src/spiral_vm_core.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

int main() {
    const int rows = 15;
    const int cols = 15;
    const int N_steps = 200;       // incremental Hadamard steps, this is my warp limit, beyond this i next
    const double tol = 1e-3;      // convergence tolerance for X_norm
    const int max_iters = 1;     // bisection iterations

    double H_min = 1;          // min candidate amplitude
    double H_max = 1;         // max candidate amplitude
    double best_H = H_min;

    // Step-by-step Hadamard display with calibrated best_H
    SpiralVM vm(rows, cols);
    vm.is_ang = true;
    vm.initialize_state("neel");
    uint32_t q0 = vm.add_qubit((int)(rows/2), (int)(cols/2));

    auto Z = [&vm](uint32_t id){ return vm.measure_logical_Z(id); };
    auto X = [&vm](uint32_t id){ return vm.measure_logical_X(id); };

    vm.logical_x_pulse(q0, 1);
    vm.logical_x_pulse(q0, 1);

    std::cout << "\n=== Step-by-step Hadamard with H_amp=" << best_H << " ===\n";
    std::cout << "Step |    Z    |    X    | X_norm\n";
    std::cout << "--------------------------------\n";
    int kick_count = 0;

    double Z_prev = Z(q0);
    double X_prev = X(q0);

    for (int step = 0; step < N_steps; ++step) {

        //vm.logical_z_rotation(q0, 2*M_PI);
        //vm.logical_x_pulse(q0, 1.0);
        // vm.logical_y_pulse(q0, 1.0);

        /*
        int K = 20; // or ramp: 1 + step/100
        for (int k = 0; k < K; ++k){
            vm.logical_x_pulse(q0, 1.0);
            vm.logical_hadamard_step(q0, best_H, 5); // 5 subharmonics per micro-step
            vm.logical_z_rotation(q0, M_PI);
        }*/


        double K_min_large = 1.0;          // minimum K for large grids
        double K_max_2x2 = 60.0;     // maximum K for 2x2
        double x0 = 25.0;                    // sigmoid midpoint
        double steepness = 0.02;              // sigmoid steepness

        double rows_cols = static_cast<double>(vm.rows * vm.cols);  // ensure floating-point
        double max_K_mult = K_min_large + (K_max_2x2 - K_min_large) / (1.0 + std::exp(steepness * (rows_cols - x0)));


        int max_K = 25 * max_K_mult;  // full-force early
        int min_K = 5;   // later-stage micro-steps

        // ramp K down as step increases
        int K = max_K - (max_K - min_K) * step / N_steps;
        K = std::max(K, min_K);

        for (int k = 0; k < K; ++k){
            vm.logical_x_pulse(q0, 0.1);           // full X kick
            vm.logical_hadamard_step(q0, best_H, 5); // 5 subharmonics
            vm.logical_z_rotation(q0, M_PI * 0.1);       // full Z correction
        }




        // measure
        double Z_now = Z(q0);
        double X_now = X(q0);


        // update history
        Z_prev = Z_now;
        X_prev = X_now;

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
