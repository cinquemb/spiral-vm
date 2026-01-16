// ---------- PHASE-BOMB EXPERIMENT ----------
#include "../src/spiral_vm_core.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <random>

int main() {
    SpiralVM vm(30, 30);
    vm.initialize_state("neel");
    uint32_t q0 = vm.add_qubit(15, 15);

    auto Z = [&vm](uint32_t id){ return vm.measure_logical_Z(id); };
    auto X = [&vm](uint32_t id){ return vm.measure_logical_X(id); };

    // Start in |0>_L
    vm.logical_x_pulse(q0, 1);
    vm.logical_x_pulse(q0, 1);

    const int N_steps = 50;
    const double H_base = 1.2;     // tiny Hadamard step amplitude
    const double Z_kick_max = 0.05; // small random Z kick

    std::mt19937 rng(1234);
    std::uniform_real_distribution<double> dist(-Z_kick_max, Z_kick_max);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Step |    φ    |    Z    |    X    | Δφ(H) | Δφ(Z)\n";
    std::cout << "---------------------------------------------------------\n";

    for (int step = 0; step < N_steps; ++step) {
        double phi_before = vm.current_orbit_phase(q0);
        double Z_before = Z(q0);
        double X_before = X(q0);

        // --- 1. Small misaligned Hadamard step ---
        double i0 = H_base * (0.5 + 0.5 * std::sin(phi_before));
        double q0_amp = H_base * (0.5 + 0.5 * std::cos(phi_before));
        vm.logical_hadamard(q0, i0, q0_amp); // new step-wise H

        double phi_after_H = vm.current_orbit_phase(q0);

        // --- 2. Random small Z kick ---
        double z_kick = dist(rng);
        vm.logical_phase_ramp(q0, z_kick, 1);   // promote to logical_phase
        double phi_after_Z = vm.current_orbit_phase(q0);

        // --- 3. Record ---
        double Z_after = Z(q0);
        double X_after = X(q0);

        std::cout << std::setw(4) << step
                  << " | " << std::setw(8) << phi_before
                  << " | " << std::setw(8) << Z_before
                  << " | " << std::setw(8) << X_before
                  << " | " << std::setw(8) << (phi_after_H - phi_before)
                  << " | " << std::setw(8) << (phi_after_Z - phi_after_H)
                  << "\n";
    }

    std::cout << "\nObserve: X now grows as the logical orbit is phase-bombed.\n";
    return 0;
}
