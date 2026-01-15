// examples/logical_hadamard_phaseaware_test.cpp
// Phase-aware Hadamard tester using orbit phase

#include "../src/spiral_vm_core.hpp"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

struct ZCalResult {
    double slope;
    double X;
    double Z;
};

// Phase-aware Z calibration
ZCalResult test_z_slope(double slope) {
    SpiralVM vm(30, 30);
    vm.initialize_state("neel");
    uint32_t q0 = vm.add_qubit(15, 15);

    auto Z = [&vm](uint32_t id){ return vm.measure_logical_Z(id); };
    auto X = [&vm](uint32_t id){ return vm.measure_logical_X(id); };

    // Prepare |0>_L
    vm.logical_x_pulse(q0, 1);
    vm.logical_x_pulse(q0, 1);

    double phi0 = vm.current_orbit_phase(q0);

    // Apply Z rotation in small increments
    int steps = 100;
    double dphi = slope * M_PI / steps;
    for(int s=0; s<steps; ++s)
        vm.logical_phase_ramp(q0, dphi, 1);

    double phi1 = vm.current_orbit_phase(q0);

    std::cout << "    phase: " << phi0
              << " -> " << phi1
              << "  Δφ=" << (phi1 - phi0) << "\n";

    return {slope, X(q0), Z(q0)};
}

int main() {
    std::cout << std::fixed << std::setprecision(9);

    // --- Step 1: Z calibration ---
    std::cout << "=== Z ROTATION CALIBRATION ===\n";
    std::vector<double> z_slopes = {1.0, 10.0, 100.0, 500.0, 1000.0,
                                    2000.0, 3500.0, 5000.0, 7500.0, 10000.0};
    double best_slope = 0.0;
    double best_err = 1e9;

    for(double s : z_slopes) {
        auto r = test_z_slope(s);
        double err = std::abs(r.X); // target X≈0
        std::cout << "Z(slope=" << std::setw(7) << s
                  << ")  X=" << std::setw(12) << r.X
                  << "  Z=" << std::setw(12) << r.Z
                  << "  |X|=" << err << "\n";

        if(err < best_err) {
            best_err = err;
            best_slope = s;
        }
    }
    std::cout << "\n*** BEST Z90 SLOPE = " << best_slope
              << " (|X|=" << best_err << ") ***\n\n";

    // --- Step 2: Hadamard test ---
    std::cout << "=== PHASE-AWARE HADAMARD TEST ===\n";
    std::cout << "Trial |    φ0    |    Z0    |    X0    |    Z1    |    X1    |    φ1\n";
    std::cout << "-----------------------------------------------------------------------\n";

    for(int trial=0; trial<5; ++trial) {
        SpiralVM vm(30, 30);
        vm.initialize_state("neel");
        uint32_t q0 = vm.add_qubit(15, 15);

        auto Z = [&vm](uint32_t id){ return vm.measure_logical_Z(id); };
        auto X = [&vm](uint32_t id){ return vm.measure_logical_X(id); };

        // Prepare |0>_L
        vm.logical_x_pulse(q0, 1);
        vm.logical_x_pulse(q0, 1);

        // Start at a different logical orbit phase
        double random_phi = (trial * M_PI)/3.0;
        int steps = 50;
        double dphi = random_phi / steps;
        for(int s=0; s<steps; ++s)
            vm.logical_phase_ramp(q0, dphi, 1);

        double phi0 = vm.current_orbit_phase(q0);
        double Z0 = Z(q0);
        double X0 = X(q0);

        vm.logical_hadamard(q0);

        double phi1 = vm.current_orbit_phase(q0);
        double Z1 = Z(q0);
        double X1 = X(q0);

        std::cout << std::setw(5) << trial << " | "
                  << std::setw(8) << phi0 << " | "
                  << std::setw(8) << Z0   << " | "
                  << std::setw(8) << X0   << " | "
                  << std::setw(8) << Z1   << " | "
                  << std::setw(8) << X1   << " | "
                  << std::setw(8) << phi1 << "\n";
    }

    std::cout << "\nExpected after H: Z≈0, X≈+1 (independent of φ)\n";

    return 0;
}
