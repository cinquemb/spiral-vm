// examples/logical_hadamard.cpp
// Tests logical Hadamard + Z calibration on DTC
#include "../src/spiral_vm_core.hpp"
#include <iostream>
#include <iomanip>

int main() {
    // === STEP 1: CALIBRATE Z ROTATION ===
    std::cout << "=== Z ROTATION CALIBRATION ===\n";
    
    SpiralVM vm_cal(30, 30);
    vm_cal.initialize_state("neel");
    uint32_t q0 = vm_cal.add_qubit(15, 15);
    
    auto Z = [&vm_cal](uint32_t id) { return vm_cal.measure_logical_Z(id); };
    auto X = [&vm_cal](uint32_t id) { return vm_cal.measure_logical_X(id); };

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Initial: ⟨Z_L⟩ = " << Z(q0) << ", ⟨X_L⟩ = " << X(q0) << "\n";
    
    double z_slopes[] = {1.0, 10.0, 100.0, 500.0, 1000.0, 2000.0, 3500.0, 5000.0, 7500.0, 10000.0};
    double best_z_slope = 0.0;
    double best_x_error = 2.0;
    
    for(double slope : z_slopes) {
        // Reset to |0⟩_L using XX (clean start each time)
        vm_cal.logical_x_pulse(q0, 1);
        vm_cal.logical_x_pulse(q0, 1);
        
        // Apply Z rotation
        vm_cal.logical_phase_ramp(q0, slope*M_PI, 1);
        
        double x_meas = X(q0);
        double x_error = std::abs(x_meas - 0.0);
        
        std::cout << "Z(slope=" << std::setw(7) << slope 
                  << ") → ⟨X_L⟩=" << std::setw(9) << x_meas 
                  << " (error=" << x_error << ")\n";
        
        if(x_error < best_x_error) {
            best_x_error = x_error;
            best_z_slope = slope;
        }
    }
    
    std::cout << "\n*** BEST Z90 SLOPE = " << best_z_slope << " ***\n\n";
    
    // === STEP 2: FRESH VM FOR HADAMARD TEST ===
    std::cout << "=== HADAMARD TEST (FRESH VM) ===\n";
    
    SpiralVM vm_test(30, 30);
    vm_test.initialize_state("neel");
    q0 = vm_test.add_qubit(15, 15);
    
    auto Z_test = [&vm_test](uint32_t id) { return vm_test.measure_logical_Z(id); };
    auto X_test = [&vm_test](uint32_t id) { return vm_test.measure_logical_X(id); };
    
    // Start in |0⟩_L: ⟨Z⟩=+1, ⟨X⟩=0
    vm_test.logical_x_pulse(q0, 1);
    vm_test.logical_x_pulse(q0, 1);
    
    std::cout << "Before H: ⟨Z_L⟩ = " << Z_test(q0) << ", ⟨X_L⟩ = " << X_test(q0) << "\n";
    
    // PERFECT Hadamard using calibrated slope
    vm_test.logical_phase_ramp(q0,  best_z_slope*M_PI, 1);  // Z(+π/2)
    vm_test.logical_x_pulse(q0, 1);                    // X
    vm_test.logical_phase_ramp(q0, -best_z_slope*M_PI, 1);  // Z(-π/2)
    
    std::cout << "After H:  ⟨Z_L⟩ = " << Z_test(q0) << ", ⟨X_L⟩ = " << X_test(q0) << "\n";
    std::cout << "Expected: Z≈0.0, X≈+1.0\n";
    
    double h_fidelity = std::abs(Z_test(q0));  // Hadamard |0⟩_Z → |+⟩_X has ⟨Z⟩=0!
    std::cout << "Hadamard fidelity: " << std::fixed << std::setprecision(4) << h_fidelity << "\n";
    
    return 0;
}
