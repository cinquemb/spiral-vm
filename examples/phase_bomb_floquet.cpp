// examples/phase_bomb_floquet.cpp
// Deliberately phase-bombs a logical qubit orbit
// Tracks φ, Z, X for each pulse step

#include "../src/spiral_vm_core.hpp"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <random>

int main() {
    std::cout << std::fixed << std::setprecision(9);

    SpiralVM vm(30, 30);
    vm.initialize_state("neel");
    uint32_t q0 = vm.add_qubit(15, 15);

    auto Z = [&vm](uint32_t id){ return vm.measure_logical_Z(id); };
    auto X = [&vm](uint32_t id){ return vm.measure_logical_X(id); };
    auto φ = [&vm](uint32_t id){ return vm.current_orbit_phase(id); };

    // Start in |0>_L
    vm.logical_x_pulse(q0,1);
    vm.logical_x_pulse(q0,1);

    // Random number generator for small misalignments
    std::mt19937 rng(1234);
    std::uniform_real_distribution<double> dphi_dist(-0.05,0.05); // ±0.05 rad
    std::uniform_real_distribution<double> dz_dist(-0.02,0.02);   // ±0.02 rad Z kick

    int N_steps = 50;

    std::cout << "Step |    φ    |    Z    |    X    | Δφ(H) | Δφ(Z)\n";
    std::cout << "---------------------------------------------------------\n";

    for(int step=0; step<N_steps; ++step) {
        double φ0 = φ(q0);
        double Z0 = Z(q0);
        double X0 = X(q0);

        // --- 1. Misaligned Hadamard ---
        double H_misalignment = dphi_dist(rng);
        double base_amp = 161.1/std::sqrt(2.0);
        double i0 = base_amp * std::cos(φ0 + M_PI/4.0 + H_misalignment);
        double q0_i = base_amp * std::sin(φ0 + M_PI/4.0 + H_misalignment);

        int wid = vm.logical_qubits[q0].waveform_id;
        if(wid >= 0 && wid < (int)vm.waveforms.size()) {
            size_t carrier_idx = static_cast<size_t>(-1);
            for(size_t i=0;i<vm.waveforms[wid].tones.size();++i){
                auto &tn = vm.waveforms[wid].tones[i];
                if(tn.logical_id == (int)q0 && std::abs(tn.freq - vm.allocated_carriers[q0])<1e-6){
                    carrier_idx = i; break;
                }
            }
            if(carrier_idx!=static_cast<size_t>(-1)){
                vm.waveforms[wid].tones[carrier_idx].set_iq(i0,q0_i);
            }
        }

        vm.compile_to_physical_waveform();
        vm.run_periods(1.0);

        // Track Hadamard phase increment (π/2)
        vm.logical_phase[q0] = std::fmod(vm.logical_phase[q0] + M_PI/2.0 + H_misalignment, 2*M_PI);

        double ΔφH = M_PI/2.0 + H_misalignment;

        // --- 2. Small random Z kick ---
        double Z_kick = dz_dist(rng);
        vm.logical_phase_ramp(q0, Z_kick, 1);

        double ΔφZ = Z_kick;

        // --- 3. Record after step ---
        double φ1 = φ(q0);
        double Z1 = Z(q0);
        double X1 = X(q0);

        std::cout << std::setw(4)<<step<<" | "
                  << std::setw(8)<<φ1<<" | "
                  << std::setw(8)<<Z1<<" | "
                  << std::setw(8)<<X1<<" | "
                  << std::setw(6)<<ΔφH<<" | "
                  << std::setw(6)<<ΔφZ<<"\n";
    }

    std::cout << "\n* Floquet phase-bomb complete! Watch φ, Z, X spiral.\n";

    return 0;
}
