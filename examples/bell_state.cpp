// examples/bell_state.cpp
// Creates |00⟩ₗ + |11⟩ₗ using spiral phase kick between two logical qubits
#include "../src/spiral_vm_core.hpp"
#include <iostream>

int main() {
    SpiralVM vm(30, 30);
    vm.is_ang = true;
    //vm.initialize_state("polarized");
    //vm.initialize_state("neel");
    vm.initialize_state("neel");

    uint32_t q0 = vm.add_qubit(12, 15);
    uint32_t q1 = vm.add_qubit(18, 15);    // two neighboring spirals

    vm.compile_to_physical_waveform();  // Now only 1 waveform, broadcast to all
    vm.dump_frequency_mapping();// dump mapping

    vm.run_periods(20);

     // MEASURE INITIAL: Both |0⟩ₗ (Néel even sites)
    std::cout << "Pre-CZ: Z0=" << vm.measure_logical_Z(q0) 
              << " Z1=" << vm.measure_logical_Z(q1) << "\n";


    // Logical Hadamard on q0
    vm.logical_phase_ramp(q0, M_PI/2.0, 1);  // Z(π/2) over 20 steps
    vm.global_pi_pulse();  // X
    vm.logical_phase_ramp(q0, -M_PI/2.0, 1);  // Z(-π/2)

    // Logical CZ via short phase-gradient kick between spirals
    vm.apply_phase_kick_between(q0, q1, 0.25, 0.15*vm.T);  // strength & duration

    vm.run_periods(10);

    // MEASURE FINAL: Perfect anticorrelation = Bell state
    std::cout << "Post-CZ: Z0=" << vm.measure_logical_Z(q0) 
              << " Z1=" << vm.measure_logical_Z(q1) << "\n";
    std::cout << "ZZ corr=" << vm.logical_zz_correlation(q0, q1) 
              << " (1.0 = perfect |00⟩+|11⟩)\n";

    return 0;
}