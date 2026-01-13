// examples/ghz_3qubit.cpp
#include "../src/spiral_vm_core.hpp"
#include <iostream>

int main() {
    SpiralVM vm(30, 30);
    vm.is_ang = true;
    vm.initialize_state("neel");

    uint32_t q0 = vm.add_qubit(15, 15);  
    uint32_t q1 = vm.add_qubit(16, 15);  
    uint32_t q2 = vm.add_qubit(15, 16);  

    vm.compile_to_physical_waveform();
    vm.dump_frequency_mapping();
    
    vm.run_periods(1);
    
    std::cout << "Pre-circuit: Z0=" << vm.measure_logical_Z(q0)
              << " Z1=" << vm.measure_logical_Z(q1)
              << " Z2=" << vm.measure_logical_Z(q2) << "\n";

    // H q0 (manual - your working sequence)
    vm.logical_phase_ramp(q0, M_PI/2.0, 1);
    vm.global_pi_pulse();
    vm.logical_phase_ramp(q0, -M_PI/2.0, 1);
    
    // CNOT q0→q1: Gate(CNOT, target=q1, angle=0, control=q0)
    vm.apply_gate(Gate(Gate::CNOT, q1, 0.0, q0), 1.0);  // control=4th param
    
    // CNOT q0→q2
    vm.apply_gate(Gate(Gate::CNOT, q2, 0.0, q0), 1.0);  // control=4th param
    
    vm.run_periods(20);

    std::cout << "Post-GHZ:   Z0=" << vm.measure_logical_Z(q0)
              << " Z1=" << vm.measure_logical_Z(q1)
              << " Z2=" << vm.measure_logical_Z(q2) << "\n";
    
    std::cout << "ZZ01=" << vm.logical_zz_correlation(q0,q1) << " "
              << "ZZ02=" << vm.logical_zz_correlation(q0,q2) << " "
              << "ZZ12=" << vm.logical_zz_correlation(q1,q2) << "\n";

    return 0;
}
