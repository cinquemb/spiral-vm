// examples/t_gate_manual.cpp
#include "../src/spiral_vm_core.hpp"
#include <iostream>
#include <cmath> // for M_PI

int main() {
    SpiralVM vm(30, 30);
    vm.is_ang = true;
    vm.initialize_state("neel");

    uint32_t q0 = vm.add_qubit(15, 15);

    vm.compile_to_physical_waveform();
    vm.dump_frequency_mapping();
        
    // Let the system settle
    vm.run_periods(1);

    // Read initial logical phase
    double phase_before = vm.get_logical_phase_frame_corrected(q0);
    std::cout << "Before T: phase = " << phase_before << "\n";

        double target_delta = M_PI/4.0;

    vm.apply_gate(Gate(Gate::T, q0), 1.0);  // control=4th param

    // Run a few periods to let ramp finish
    vm.run_periods(1);

    // Read final logical phase
    double phase_after = vm.get_logical_phase_frame_corrected(q0);
    std::cout << "After T (manual ramp): phase = " << phase_after
              << " (expect ~ " << phase_before + target_delta << ")\n";

    return 0;
}
