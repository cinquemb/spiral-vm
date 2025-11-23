// examples/bell_state.cpp
// Creates |00⟩ₗ + |11⟩ₗ using spiral phase kick between two logical qubits
#include "../src/spiral_vm.hpp"
#include <iostream>

int main() {
    SpiralVM vm(30, 30);
    uint32_t q0 = vm.add_qubit(12, 15);
    uint32_t q1 = vm.add_qubit(18, 15);    // two neighboring spirals

    vm.run_periods(20);

    // Logical Hadamard on q0
    vm.apply_global_pi_pulse_on_even_cycles();     // X
    vm.apply_phase_shift(M_PI/2);                  // S gate → H = S·X·S†

    // Logical CZ via short phase-gradient kick between spirals
    vm.apply_phase_kick_between(q0, q1, 0.18, 0.3*T);  // strength & duration

    vm.run_periods(10);

    double corr = vm.logical_zz_correlation(q0, q1);
    std::cout << "Logical ZZ correlation = " << corr
              << "  (≈1.0 = perfect Bell state)\n";

    return 0;
}