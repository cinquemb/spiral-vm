// examples/t_gate.cpp
// Implements non-Clifford T-gate via asymmetric twist ramp
#include "../src/spiral_vm_core.hpp"
#include <iostream>

int main() {
    SpiralVM vm(30, 30);
    uint32_t q0 = vm.add_qubit(15, 15);

    vm.run_periods(15);

    std::cout << "Before T: phase = " << vm.get_logical_phase(q0) << "\n";

    // T-gate = π/4 phase → ramp omega_ang asymmetrically over 8 periods
    vm.ramp_omega_ang(126.0, 138.0, 8*T);   // +12 rad/s over 8 cycles

    vm.run_periods(10);

    std::cout << "After T:  phase = " << vm.get_logical_phase(q0)
              << "  (expected +π/4)\n";

    return 0;
}