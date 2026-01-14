// examples/logical_t.cpp
// Tests logical T gate (phase π/4 on |1⟩ branch) using superposition
#include "../src/spiral_vm_core.hpp"
#include <iostream>
#include <iomanip>

int main() {
    SpiralVM vm(30, 30);
    vm.initialize_state("polarized");

    uint32_t q0 = vm.add_qubit(15, 15);

    std::cout << "Stabilizing...\n";
    vm.run_periods(1);

    // Create |+⟩ = H|0⟩
    vm.logical_hadamard(q0);
    vm.run_periods(1);

    auto Z = [&vm](uint32_t id) { return vm.measure_logical_Z(id); };
    auto Y = [&vm](uint32_t id) { return vm.measure_logical_Y(id); };

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Before T       → ⟨Z_L⟩ = " << Z(q0) << ", ⟨Y_L⟩ = " << Y(q0) << "\n";

    vm.logical_phase_ramp(q0, M_PI / 4.0, 1);  // T gate
    vm.run_periods(1);

    std::cout << "After T        → ⟨Z_L⟩ = " << Z(q0) << ", ⟨Y_L⟩ = " << Y(q0) << "\n";

    // Expect: on |+⟩, T should rotate phase → ⟨Y_L⟩ should become positive (~0.707 ideal)
    //For T gate: Check ⟨Y⟩ after H+T+H (should give S gate behavior).
    return 0;
}