// examples/logical_x.cpp
// Demonstrates logical X gate via timed global π-pulse (even cycles only)
#include "../src/spiral_vm_core.hpp"
#include <iostream>

int main() {
    SpiralVM vm(30, 30);                    // 30×30 lattice = 900 sites
    uint32_t q0 = vm.add_qubit(15, 15);    // logical qubit at center

    vm.run_periods(10);                    // settle into steady-state cat

    std::cout << "Before X: logical |0⟩ₗ (even cycle pop = "
              << vm.measure_even_population(q0) << ")\n";

    vm.apply_global_pi_pulse_on_even_cycles();  // This is logical X
    vm.run_periods(5);

    std::cout << "After X:  logical |1⟩ₗ (even cycle pop = "
              << vm.measure_even_population(q0) << ")\n";

    return 0;
}