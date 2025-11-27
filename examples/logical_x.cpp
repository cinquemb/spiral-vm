// examples/logical_x.cpp
// Demonstrates logical X gate via timed global π-pulse (even cycles only)
#include "../src/spiral_vm_core.hpp"
#include <iostream>
#include <iomanip>

int main() {
    SpiralVM vm(30, 30);
    vm.is_ang = true;
    //vm.initialize_state("polarized");
    //vm.initialize_state("neel");
    vm.initialize_state("disordered");

    uint32_t q0 = vm.add_qubit(15, 15);

    std::cout << "Stabilizing DTC...\n";
    vm.run_periods(25);   // let the cat form

    auto Z = [&vm](uint32_t id) { return vm.measure_logical_Z(id); };

    std::cout << std::fixed << std::setprecision(9);
    std::cout << "Before logical X → ⟨Z_L⟩ = " << Z(q0) << "\n";

    vm.global_pi_pulse();        // ← LOGICAL X ON THE ENTIRE LATTICE
    vm.run_periods(5);           // relax

    std::cout << "After logical X  → ⟨Z_L⟩ = " << Z(q0) << "\n";

    return 0;
}