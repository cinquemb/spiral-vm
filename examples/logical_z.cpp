// examples/logical_z.cpp
// Tests logical Z rotation (phase gate) on a single logical qubit
#include "../src/spiral_vm_core.hpp"
#include <iostream>
#include <iomanip>

int main() {
    SpiralVM vm(30, 30);
    vm.initialize_state("neel");

    uint32_t q0 = vm.add_qubit(15, 15);

    std::cout << "Stabilizing...\n";
    vm.run_periods(1);

    auto Z = [&vm](uint32_t id) { return vm.measure_logical_Z(id); };

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Before Z(π/2)   → ⟨Z_L⟩ = " << Z(q0) << "\n";

    vm.logical_z_rotation(q0, M_PI/2.0);

    std::cout << "After Z(π/2)    → ⟨Z_L⟩ = " << Z(q0) << "\n";

    // Expect: Z should stay roughly the same (phase gate doesn't change populations)
    return 0;
}