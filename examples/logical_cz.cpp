// examples/logical_cz.cpp
// Basic test of controlled-Z: apply to Bell-like state
#include "../src/spiral_vm_core.hpp"
#include <iostream>
#include <iomanip>

int main() {
    SpiralVM vm(30, 30);
    vm.initialize_state("polarized");

    uint32_t ctrl = vm.add_qubit(12, 15);
    uint32_t targ = vm.add_qubit(18, 15);  // separated enough

    std::cout << "Stabilizing...\n";
    vm.run_periods(25);

    auto ZZ = [&vm](uint32_t a, uint32_t b) { return vm.logical_zz_correlation(a, b); };

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Before CZ      → ⟨Z_c Z_t⟩ = " << ZZ(ctrl, targ) << "\n";

    // Create |Φ+⟩-like: H on target, then CZ
    vm.logical_hadamard(targ);
    vm.run_periods(3);

    vm.logical_cz(ctrl, targ);          // apply CZ
    vm.run_periods(5);

    std::cout << "After CZ       → ⟨Z_c Z_t⟩ = " << ZZ(ctrl, targ) << "\n";

    // Expect: if ctrl in |1⟩ subspace, should see anti-correlation (negative ZZ)
    // For full Bell test you'd need more prep, but this shows conditional phase
    //For CZ: Tune zz_strength scaling in logical_controlled_phase until you get strong negative ZZ correlation in |11⟩ subspace.
    return 0;
}