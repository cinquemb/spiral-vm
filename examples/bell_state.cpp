// examples/bell_state.cpp
// Creates (|00⟩ₗ + |11⟩ₗ)/√2 Bell state using DTC universal gates
#include "../src/spiral_vm_core.hpp"
#include <iostream>
#include <iomanip>

int main() {
    SpiralVM vm(30, 30);
    vm.is_ang = true;
    vm.initialize_state("neel");

    uint32_t q0 = vm.add_qubit(12,15);
    uint32_t q1 = vm.add_qubit(15,15);  // farther apart to reduce overlap crosstalk (have to stay more than 2 way if sharing 1 dim index) TODO: NEED TO TUNE SIDE BANDS LATER

    vm.run_periods(5);  // short stabilization
    vm.compile_to_physical_waveform();

    auto Z0 = [&vm]() { return vm.measure_logical_Z(0); };
    auto Z1 = [&vm]() { return vm.measure_logical_Z(1); };
    auto ZZ = [&vm]() { return vm.logical_zz_correlation(0,1); };

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Initial: Z0=" << Z0() << ", Z1=" << Z1() << ", ZZ=" << ZZ() << "\n";


    vm.run_periods(5);

    // Standard Bell prep: H on q1 → CZ → X on q0
    vm.logical_hadamard(q0);       // |+⟩ on q1, |0⟩ on q0
    vm.run_periods(3);

    //vm.logical_cz(q0, q1);      // CZ (controlled from q0, target q1)
    vm.apply_phase_kick_between_full(q0,  q1,
                                             vm.LOGICAL_X_AMPLITUDE, 20); 
    vm.run_periods(3);



    std::cout << "After Bell prep:\n";
    std::cout << "  Z0 = " << Z0() << "\n";
    std::cout << "  Z1 = " << Z1() << "\n";
    std::cout << "  ZZ = " << ZZ() << " (should be negative, ~ -1 for Bell)\n";

    // Optional: flip q0 again and check anti-correlation persists
    vm.logical_x_pulse(q0, 1.0);
    vm.run_periods(3);
    std::cout << "After extra X on q0: Z0=" << Z0() << ", Z1=" << Z1() << ", ZZ=" << ZZ() << "\n";

    return 0;
}