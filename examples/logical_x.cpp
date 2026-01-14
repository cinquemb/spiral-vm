// examples/logical_x_loop.cpp
#include "../src/spiral_vm_core.hpp"
#include <iostream>
#include <iomanip>

int main() {
    SpiralVM vm(30, 30);
    vm.is_ang = true;
    vm.initialize_state("neel");

    uint32_t q0 = vm.add_qubit(15, 15);

    std::cout << "Stabilizing...\n";
    vm.run_periods(1);  // extra stabilization

    auto Z = [&vm](uint32_t id) { return vm.measure_logical_Z(id); };

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\n=== LOGICAL X GATE LOOP TEST ===\n";
    
    // Flip flop 10 times
    for(int cycle = 0; cycle < 10; cycle++) {
        double z_before = Z(q0);
        std::cout << "\nCycle " << cycle << " START → ⟨Z_L⟩ = " << z_before << "\n";
        
        vm.logical_x_pulse(q0, 1);  // X gate
        vm.run_periods(1);          // relax
        
        double z_after = Z(q0);
        std::cout << "Cycle " << cycle << "  X  → ⟨Z_L⟩ = " << z_after << "\n";
        
        // Check flip quality
        double flip_fidelity = std::abs(z_before * z_after);
        std::cout << "Flip fidelity: " << std::fixed << std::setprecision(4) 
                  << flip_fidelity << " (target: 1.0000)\n";
        
        if(flip_fidelity < 0.99) {
            std::cout << "WARNING: Poor flip fidelity!\n";
        }
    }
    
    std::cout << "\n=== LOOP COMPLETE ===\n";
    return 0;
}
