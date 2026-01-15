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

    auto Z = [&vm](uint32_t id) { return vm.measure_logical_Z(id); };

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\n=== ADAPTIVE LOGICAL X (100 cycles) ===\n";
    
    double k = 0.1;  // Conservative learning rate
    
    std::cout << "Initial amp=" << vm.LOGICAL_X_AMPLITUDE << "\n";
    //vm.LOGICAL_X_AMPLITUDE *= 0.5;
    for(int cycle = 0; cycle < 50000; cycle++) {
        double z_before = Z(q0);
        
        vm.logical_x_pulse(q0, 1);  // Uses vm.LOGICAL_X_AMPLITUDE/2
        //vm.logical_phase_ramp(q0, 2.0*M_PI, 1);
        //vm.logical_x_pulse(q0, 1);  // Uses vm.LOGICAL_X_AMPLITUDE/2
        
        double z_after = Z(q0);
        double flip_fidelity = std::abs(z_before * z_after);
        
        // Adaptive update: target = -z_before
        double target_z = -z_before;
        double err = z_after - target_z;
        
        // UPDATE THE VM'S AMPLITUDE
        //vm.LOGICAL_X_AMPLITUDE += (vm.get_period() % 2) ? 0 : -20;
        
        std::cout << "Cycle " << cycle 
                  << ": " << z_before << " → " << z_after 
                  << " fid=" << std::fixed << std::setprecision(4) << flip_fidelity
                  << " amp=" << vm.LOGICAL_X_AMPLITUDE << "\n";
        
        if(flip_fidelity < 0.97) {
            std::cout << "WARNING: Poor fidelity!\n";
            exit(0);
        }
        
        /*// Converged?
        if(std::abs(err) < 0.005 && flip_fidelity > 0.9995) {
            std::cout << "\n*** CONVERGED! Final amp = " << vm.LOGICAL_X_AMPLITUDE << " ***\n";
            break;
        }*/
    }
    
    return 0;
}
