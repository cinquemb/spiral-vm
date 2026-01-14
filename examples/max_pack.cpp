// examples/max_pack.cpp
// Pack MAX logical qubits
#include "../src/spiral_vm_core.hpp"
#include <iostream>
#include <vector>

int main() {
    const int max_n = 100;
    SpiralVM vm(max_n, max_n);
    vm.is_ang = true;
    vm.initialize_state("neel");

    // ULTRA-DENSE
    std::vector<uint32_t> qubits;
    int x_base = max_n/2, y_base = max_n/2;
    
    // Fill every site in 5×4 block
    for(int dy = 0; dy < max_n; dy++) {
        for(int dx = 0; dx < max_n; dx++) {
            uint32_t q = vm.add_qubit(x_base+dx, y_base+dy);
            qubits.push_back(q);
        }
    }

    std::cout << "Packed " << qubits.size() << " logical qubits in " << max_n <<" x " << max_n << "phys block\n";
    
    vm.compile_to_physical_waveform();
    vm.dump_frequency_mapping();
    
    vm.run_periods(1);
    
    // Test all-to-all control: X all qubits
    // SINGLE GLOBAL PULSE - flips ALL logical qubits instantly
    vm.global_pi_pulse();  
    
    vm.run_periods(1);
    
   double avg_z = 0, neel_order = 0;
    int counted = 0;
    size_t qidx = 0;
    for(int dy = 0; dy < max_n; dy++) {
        for(int dx = 0; dx < max_n; dx++) {
            if(qidx >= qubits.size()) break;
            uint32_t q = qubits[qidx++];
            double z = vm.measure_logical_Z(q);
            if(std::abs(z) < 0.1) continue;  // skip unstable/uninitialized

            avg_z += z;

            int expected_sign = ((x_base + dx + y_base + dy) % 2) ? -1 : 1;
            neel_order += z * expected_sign;
            counted++;
        }
    }
    std::cout << "Avg Z (stable qubits) = " << avg_z/counted 
              << ", Néel order = " << neel_order/counted 
              << " over " << counted << " stable qubits\n";
    return 0;
}
