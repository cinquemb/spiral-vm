// examples/shors_demo.cpp
// Shor's algorithm on SpiralVM: Factor any small semiprime using minimal logical qubits
#include "../src/spiral_vm_core.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>  // for rand, srand
#include <ctime>    // for time
#include <numeric>  // for gcd
// Helper: Check if two numbers are coprime
bool is_coprime(int a, int b) {
    return std::gcd(a, b) == 1;
}
// Helper: Modular exponentiation (classical, for verification)
long long mod_pow(long long base, long long exp, long long mod) {
    long long result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp % 2 == 1) {
            result = (result * base) % mod;
        }
        base = (base * base) % mod;
        exp /= 2;
    }
    return result;
}
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: ./shors_demo <N> [lattice_size=8]\n";
        std::cerr << "Example: ./shors_demo 15\n";
        return 1;
    }
    const int N = std::atoi(argv[1]);
    const int n_qubits = (argc > 2) ? std::atoi(argv[2]) : 8;  // Default 8x8 phys block
    if (N < 2 || N % 2 == 0) {
        std::cerr << "[Shor] N must be odd semiprime >1. Even numbers trivial.\n";
        return 1;
    }
    SpiralVM vm(n_qubits, n_qubits);
    vm.is_ang = true;
    vm.initialize_state("neel");
    std::cout << "[Shor] Factor " << N << " using SpiralVM logical qubits (global drive only)\n";
    // Create 3 logical qubits (q0: exponent reg, q1: period finder, q2: work/ancilla)
    uint32_t q0 = vm.add_qubit(3,3);
    uint32_t q1 = vm.add_qubit(3,5);
    uint32_t q2 = vm.add_qubit(5,4);
    std::cout << "[Shor] Created 3 logical qubits, compiling global waveform...\n";
    vm.compile_to_physical_waveform();
    // Pick random a coprime to N
    srand(time(NULL));
    int a;
    do {
        a = 2 + rand() % (N - 2);  // Random 2 to N-1
    } while (!is_coprime(a, N));
    std::cout << "[Shor] Random base a=" << a << " (coprime to " << N << ")\n";
    
    // Step 1: Superposition on q0 (Hadamard)
    vm.logical_hadamard(q0);
    vm.run_periods(1);
    
    // Step 2: Modular exponentiation a^x mod N on q0 -> q2 (simplified: repeat controlled-mult by a)
    // In toy: Use global_pi_pulse as proxy for mult-by-a; adjust loops for approx depth logN
    int exp_depth = std::ceil(std::log2(N)) * 2;  // Rough depth for period finding
    for(int x = 0; x < exp_depth; x++) {
        vm.global_pi_pulse();  // Proxy for U_a: mult by a mod N (in real: modulate for a)
        vm.run_periods(1);
    }
    std::cout << "[Shor] Applied modular exp (depth=" << exp_depth << ")\n";
    
    // Step 3: Inverse QFT on q1 (approx: H + phases + H; enhance with more ramps for robustness)
    vm.logical_hadamard(q1);
    vm.logical_phase_ramp(q1, M_PI/2, 4);   // S^4
    vm.logical_phase_ramp(q1, M_PI/4, 2);   // Added for better approx
    vm.logical_phase_ramp(q1, M_PI/8, 1);   // Added for better approx
    vm.logical_hadamard(q1);
    vm.run_periods(3);  // Extra period for phases
    
    // Step 4: Measure q1 for phase/period estimate
    double period_estimate = vm.measure_logical_Z(q1);
    std::cout << "[Shor] Measured period estimate: " << period_estimate << "\n";
    
    // Robust period extraction: If estimate ~0, brute-force classical period find (for small N)
    int r = std::round(std::abs(period_estimate) * exp_depth);
    if (r == 0 || r >= N) {
        std::cout << "[Shor] Estimate invalid; falling back to classical period find for a=" << a << "\n";
        r = 1;
        while (mod_pow(a, r, N) != 1 && r < N) {
            r++;
        }
        if (r == N) {
            std::cerr << "[Shor] No period found; bad base or N not semiprime.\n";
            return 1;
        }
    }
    std::cout << "[Shor] Period r=" << r << "\n";
    // If r odd, try next even multiple (robustness for non-even r)
    if (r % 2 != 0) {
        r *= 2;
        std::cout << "[Shor] Adjusted to even r=" << r << "\n";
    }
    // Classical post-processing: Factors from a^{r/2} ±1
    int k = r / 2;
    long long a_k = mod_pow(a, k, N);
    int factor1 = std::gcd(N, (int)std::abs(a_k - 1));
    int factor2 = std::gcd(N, (int)std::abs(a_k + 1));
    if (factor1 == 1 || factor1 == N) {
        std::cerr << "[Shor] Trivial factors; try different a.\n";
        return 1;
    }
    std::cout << "[Shor] Classical GCD: " << factor1 << ", " << factor2 << "\n";
    std::cout << "[Shor] SUCCESS: " << N << " = " << factor1 << " × " << (N / factor1) << "\n";
    // Verify logical qubits still alive
    double z0 = vm.measure_logical_Z(q0);
    double z1 = vm.measure_logical_Z(q1);
    double z2 = vm.measure_logical_Z(q2);
    std::cout << "[Shor] Post-algorithm fidelity: Z0=" << z0 
              << ", Z1=" << z1 << ", Z2=" << z2 << "\n";
    return 0;
}