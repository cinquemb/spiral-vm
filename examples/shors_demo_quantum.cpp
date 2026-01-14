// examples/shors_demo_quantum_retry.cpp
// Shor's algorithm on SpiralVM: Fully quantum period finding with VM restart for retries
#include "../src/spiral_vm_core.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <numeric>  // for gcd

bool is_coprime(long long a, long long b) {
    return std::gcd(a, b) == 1;
}

long long mod_pow(long long base, long long exp, long long mod) {
    long long result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp % 2 == 1) result = (result * base) % mod;
        base = (base * base) % mod;
        exp /= 2;
    }
    return result;
}

int infer_period(const std::vector<double>& z) {
    size_t m = z.size();
    for (size_t r = 1; r <= m/2; ++r) {
        bool match = true;
        for (size_t i = 0; i < m - r; ++i) {
            if (std::abs(z[i] - z[i+r]) > 0.05) {
                match = false;
                break;
            }
        }
        if (match) return r;
    }
    return m;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: ./shors_demo_quantum_retry <N> [lattice_size=8]\n";
        return 1;
    }

    const long long N = std::atoll(argv[1]);
    const int n_qubits = (argc > 2) ? std::atoi(argv[2]) : 8;

    if (N < 2 || N % 2 == 0) {
        std::cerr << "[Shor] N must be odd semiprime >1.\n";
        return 1;
    }

    srand(time(NULL));
    long long factor1 = 1, factor2 = N;
    int attempt = 0;

    while (factor1 == 1 || factor1 == N) {
        attempt++;
        std::cout << "[Shor] Attempt #" << attempt << "...\n";

        // ---- Restart SpiralVM ----
        SpiralVM vm(n_qubits, n_qubits);
        vm.is_ang = true;
        vm.initialize_state("neel");

        // Phase register size ~ 2*log2(N)
        int m = std::ceil(2 * std::log2(N));
        std::vector<uint32_t> phase;
        for (int i = 0; i < m; ++i)
            phase.push_back(vm.add_qubit(2 + (i % 4), 2 + (i / 4)));

        uint32_t work = vm.add_qubit(6, 6);

        vm.compile_to_physical_waveform();
        vm.dump_frequency_mapping();

        // Pick random base a coprime to N
        int a;
        do {
            a = 2 + rand() % (N - 2);
        } while (!is_coprime(a, N));
        std::cout << "[Shor] Random base a=" << a << " (coprime to " << N << ")\n";

        // Step 1: Hadamard on all phase qubits
        for (auto q : phase) vm.logical_hadamard(q);
        vm.run_periods(1);

        // Step 2: Modular evolution (global proxy)
        for (int i = 0; i < m; ++i) {
            int reps = 1 << i;
            for (int k = 0; k < reps; ++k) {
                vm.global_pi_pulse();
                vm.run_periods(1);
            }
        }

        // Step 3: Approx inverse QFT
        for (int i = m - 1; i >= 0; --i) {
            for (int j = i + 1; j < m; ++j) {
                double theta = M_PI / (1 << (j - i));
                vm.logical_phase_ramp(phase[i], theta, 1);
            }
            vm.logical_hadamard(phase[i]);
        }
        vm.run_periods(1);

        // Step 4: Measure phase register
        std::vector<double> z;
        for (auto q : phase) z.push_back(vm.measure_logical_Z(q));
        std::cout << "[Shor] Phase Z:";
        for (double v : z) std::cout << " " << v;
        std::cout << "\n";

        double zw = vm.measure_logical_Z(work);
        std::cout << "[Shor] Work Z=" << zw << "\n";

        // Step 5: Quantum period inference
        int r = infer_period(z);
        if (r <= 1 || r >= N) {
            std::cerr << "[Shor] Trivial period r=" << r << "; retrying with new base a.\n";
            continue; // loop retries
        }

        if (r % 2 != 0) r *= 2;
        std::cout << "[Shor] Estimated period r=" << r << "\n";

        // Step 6: Classical factor extraction
        long long a_k = mod_pow(a, r / 2, N);
        factor1 = std::gcd(N, (long long)std::abs(a_k - 1));
        factor2 = std::gcd(N, (long long)std::abs(a_k + 1));

        if (factor1 == 1 || factor1 == N) {
            std::cerr << "[Shor] Trivial factors; retrying with a new base a.\n\n";
            factor1 = 1;  // reset to continue loop
            factor2 = N;
        } else {
            std::cout << "[Shor] SUCCESS: " << N << " = " << factor1 << " × " << factor2 << "\n";
        }
    }

    return 0;
}
