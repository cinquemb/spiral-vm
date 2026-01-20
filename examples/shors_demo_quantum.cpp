// examples/shors_demo_quantum_retry.cpp
// Shor's algorithm on SpiralVM: Fully quantum period finding with VM restart for retries
#include "../src/spiral_vm_core.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <numeric>  // for std::gcd (C++17+)

bool is_coprime(long long a, long long b) {
    return std::gcd(static_cast<unsigned long long>(a), static_cast<unsigned long long>(b)) == 1;
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
            if (std::abs(z[i] - z[i+r]) > 0.05) {  // tolerance for noise
                match = false;
                break;
            }
        }
        if (match) return static_cast<int>(r);
    }
    return static_cast<int>(m);  // fallback: no clear period
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: ./shors_demo_quantum_retry <N> [lattice_size=8] [use_phi_direct=1]\n";
        std::cerr << "  N = number to factor (odd semiprime >1)\n";
        std::cerr << "  lattice_size = square lattice side (default 8)\n";
        std::cerr << "  use_phi_direct = 1 (fast phi mode) or 0 (waveform mode), default 1\n";
        return 1;
    }

    const long long N = std::atoll(argv[1]);
    const int lattice_size = (argc > 2) ? std::atoi(argv[2]) : 8;
    const bool use_phi = (argc > 3) ? (std::atoi(argv[3]) != 0) : true;

    if (N < 3 || N % 2 == 0) {
        std::cerr << "[Shor] N must be odd semiprime >2.\n";
        return 1;
    }

    srand(static_cast<unsigned>(time(nullptr)));

    long long factor1 = 1, factor2 = N;
    int attempt = 0;

    while (factor1 == 1 || factor1 == N) {
        attempt++;
        std::cout << "\n[Shor] Attempt #" << attempt << " on N=" << N << "\n";

        // ---- Fresh VM instance ----
        SpiralVM vm(lattice_size, lattice_size);
        vm.is_ang = true;
        vm.auto_compile_enabled = false;
        vm.use_phi_direct = use_phi;  // <--- the toggle you wanted
        vm.initialize_state("neel");

        // Phase register size ~ 2*log2(N)
        int m = static_cast<int>(std::ceil(2 * std::log2(N)));
        std::vector<uint32_t> phase;
        for (int i = 0; i < m; ++i) {
            phase.push_back(vm.add_qubit(2 + (i % 4), 2 + (i / 4)));
        }

        uint32_t work = vm.add_qubit(6, 6);

        vm.compile_to_physical_waveform();
        vm.dump_frequency_mapping("freq_map_attempt_" + std::to_string(attempt) + ".json");

        // Pick random base a coprime to N
        long long a;
        do {
            a = 2 + (rand() % (N - 2));
        } while (!is_coprime(a, N));
        std::cout << "[Shor] Random base a=" << a << " (coprime to " << N << ")\n";

        // Step 1: Hadamard on all phase qubits (superposition)
        for (auto q : phase) {
            vm.logical_hadamard(q, 200, 1.0);  // full Hadamard calibration
        }
        vm.run_periods(1);  // let the state settle

        // Step 2: Modular exponentiation proxy (repeated X pulses scaled by power)
        // This is a very crude approximation — real Shor needs controlled-U^{2^k}
        // Here we just apply global X pulses proportional to 2^k (not controlled)
        for (int k = 0; k < m; ++k) {
            int reps = 1 << k;  // 2^k
            std::cout << "[Shor] Applying 2^" << k << " = " << reps << " X pulses\n";
            for (int r = 0; r < reps; ++r) {
                for (auto q : phase) {
                    vm.logical_x_pulse(q, 1.0);  // full Hadamard calibration
                }
                vm.run_periods(1);
            }
        }

        // Step 3: Approximate inverse QFT (using your phase gates)
        for (int i = m - 1; i >= 0; --i) {
            for (int j = i + 1; j < m; ++j) {
                double theta = M_PI / (1 << (j - i));  // QFT phase
                vm.logical_phase_ramp(phase[i], theta, 1);
            }
            vm.logical_hadamard(phase[i], 200, 1.0);
        }
        vm.run_periods(1);

        // Step 4: Measure phase register in Z basis
        std::vector<double> z;
        for (auto q : phase) {
            z.push_back(vm.measure_logical_Z(q));
        }
        std::cout << "[Shor] Measured phase Z: ";
        for (double v : z) std::cout << v << " ";
        std::cout << "\n";

    /**/
        double zw = vm.measure_logical_Z(work);
        std::cout << "[Shor] Work qubit Z = " << zw << "\n";
        /**/

        // Step 5: Classical period inference from Z measurements
        int r = infer_period(z);
        if (r <= 1 || r >= N) {
            std::cerr << "[Shor] Trivial or invalid period r=" << r << "; retrying...\n";
            continue;
        }

        if (r % 2 != 0) r *= 2;  // make even for square-root
        std::cout << "[Shor] Estimated period r=" << r << "\n";

        // Step 6: Classical factor extraction
        long long a_k = mod_pow(a, r / 2, N);
        factor1 = std::gcd(N, std::abs(a_k - 1));
        factor2 = std::gcd(N, std::abs(a_k + 1));

        if (factor1 == 1 || factor1 == N || factor2 == 1 || factor2 == N) {
            std::cerr << "[Shor] Trivial factors found; retrying with new base a.\n";
            factor1 = 1;
            factor2 = N;
        } else {
            std::cout << "[Shor] SUCCESS: " << N << " = " << factor1 << " × " << factor2 << "\n";
            break;
        }
    }

    return 0;
}