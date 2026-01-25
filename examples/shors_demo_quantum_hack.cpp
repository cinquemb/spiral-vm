// examples/shors_demo_spiral_shor.cpp
// SpiralVM-native Shor: a-dependent entangled period finding

#include "../src/spiral_vm_core.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <numeric>
#include <unordered_set>

double mean_abs(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    double s = 0.0;
    for (double x : v) s += std::abs(x);
    return s / v.size();
}

double variance(const std::vector<double>& v) {
    if (v.size() < 2) return 0.0;
    double mu = 0.0;
    for (double x : v) mu += x;
    mu /= v.size();

    double s2 = 0.0;
    for (double x : v) s2 += (x - mu) * (x - mu);
    return s2 / (v.size() - 1);
}


bool is_coprime(long long a, long long b) {
    return std::gcd(static_cast<unsigned long long>(a),
                     static_cast<unsigned long long>(b)) == 1;
}

uint64_t modmul(uint64_t a, uint64_t b, uint64_t mod) {
    __int128 res = ( __int128 )a * b;
    return (uint64_t)(res % mod);
}



long long mod_pow(long long base, long long exp, long long mod) {
    long long result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) result = modmul(result, base, mod);
        base = modmul(base, base, mod);
        exp >>= 1;
    }
    return result;
}

double phase_estimate(const std::vector<double>& z) {
    double acc = 0.0;
    for (size_t i = 0; i < z.size(); ++i) {
        acc += (z[i] > 0 ? 1.0 : 0.0) / (1ULL << (i + 1));
    }
    return acc;
}

long long continued_fraction_denominator(double x, long long max_den = 1e6) {
    long long h1 = 1, h2 = 0;
    long long k1 = 0, k2 = 1;
    double b = x;

    while (true) {
        long long a = static_cast<long long>(std::floor(b));
        long long h = a * h1 + h2;
        long long k = a * k1 + k2;

        if (k > max_den) return k1;

        if (std::abs(x - double(h)/k) < 1e-6)
            return k;

        h2 = h1; h1 = h;
        k2 = k1; k1 = k;
        b = 1.0 / (b - a);
    }
}

bool should_zoom(SpiralVM& vm, uint32_t q) {
    return vm.winding[q].crossings >= 3;
}


// Spiral analogue of "multiply by a mod N"
void apply_Ua(SpiralVM& vm, uint32_t work, long long ak, long long N, int m, double zoom) {
    // Direct continuous embedding of ak mod N as phase kick
    // ak is already a^{2^k} mod N from the caller
    double frac = double(ak % N) / double(N);           // [0,1) fraction
    double phase_resolution = 1.0 / (1LL << m);   // m = # phase qubits
    //double theta = 2.0 * M_PI * frac / phase_resolution; // full 2π range for one full cycle
    double theta = zoom * 2.0 * M_PI * frac / phase_resolution;


   // double theta = 2.0 * M_PI * frac;                   

    // Optional: scale theta slightly to prevent over-rotation on large ak
    // (tune factor 0.5–2.0 depending on your lattice stability)
   // theta *= 0.8;  // start conservative, increase if phase spread is too weak
    
    // keep it bounded
    theta = std::fmod(theta, 2.0 * M_PI);

    vm.logical_phase_ramp(work, theta, 1);              // one period ramp
    vm.run_periods(1);
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: ./shors_demo_quantum_retry <N> [lattice_size=8] [use_phi_direct=1]\n";
        std::cerr << "  N = number to factor (odd semiprime >1)\n";
        std::cerr << "  lattice_size = square lattice side (default 8)\n";
        std::cerr << "  use_phi_direct = 1 (fast phi mode) or 0 (waveform mode), default 1\n";
        return 1;
    }

    std::unordered_set<long long> seen;


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

    int m_min = std::ceil(std::log2(N));
    int m = m_min + 2;   // start modest

    std::vector<long long> best_as;
    double best_score = 1e18;
    double best_a = 0;

    while (factor1 == 1 || factor1 == N) {
        // ===== Reset per-a winding diagnostics =====
        std::vector<int> crossings_history;
        std::vector<double> delta_history;
        std::vector<double> accel_history;

        std::vector<double> z;
        long long a_k;
        double phi;
        long long r;

        double zoom =1.0;
        double score = 0.0;

        attempt++;
        std::cout << "\n[Shor] Attempt #" << attempt << " on N=" << N << "\n";

        SpiralVM vm(lattice_size, lattice_size);
        vm.is_ang = true;
        vm.auto_compile_enabled = false;
        vm.use_phi_direct = use_phi;  // <--- the toggle you wanted
        vm.initialize_state("neel");

        //int m = static_cast<int>(std::ceil(2 * std::log2(N)));

        


        std::vector<uint32_t> phase;
        for (int i = 0; i < m; ++i) {
            phase.push_back(vm.add_qubit(2 + (i % 4), 2 + (i / 4)));
        }

        uint32_t work = vm.add_qubit(6, 6);
        vm.compile_to_physical_waveform();

        long long a;
        do {
            if(!best_as.empty() && rand() % 3 == 0) {  // 33% guided
                a = best_as[rand() % best_as.size()];
                long long g = 2 + (rand() % 11);
                a = mod_pow(a, g, N);
                std::cout << "[Guide] Sampling from top-K pool\n";
            } else {
                a = 2 + (rand() % (N - 3));  // 67% explore
            }

            if (seen.count(a)) continue;
            if (std::gcd(a, N) != 1) continue;
            if (mod_pow(a, 2, N) == 1) continue;

            seen.insert(a);
            break;

        } while (true);



        std::cout << "[Shor] Random base a=" << a << "\n";

        // Step 1: Hadamard on phase register
        for (auto q : phase)
            vm.logical_hadamard(q, 200, 1.0);
        vm.run_periods(1);

        // Step 2: Entangled controlled-Ua^x analogue
        for (int k = 0; k < m; ++k) {
            static long long ak = 1;
            if (k == 0) ak = a % N;
            else ak = modmul(ak, ak, N);

            std::cout << "[Shor] k=" << k << "  a^(2^k) mod N=" << ak << "\n";

            // Apply U_a^{2^k} as single scaled phase ramp (no repetition needed)
            apply_Ua(vm, work, ak, N, m, zoom);

            int Ck = vm.winding[work].crossings;
            crossings_history.push_back(Ck);

            int n = crossings_history.size();

            // First difference: ΔC
            if (n >= 2) {
                double delta = crossings_history[n - 1] - crossings_history[n - 2];
                delta_history.push_back(delta);
            }

            // Second difference: Δ²C (phase acceleration)
            if (n >= 3) {
                double delta1 = crossings_history[n - 1] - crossings_history[n - 2];
                double delta2 = crossings_history[n - 2] - crossings_history[n - 3];
                double accel  = delta1 - delta2;
                accel_history.push_back(accel);
            }


            // ---- Early abort: chaotic winding ----
            if (accel_history.size() >= 4) {

                double mean_accel = mean_abs(accel_history);
                double var_delta  = variance(delta_history);

                // Tunable thresholds (start conservative)
                if (mean_accel > 1.5 || var_delta > 4.0) {
                    std::cout << "[Guide] Early abort: chaotic phase "
                              << "(|Δ²C|=" << mean_accel
                              << ", var(ΔC)=" << var_delta << ")\n";
                    goto NEXT_A;  // jump to next a
                }
            }


            // Imprint a scaled version of the exponent onto phase[k]
            // (use log scale or direct fraction — log2 gives smoother gradient)
            double log_frac = double(k) / double(m);
            double imprint = 2.0 * M_PI * log_frac;  // [0, 2π] ramp across register
            vm.logical_phase_ramp(phase[k], imprint, 1);
            vm.run_periods(1);
        }

        // ---- Score this a ----
        score =
            mean_abs(accel_history)
            + 0.25 * variance(delta_history);

        if (score < best_score) {
            best_score = score;
            best_a = a;
            best_as.push_back(a);
            std::sort(best_as.begin(), best_as.end(), [&](long long x, long long y) {
                return mod_pow(x, 2, N) != 1 && mod_pow(y, 2, N) == 1;  // heuristic pre-filter
            });
            if(best_as.size() > 5) best_as.resize(5);

            std::cout << "[Guide] New best a=" << a << " score=" << score 
              << " (top-" << best_as.size() << ")\n";
        }

        // Step 3: inverse QFT (approximate)
        for (int i = m - 1; i >= 0; --i) {
            for (int j = i + 1; j < m; ++j) {
                double theta = M_PI / (1 << (j - i));
                vm.logical_phase_ramp(phase[i], theta, 1);
            }
            vm.logical_hadamard(phase[i], 200, 1.0);
        }
        vm.run_periods(1);

        // Step 4: measure phase register
        for (auto q : phase)
            z.push_back(vm.measure_logical_Z(q));

        std::cout << "[Shor] Phase Z: ";
        for (double v : z) std::cout << v << " ";
        std::cout << "\n";

        phi = phase_estimate(z);
        r = continued_fraction_denominator(phi, N);

        if (r % 2 != 0) r *= 2;
        if (r <= 1 || r >= N) {
            continue;
        }
        if (r <= 1 || r >= N || mod_pow(a, r, N) != 1) {

            uint32_t dominant = work;
            bool s_zoom = false;
            for (uint32_t q = 0; q < vm.winding.size(); ++q) {
                s_zoom = should_zoom(vm, q);
                if (vm.winding[q].crossings > vm.winding[dominant].crossings)
                    dominant = q;
            }


            //TODO, loop over recent detected phase wraps for every qubit    
            double phi0 = vm.logical_phase[dominant];
            if(s_zoom) {
                std::cout << "[Zoom] Recentering and increasing phase gain\n";
                vm.recenter_phase(dominant, phi0);   
                zoom *= 2.0;   // exponential zoom       
                continue;      
            }
            
        }

        std::cout << "[Shor] Estimated period r=" << r << "\n";

        a_k = mod_pow(a, r / 2, N);
        factor1 = std::gcd(N, std::abs(a_k - 1));
        factor2 = std::gcd(N, std::abs(a_k + 1));

        if (factor1 == 1 || factor1 == N ||
            factor2 == 1 || factor2 == N) {
           // std::cerr << "[Shor] Trivial factors; retrying...\n";
            factor1 = 1;
            factor2 = N;
        } else {
            std::cout << "[Shor] SUCCESS: "
                      << N << " = " << factor1
                      << " × " << factor2 << "\n";
            break;
        }


        NEXT_A:
            continue;
    }

    return 0;
}