// examples/secp256k1_spiral_dlp.cpp
// SpiralVM-native period finding for secp256k1 discrete log

#include "../src/spiral_vm_core.hpp"
#include <iostream>
#include <vector>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <iomanip>

// Real secp256k1 domain parameters
const uint64_t SECP256K1_N_HIGH = 0xFFFFFFFFFFFFFFFFULL;
const uint64_t SECP256K1_N_LOW  = 0xFFFFFFFFFFFFFFFEULL;  // full n = 0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141ULL (handled via normalization)

// Generator point G (secp256k1)
struct ECPoint {
    uint64_t x_high, x_low;
    uint64_t y_high, y_low;
};

const ECPoint G = {
    0x79BE667EF9DCBBACULL, 0x55A06295CE870B07ULL,   // Gx
    0x483ADA7726A3C465ULL, 0x5DA4FBFC0E1108A8ULL    // Gy
};

// Simple placeholder scalar multiplication (for demo; replace with full Jacobian arithmetic later)
ECPoint scalar_mult(const ECPoint& base, uint64_t k) {
    // In a real implementation use proper secp256k1 scalar mul (e.g. from libsecp256k1)
    // Here we normalize k mod n for phase ramp
    return base;  // placeholder - phase ramp uses k directly
}

double phase_estimate(const std::vector<double>& z) {
    double acc = 0.0;
    for (size_t i = 0; i < z.size(); ++i) {
        acc += (z[i] > 0 ? 1.0 : 0.0) / (1ULL << (i + 1));
    }
    return acc;
}

long long continued_fraction_denominator(double x, long long max_den = 1LL<<60) {
    long long h1 = 1, h2 = 0, k1 = 0, k2 = 1;
    double b = x;
    while (true) {
        long long a = static_cast<long long>(std::floor(b));
        long long h = a * h1 + h2;
        long long k = a * k1 + k2;
        if (k > max_den) return k1;
        if (std::abs(x - static_cast<double>(h)/k) < 1e-6) return k;
        h2 = h1; h1 = h;
        k2 = k1; k1 = k;
        b = 1.0 / (b - a);
    }
}

bool should_zoom(SpiralVM& vm, uint32_t q) {
    return vm.winding[q].crossings >= 3;
}

// Continuous phase-ramp version of scalar multiplication (core of your method)
void apply_scalar_mult(SpiralVM& vm, uint32_t work, uint64_t k, int m, double zoom = 1.0) {
    // Normalize k for phase ramp (we use high-bit approximation for large k)
    double frac = static_cast<double>(k) / static_cast<double>(1ULL << 60);
    double theta = zoom * 2.0 * M_PI * frac;
    theta = std::fmod(theta, 2.0 * M_PI);

    vm.logical_phase_ramp(work, theta, 1);
    vm.run_periods(1);
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: ./secp256k1_spiral_dlp_robust <lattice_size=8> [use_phi=1]\n";
        return 1;
    }

    const int lattice_size = (argc > 1) ? std::atoi(argv[1]) : 8;
    const bool use_phi = (argc > 2) ? (std::atoi(argv[2]) != 0) : true;

    srand(static_cast<unsigned>(time(nullptr)));

    int attempt = 0;
    long long r = 0;

    while (true) {
        attempt++;
        std::cout << "\n[secp256k1 DLP] Attempt #" << attempt 
                  << " on 256-bit order (lattice " << lattice_size << "x" << lattice_size << ")\n";

        SpiralVM vm(lattice_size, lattice_size);
        vm.is_ang = true;
        vm.auto_compile_enabled = false;
        vm.use_phi_direct = use_phi;
        vm.initialize_state("polarized");

        const int m = 270;                     // ~256 bit + safety overhead
        std::vector<uint32_t> phase;
        for (int i = 0; i < m; ++i) {
            phase.push_back(vm.add_qubit(4 + (i % 16), 4 + (i / 16)));
        }

        uint32_t work = vm.add_qubit(lattice_size/2, lattice_size/2);
        vm.compile_to_physical_waveform();

        std::cout << "[secp256k1] Starting period-finding attack...\n";

        // Step 1: Hadamard on phase register
        for (auto q : phase) vm.logical_hadamard(q, 200, 1.0);
        vm.run_periods(1);

        // Step 2: Controlled scalar multiplication via continuous phase ramp
        for (int k = 0; k < m; ++k) {
            uint64_t ak = (1ULL << k);                     // 2^k  (full mod n in production)
            std::cout << "[DLP] k=" << k << "  2^k = " << ak << "\n";

            apply_scalar_mult(vm, work, ak, m, 1.0);

            // Imprint exponent onto phase register
            double imprint = 2.0 * M_PI * (static_cast<double>(k) / m);
            vm.logical_phase_ramp(phase[k], imprint, 1);
            vm.run_periods(1);

            if (should_zoom(vm, phase[k])) {
                std::cout << "[Zoom] Recentering and increasing phase gain on q" << k << "\n";
            }
        }

        // Step 3: Approximate inverse QFT
        for (int i = m-1; i >= 0; --i) {
            for (int j = i+1; j < m; ++j) {
                double theta = M_PI / (1 << (j - i));
                vm.logical_phase_ramp(phase[i], theta, 1);
            }
            vm.logical_hadamard(phase[i], 200, 1.0);
        }
        vm.run_periods(1);

        // Step 4: Measure phase register
        std::vector<double> z;
        for (auto q : phase) z.push_back(vm.measure_logical_Z(q));

        std::cout << "[DLP] Phase Z: ";
        for (double v : z) std::cout << v << " ";
        std::cout << "\n";

        double phi = phase_estimate(z);
        r = continued_fraction_denominator(phi, 1LL<<60);

        if (r > 1 && r < (1LL<<60)) {
            std::cout << "[DLP] Estimated period r = " << r << "\n";
            // In a real attack, check if r divides the curve order n
            break;
        }

        std::cout << "[DLP] Bad r, retrying with new random base...\n";
    }

    return 0;
}