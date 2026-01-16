// Multi-resonance wrap-regime probe with twin-trajectory divergence
// Adaptive carpet-bomber: switches transverse axis based on λ_loc

#include "../src/spiral_vm_core.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

static constexpr double TWO_PI = 6.283185307179586;

// Wrap-seam resonances (pre- and post-wrap)
static constexpr double RESONANCES[] = {5.98, 6.01, 6.03, 0.02, 0.04};
static constexpr int    N_RES = sizeof(RESONANCES) / sizeof(double);

static constexpr double WINDOW  = 0.04;   // width per resonance
static constexpr double EPS     = 1e-4;   // initial separation
static constexpr double H_AMP   = 0.2;
static constexpr double Z_SLOPE = 0.03;
static constexpr double BOMB    = 0.03;

static inline double wrap_dist(double a, double b) {
    double d = std::fabs(a - b);
    return std::min(d, TWO_PI - d);
}

int main() {
    SpiralVM vm(30, 30);
    vm.initialize_state("neel");

    uint32_t q0 = vm.add_qubit(15, 15);

    // Prepare |0>_L
    vm.logical_x_pulse(q0, 1);
    vm.logical_x_pulse(q0, 1);

    std::cout << "step |   φ    |   Z    |   X    |  |δ|   |  λ_loc\n";
    std::cout << "----------------------------------------------------\n";

    double last_lambda = 0.0;
    int step = 0;

    while (step < 5000) {
        double phi = std::fmod(vm.logical_phase[q0], TWO_PI);
        if (phi < 0) phi += TWO_PI;

        // Composite window weight over all resonances
        double w = 0.0;
        bool in_window = false;
        for (int i = 0; i < N_RES; ++i) {
            double d = wrap_dist(phi, RESONANCES[i]);
            if (d < WINDOW) in_window = true;
            w += std::exp(-(d * d) / (WINDOW * WINDOW));
        }

        // If not near any seam, just advance
        if (!in_window) {
            vm.logical_hadamard_step(q0, H_AMP);
            vm.logical_phase_ramp(q0, Z_SLOPE, 1);
            step++;
            continue;
        }

        // Twin systems
        SpiralVM a = vm;
        SpiralVM b = vm;

        // Tiny transverse offset
        b.logical_phase_ramp(q0, EPS, 1);

        // Identical base evolution
        a.logical_hadamard_step(q0, H_AMP);
        b.logical_hadamard_step(q0, H_AMP);

        a.logical_phase_ramp(q0, Z_SLOPE, 1);
        b.logical_phase_ramp(q0, Z_SLOPE, 1);

        // Choose transverse axis from last λ
        bool expand = (step > 0 && last_lambda > 0.0);

        if (expand) {
            a.logical_x_pulse(q0, BOMB * w);
            b.logical_x_pulse(q0, BOMB * w);
        } else {
            a.logical_y_pulse(q0, BOMB * w);
            b.logical_y_pulse(q0, BOMB * w);
        }

        // Measure divergence
        double Za = a.measure_logical_Z(q0);
        double Xa = a.measure_logical_X(q0);
        double Zb = b.measure_logical_Z(q0);
        double Xb = b.measure_logical_X(q0);

        double d0 = EPS;
        double d1 = std::hypot(Zb - Za, Xb - Xa);
        double lambda = std::log(d1 / d0);
        last_lambda = lambda;

        std::cout << std::setw(4) << step << " | "
                  << std::setw(6) << phi << " | "
                  << std::setw(6) << Za << " | "
                  << std::setw(6) << Xa << " | "
                  << std::setw(7) << d1 << " | "
                  << std::setw(8) << lambda << "\n";

        // Evolve the *real* system the same way
        vm.logical_hadamard_step(q0, H_AMP);
        vm.logical_phase_ramp(q0, Z_SLOPE, 1);

        if (expand)
            vm.logical_x_pulse(q0, BOMB * w);
        else
            vm.logical_y_pulse(q0, BOMB * w);

        step++;
    }

    return 0;
}
