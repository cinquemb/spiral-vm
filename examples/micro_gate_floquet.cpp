// Floquet seam macro-gate probe
// Chaos-mediated bifurcation + sheet-commit + recombination
// Tests whether a nontrivial (X,Z)->(X',Z') map can emerge

#include "../src/spiral_vm_core.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

static constexpr double TWO_PI = 6.283185307179586;

// Wrap-seam resonances
static constexpr double RESONANCES[] = {5.98, 6.01, 6.03, 0.02, 0.04};
static constexpr int    N_RES = sizeof(RESONANCES) / sizeof(double);

static constexpr double WINDOW  = 0.04;
static constexpr double EPS     = 1e-4;
static constexpr double Z_SLOPE = 0.03;

static constexpr double H_AMP   = 0.35;
static constexpr double BOMB    = 0.08;
static constexpr int    RELAX_STEPS = 160;

static inline double wrap_dist(double a, double b) {
    double d = std::fabs(a - b);
    return std::min(d, TWO_PI - d);
}

// One chaos-mediated macro step
void macro_step(SpiralVM& vm, uint32_t q0, double& outX, double& outZ) {
    double last_lambda = 0.0;

    while (true) {
        double phi = std::fmod(vm.logical_phase[q0], TWO_PI);
        if (phi < 0) phi += TWO_PI;

        double w = 0.0;
        bool in_window = false;
        for (int i = 0; i < N_RES; ++i) {
            double d = wrap_dist(phi, RESONANCES[i]);
            if (d < WINDOW) in_window = true;
            w += std::exp(-(d * d) / (WINDOW * WINDOW));
        }

        if (!in_window) {
            vm.logical_hadamard_step(q0, H_AMP);
            vm.logical_phase_ramp(q0, Z_SLOPE, 1);
            continue;
        }

        // Twin systems for local Lyapunov probe
        SpiralVM a(vm);
        SpiralVM b(vm);
        b.logical_phase_ramp(q0, EPS, 1);

        a.logical_hadamard_step(q0, H_AMP);
        b.logical_hadamard_step(q0, H_AMP);
        a.logical_phase_ramp(q0, Z_SLOPE, 1);
        b.logical_phase_ramp(q0, Z_SLOPE, 1);

        bool expand = (last_lambda > 0.0);
        if (expand) {
            a.logical_x_pulse(q0,  BOMB * w);
            b.logical_x_pulse(q0, -BOMB * w);
        } else {
            a.logical_y_pulse(q0,  BOMB * w);
            b.logical_y_pulse(q0, -BOMB * w);
        }

        double Za = a.measure_logical_Z(q0);
        double Xa = a.measure_logical_X(q0);
        double Zb = b.measure_logical_Z(q0);
        double Xb = b.measure_logical_X(q0);

        double d1 = std::hypot(Zb - Za, Xb - Xa);
        last_lambda = std::log(d1 / EPS);

        // Apply same evolution to real system
        vm.logical_hadamard_step(q0, H_AMP);
        vm.logical_phase_ramp(q0, Z_SLOPE, 1);
        if (expand)
            vm.logical_x_pulse(q0,  BOMB * w);
        else
            vm.logical_y_pulse(q0, -BOMB * w);

        // --- Sheet-commit step ---
        double phi_post = std::fmod(vm.logical_phase[q0], TWO_PI);
        if (phi_post < 0) phi_post += TWO_PI;

        if (phi_post < M_PI) {
            vm.logical_phase_ramp(q0, +0.25, 1);
        } else {
            vm.logical_phase_ramp(q0, -0.25, 1);
        }

        break;
    }

    // Recombination / relaxation
    for (int i = 0; i < RELAX_STEPS; ++i) {
        vm.logical_hadamard_step(q0, H_AMP);
        vm.logical_phase_ramp(q0, Z_SLOPE, 1);
    }

    outX = vm.measure_logical_X(q0);
    outZ = vm.measure_logical_Z(q0);
}

int main() {
    std::cout << "X_in   Z_in   ->   X_out   Z_out\n";
    std::cout << "--------------------------------\n";

    for (int mode = 0; mode < 4; ++mode) {
        SpiralVM vm(30, 30);
        vm.is_ang = true;
        vm.initialize_state("neel");
        uint32_t q0 = vm.add_qubit(15, 15);

        if (mode == 1) {
            vm.logical_x_pulse(q0, M_PI);          // |1>
        } else if (mode == 2) {
            vm.logical_hadamard_step(q0, 0.6);    // |+>
        } else if (mode == 3) {
            vm.logical_hadamard_step(q0, 0.6);    // |->
            vm.logical_z_rotation(q0, M_PI);
        }

        double Xin = vm.measure_logical_X(q0);
        double Zin = vm.measure_logical_Z(q0);

        double Xout, Zout;
        macro_step(vm, q0, Xout, Zout);

        std::cout << std::setw(7) << Xin << " "
                  << std::setw(7) << Zin
                  << "  ->  "
                  << std::setw(7) << Xout << " "
                  << std::setw(7) << Zout << "\n";
    }

    return 0;
}
