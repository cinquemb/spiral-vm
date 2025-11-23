/*

Requirements: Armadillo installed (e.g., via libarmadillo-dev on Ubuntu) and linked with LAPACK/BLAS.
Compile: g++ -O2 spiral_vm_core.cpp -o spiral_vm -larmadillo `pkg-config lapack --libs` `pkg-config blas --libs`
Run: ./spiral_vm

sudo apt-get install libopenblas-openmp-dev libarmadillo-dev libblas-dev liblapack-dev gfortran

*/

// spiral_vm.cpp

#include <armadillo>
#include <complex>
#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <sstream>
#include <cmath>

#include "../src/spiral_vm_core.hpp"

using namespace arma;
using namespace std;

const int N = 900;
const int D = 2;
const int ROWS = 30;
const int COLS = 30;

struct LogicalQubit {
    uint32_t center_x;
    uint32_t center_y;
};

class SpiralVM {
public:
    static constexpr int D = 2;
    const int R = 1; // Physical neighborhood radius around each logical qubit's center
    const int rows, cols;         // Lattice dimensions
    const int N;                  // Number of sites
    double J, h0, h1, omega, T;  // Hamiltonian / Floquet parameters
    bool is_ang;                 // Spiral angle flag
    bool overlap_enabled = false; // Overlap mode toggle

private:
    cx_mat phi;                  // Quantum state vector (2*N x 1)
    cx_mat phi_in;               // Initial state for fidelity measurement
    int steps;                   // RK4 steps per period
    int current_period;          // Tracks Floquet periods elapsed
    double sx_gain;              // h1 gain parameter

    mt19937 rng;
    uniform_real_distribution<double> dist;

    vector<double> fidelities;
    vector<double> fidelity_window;
    std::vector<LogicalQubit> logical_qubits; // Logical qubit list

    // Mapping physical qubits to logical qubits (allows many logicals per physical)
    std::vector<std::vector<uint32_t>> phys_to_logicals;

public:
    SpiralVM(int r, int c) : rows(r), cols(c), N(r * c),
        J(0.3), h0(0), h1(2.5 / 4.4),
        omega(20 * 2 * datum::pi), T(2 * datum::pi / omega),
        is_ang(true), steps(200),
        sx_gain(1900.0), current_period(0),
        rng(777), dist(0.0, 1.0),
        fidelities(5001, 0.0), fidelity_window(32, 0.0),
        phys_to_logicals(N) // resize to number of physical qubits
    {
        phi = zeros<cx_mat>(N * D, 1);
        phi_in = zeros<cx_mat>(N * D, 1);
    }

    // Initialize quantum state: "neel", "polarized", or "disordered"
    void initialize_state(const string& initial_state = "neel") {
        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < cols; col++) {
                int i = row * cols + col;
                if (initial_state == "neel") {
                    if ((row + col) % 2 == 0) {
                        phi(i * D, 0) = 1.0;       // |0>
                        phi(i * D + 1, 0) = 0.0;
                    } else {
                        phi(i * D, 0) = 0.0;
                        phi(i * D + 1, 0) = 1.0;   // |1>
                    }
                } else if (initial_state == "polarized") {
                    phi(i * D, 0) = 1.0;           // |↑>
                    phi(i * D + 1, 0) = 0.0;
                } else if (initial_state == "disordered") {
                    if (dist(rng) < 0.5) {
                        phi(i * D, 0) = 1.0;       // |↑>
                        phi(i * D + 1, 0) = 0.0;
                    } else {
                        phi(i * D, 0) = 0.0;
                        phi(i * D + 1, 0) = 1.0;   // |↓>
                    }
                }
            }
        }
        double norm0 = sqrt(inner_product_cl10(phi, phi));
        phi /= norm0;
        phi_in = phi;
        fidelities[0] = 1.0;

        current_period = 0;

        cout << "Initial norm: " << norm0 << "\n";
    }

    void run_floquet(int N_max, const string& initial_state) {
        ofstream fout;
        stringstream fname;
        fname << "dtc_floquet_with_" << initial_state << "_state_" << (is_ang ? "_spiral_" : "_no_spiral_") << (int)sx_gain << "_omega_" << omega << "_T_" << T << ".txt";
        fout.open(fname.str());
        fout << "Period,Fidelity,Stabilizer,Energy,sx_energy,zz_energy,Delta_F,ht_eff_end,sx_avg\n";

        double delta_F = 0.0;

        sp_cx_mat H_init = hamiltonian_cl10_90_spiral_twist(J, 0, 0);
        cx_mat Hphi_init = mat_vec_mult_cl10(H_init, phi);
        double energy_init = real(inner_product_cl10(phi, Hphi_init)) + compute_zz_energy(phi, J, 0, 0);
        double zz_energy_init = compute_zz_energy(phi, J, 0, 0);
        double sx_energy_init = 0.0;
        for (int i = 0; i < N; i++) {
            sx_energy_init -= h1 * 2.0 * real(phi(i * D, 0) * conj(phi(i * D + 1, 0)));
        }
        fout << "0,1.0," << compute_avg_stabilizer(phi) << "," << energy_init << "," << sx_energy_init << "," << zz_energy_init << "," << h1 << ",0,0\n";

        for (int n = 0; n < N_max; n++) {
            step_period(n, delta_F);

            fout << n + 1 << "," << fidelities[n + 1] << "," << compute_avg_stabilizer(phi) << ",";
            double sx_energy = 0.0;
            for (int i = 0; i < N; i++) {
                sx_energy -= phi(i * D, 0) * conj(phi(i * D + 1, 0)) * 2.0 * h1;
            }
            double zz_energy = is_ang ? compute_zz_energy(phi, J, omega_ang_end(n), n + 1, true) : compute_zz_energy(phi, J, 0, 0);
            double energy = sx_energy + zz_energy;

            fout << energy << "," << sx_energy << "," << zz_energy << "," << delta_F << "," << h_effective_end(n) << "," << sx_avg(n) << "\n";

            cout << "Period " << n + 1 << ": Fidelity = " << fidelities[n + 1] << ", Energy = " << energy << "\n";

            if (n % 2 == 1 && fidelities[n + 1] > 0.9) {
                cout << "Period-doubling at " << n + 1 << "T\n";
            }
        }
        fout.close();
    }

    void step_period(int n, double &delta_F) {
        double dt = T / steps;
        double delta_F_target = 1.0;

        cx_mat phi_new = phi;

        double k_bA = 0.0;
        double h1_limit = 10000.0;
        double h1_f_gain = 0.0;
        double k_oA = 1.0;
        double alpha_ang = is_ang ? 1.0 : 0.0;
        double beta_ang = 0.0;
        double k_angular = 0.0;
        double omega_omega_angA = omega;
        double omega_omega_angB = 2 * omega_omega_angA;
        double omega_ang_base = 0.0;

        for (int k = 0; k < steps; k++) {
            double t = n * T + k * dt;

            double angular_freq_quasi = alpha_ang * (sin(omega_omega_angA * M_PI * t / T) + sin(omega_omega_angB * M_PI * t / T));
            double angular_freq_feedback = beta_ang * (delta_F_target - delta_F);
            double omega_ang_mod = omega_ang_base + k_oA * sin(omega_omega_angA * M_PI * t / T);
            double omega_ang = omega_ang_mod + angular_freq_quasi + k_angular * angular_freq_feedback;

            double h1_feedback = h1_f_gain * (k_bA * (delta_F_target - delta_F));
            double ht_eff = h0 + h1 * cos(omega * (t / 2) + M_PI / 4.0) + h1_feedback;

            if (ht_eff / J > h1_limit) ht_eff = h1_limit * J;
            else if (ht_eff / J < -h1_limit) ht_eff = -h1_limit * J;

            sp_cx_mat H_sx = hamiltonian_cl10_90_spiral_twist(J, ht_eff, omega_ang);

            if (!is_ang) {
                cx_mat Hzz = compute_zz_energy_vector(phi_new, J, omega_ang, 0);
                cx_mat k1 = mat_vec_mult_cl10(H_sx, phi_new) + Hzz;
                cx_mat phi_temp = phi_new + (-cx_double(0, 1) * dt / 2.0) * k1;
                cx_mat k2 = mat_vec_mult_cl10(H_sx, phi_temp) + compute_zz_energy_vector(phi_temp, J, 0, 0);
                phi_temp = phi_new + (-cx_double(0, 1) * dt / 2.0) * k2;
                cx_mat k3 = mat_vec_mult_cl10(H_sx, phi_temp) + compute_zz_energy_vector(phi_temp, J, 0, 0);
                phi_temp = phi_new + (-cx_double(0, 1) * dt) * k3;
                cx_mat k4 = mat_vec_mult_cl10(H_sx, phi_temp) + compute_zz_energy_vector(phi_temp, J, 0, 0);
                phi_new += (-cx_double(0, 1) * dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
                phi_new /= sqrt(real(inner_product_cl10(phi_new, phi_new)));
            } else {
                cx_mat Hzz = compute_zz_energy_vector(phi_new, J, omega_ang, t, true);
                cx_mat k1 = mat_vec_mult_cl10(H_sx, phi_new) + Hzz;
                cx_mat phi_temp = phi_new + (-cx_double(0, 1) * dt / 2.0) * k1;
                cx_mat k2 = mat_vec_mult_cl10(H_sx, phi_temp) + compute_zz_energy_vector(phi_temp, J, omega_ang, t, true);
                phi_temp = phi_new + (-cx_double(0, 1) * dt / 2.0) * k2;
                cx_mat k3 = mat_vec_mult_cl10(H_sx, phi_temp) + compute_zz_energy_vector(phi_temp, J, omega_ang, t, true);
                phi_temp = phi_new + (-cx_double(0, 1) * dt) * k3;
                cx_mat k4 = mat_vec_mult_cl10(H_sx, phi_temp) + compute_zz_energy_vector(phi_temp, J, omega_ang, t, true);
                phi_new += (-cx_double(0, 1) * dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
                double normrk4 = sqrt(real(inner_product_cl10(phi_new, phi_new)));
                phi_new /= normrk4;
            }
        }
        phi = phi_new;

        fidelities[n + 1] = abs(inner_product_cl10(phi_in, phi));
        delta_F = 1.0 - fidelities[n + 1];
        current_period++;
    }

    double omega_ang_end(int n) {
        double t_end = n * T;
        double base = omega_ang_base;
        double omega_omega_angA = omega;
        double omega_omega_angB = 2 * omega;
        double alpha_ang = is_ang ? 1.0 : 0.0;
        double quasi = alpha_ang * (sin(omega_omega_angA * M_PI * t_end / T) +
                                    sin(omega_omega_angB * M_PI * t_end / T));
        return base + quasi;
    }

    double h_effective_end(int n) {
        double t_end = n * T;
        double h1_osc = h1 * cos(omega * (t_end / 2.0) + M_PI / 4.0);
        double h_eff = h0 + h1_osc;
        double h1_limit = 10000.0;

        if (h_eff / J > h1_limit) h_eff = h1_limit * J;
        else if (h_eff / J < -h1_limit) h_eff = -h1_limit * J;

        return h_eff;
    }

    double sx_avg(int n) {
        double sx_sum = 0.0;
        for (int i = 0; i < N; i++) {
            sx_sum += 2.0 * real(phi(i * D, 0) * conj(phi(i * D + 1, 0)));
        }
        return sx_sum / N;
    }

    //------------------------------------------------------------------------
    // Add logical qubit and update phys->logical mapping if overlap enabled
    uint32_t add_qubit(uint32_t x, uint32_t y) {
        uint32_t qid = logical_qubits.size();
        logical_qubits.push_back({x, y});

        if (overlap_enabled) {
            for (int row = (int)y - R; row <= (int)y + R; ++row) {
                if (row < 0 || row >= rows) continue;
                for (int col = (int)x - R; col <= (int)x + R; ++col) {
                    if (col < 0 || col >= cols) continue;
                    int phys_idx = row * cols + col;
                    phys_to_logicals[phys_idx].push_back(qid);
                }
            }
        }

        return qid;
    }

    //------------------------------------------------------------------------
    void run_periods(uint32_t N_periods) {
        double dummy_deltaF = 0.0;
        for (uint32_t i = 0; i < N_periods; i++) {
            step_period(current_period, dummy_deltaF);
            current_period++;
        }
    }

    void apply_global_pi_pulse_on_even_cycles() {
        if (current_period % 2 == 0) {
            std::cout << "Applying logical X (global pi pulse) at period " << current_period << "\n";
            global_pi_pulse();
        } else {
            std::cout << "Skipping logical X on odd period " << current_period << "\n";
        }
    }

    void global_pi_pulse() {
        for (int i = 0; i < N; i++) {
            cx_double temp = phi(i * D, 0);
            phi(i * D, 0) = phi(i * D + 1, 0);
            phi(i * D + 1, 0) = temp;
        }
    }

    double measure_even_population(uint32_t qid) {
        if (qid >= logical_qubits.size()) {
            std::cerr << "Invalid logical qubit ID\n";
            return 0.0;
        }
        const LogicalQubit& q = logical_qubits[qid];

        double pop_sum = 0.0;
        int count = 0;

        for (int row = 0; row < rows; ++row) {
            for (int col = 0; col < cols; ++col) {
                int dist = abs((int)row - (int)q.center_y) + abs((int)col - (int)q.center_x);
                if (dist <= R && ((row + col) % 2 == 0)) {
                    int i = row * cols + col;
                    pop_sum += std::norm(phi(i * D, 0));
                    count++;
                }
            }
        }

        if (count == 0) return 0.0;
        return pop_sum / count;
    }

    void apply_phase_shift(double angle) {
        cx_double phase_factor = std::exp(cx_double(0, angle));
        for (int i = 0; i < N * D; i++) {
            phi(i, 0) *= phase_factor;
        }
    }

    void apply_phase_kick_between(uint32_t qid1, uint32_t qid2, double strength, double duration_fraction) {
        if (qid1 >= logical_qubits.size() || qid2 >= logical_qubits.size()) {
            std::cerr << "Invalid logical qubit IDs for phase kick.\n";
            return;
        }
        const LogicalQubit& q1 = logical_qubits[qid1];
        const LogicalQubit& q2 = logical_qubits[qid2];

        // Identify physical qubits in neighborhoods
        for (int row1 = q1.center_y - R; row1 <= q1.center_y + R; ++row1) {
            if (row1 < 0 || row1 >= rows) continue;
            for (int col1 = q1.center_x - R; col1 <= q1.center_x + R; ++col1) {
                if (col1 < 0 || col1 >= cols) continue;

                for (int row2 = q2.center_y - R; row2 <= q2.center_y + R; ++row2) {
                    if (row2 < 0 || row2 >= rows) continue;
                    for (int col2 = q2.center_x - R; col2 <= q2.center_x + R; ++col2) {
                        if (col2 < 0 || col2 >= cols) continue;

                        int i1 = row1 * cols + col1;
                        int i2 = row2 * cols + col2;

                        cx_double phase = std::exp(cx_double(0, strength * duration_fraction));

                        phi(i1 * D + 1, 0) *= phase;
                        phi(i2 * D + 1, 0) *= phase;
                    }
                }
            }
        }

        // Normalize state vector
        double norm = sqrt(real(inner_product_cl10(phi, phi)));
        phi /= norm;
    }

    double logical_zz_correlation(uint32_t qid1, uint32_t qid2) {
        if (qid1 >= logical_qubits.size() || qid2 >= logical_qubits.size()) {
            std::cerr << "Invalid logical qubit IDs for ZZ correlation measurement.\n";
            return 0.0;
        }
        const LogicalQubit& q1 = logical_qubits[qid1];
        const LogicalQubit& q2 = logical_qubits[qid2];

        double sum_correlation = 0.0;
        int count = 0;

        for (int row1 = q1.center_y - R; row1 <= q1.center_y + R; ++row1) {
            if (row1 < 0 || row1 >= rows) continue;
            for (int col1 = q1.center_x - R; col1 <= q1.center_x + R; ++col1) {
                if (col1 < 0 || col1 >= cols) continue;

                for (int row2 = q2.center_y - R; row2 <= q2.center_y + R; ++row2) {
                    if (row2 < 0 || row2 >= rows) continue;
                    for (int col2 = q2.center_x - R; col2 <= q2.center_x + R; ++col2) {
                        if (col2 < 0 || col2 >= cols) continue;

                        int i1 = row1 * cols + col1;
                        int i2 = row2 * cols + col2;

                        double z_expect_i1 = std::norm(phi(i1 * D, 0)) - std::norm(phi(i1 * D + 1, 0));
                        double z_expect_i2 = std::norm(phi(i2 * D, 0)) - std::norm(phi(i2 * D + 1, 0));
                        sum_correlation += z_expect_i1 * z_expect_i2;
                        count++;
                    }
                }
            }
        }
        if (count == 0) return 0.0;

        return sum_correlation / count;
    }

    double get_logical_phase(uint32_t qid) {
        if (qid >= logical_qubits.size()) {
            std::cerr << "Invalid logical qubit ID for get_logical_phase.\n";
            return 0.0;
        }
        const LogicalQubit& q = logical_qubits[qid];

        std::complex<double> phase_sum(0.0, 0.0);
        int count = 0;

        for (int row = q.center_y - R; row <= q.center_y + R; ++row) {
            if (row < 0 || row >= rows) continue;
            for (int col = q.center_x - R; col <= q.center_x + R; ++col) {
                if (col < 0 || col >= cols) continue;
                int i = row * cols + col;
                std::complex<double> amp1 = phi(i * D + 1, 0);
                if (std::abs(amp1) > 1e-12) {
                    std::complex<double> relative_phase = amp1 / std::abs(amp1);
                    phase_sum += relative_phase;
                    count++;
                }
            }
        }
        if (count == 0) return 0.0;

        std::complex<double> avg_phase = phase_sum / static_cast<double>(count);
        return std::arg(avg_phase);
    }

    void ramp_omega_ang(double start, double end, double duration_seconds) {
        double total_periods = duration_seconds / T;
        if (total_periods < 1) total_periods = 1;

        double omega_ang_delta = end - start;
        double omega_ang_step = omega_ang_delta / total_periods;

        for (int step = 0; step < static_cast<int>(total_periods); ++step) {
            omega_ang_base = start + omega_ang_step * step;
            double dummy_deltaF = 0.0;
            step_period(current_period, dummy_deltaF);
            current_period++;
        }
        omega_ang_base = end;
    }

    void print_overlap_stats() {
        if (!overlap_enabled) {
            std::cout << "Overlap mode is disabled.\n";
            return;
        }
        for (int i = 0; i < N; i++) {
            int count = phys_to_logicals[i].size();
            if (count > 1) {
                std::cout << "Physical qubit " << i << " belongs to " << count << " logical qubits.\n";
            }
        }
    }

private:
    cx_mat mat_vec_mult_cl10(const sp_cx_mat& H, const cx_mat& phi) {
        return H * phi;
    }

    double inner_product_cl10(const cx_mat& phi1, const cx_mat& phi2) {
        cx_double prod = as_scalar(phi1.t() * phi2);
        return real(prod);
    }

    double compute_zz_energy(const cx_mat& phi, double J, double omega_ang, double period, bool is_ang = false) {
        double energy = 0.0;
        double theta_max_base = M_PI / 512;
        double theta_max = theta_max_base * (1.0 - cos(M_PI * period / 2.0)) / 2.0;

        double center_x = cols / 2.0, center_y = rows / 2.0;
        double r_max = sqrt(center_x * center_x + center_y * center_y);
        double norm_factor = sqrt(rows * cols);

        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < cols; col++) {
                int i = row * cols + col;
                int j_right = row * cols + ((col + 1) % cols);
                int j_down = ((row + 1) % rows) * cols + col;

                double r = sqrt(pow(row - center_y, 2) + pow(col - center_x, 2));
                double phi_ang = atan2(row - center_y, col - center_x);
                double theta = (omega_ang == 0) ? 0 : theta_max * (r / r_max) * cos(omega_ang * phi_ang);

                cx_double J_twist;
                if (is_ang) {
                    J_twist = J * cx_double(0, sin(theta));
                    energy += imag(J_twist * (phi(i * D, 0) - phi(i * D + 1, 0)) * norm_factor * (phi(j_right * D, 0) - phi(j_right * D + 1, 0)) * norm_factor);
                    energy += imag(J_twist * (phi(i * D, 0) - phi(i * D + 1, 0)) * norm_factor * (phi(j_down * D, 0) - phi(j_down * D + 1, 0)) * norm_factor);
                } else {
                    J_twist = J;
                    energy += real(J_twist * (phi(i * D, 0) - phi(i * D + 1, 0)) * norm_factor * (phi(j_right * D, 0) - phi(j_right * D + 1, 0)) * norm_factor);
                    energy += real(J_twist * (phi(i * D, 0) - phi(i * D + 1, 0)) * norm_factor * (phi(j_down * D, 0) - phi(j_down * D + 1, 0)) * norm_factor);
                }
            }
        }
        return energy;
    }

    cx_mat compute_zz_energy_vector(const cx_mat& phi, double J, double omega_ang, double period, bool is_ang = false) {
        cx_mat Hzz_phi = zeros<cx_mat>(N * D, 1);

        double theta_max_base = M_PI / 512;
        double theta_max = theta_max_base * (1.0 - cos(M_PI * period / 2.0)) / 2.0;

        double center_x = cols / 2.0, center_y = rows / 2.0;
        double r_max = sqrt(center_x * center_x + center_y * center_y);
        double norm_factor = sqrt(N);

        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < cols; col++) {
                int i = row * cols + col;
                int j_right = row * cols + ((col + 1) % cols);
                int j_down = ((row + 1) % rows) * cols + col;

                double r = sqrt(pow(row - center_y, 2) + pow(col - center_x, 2));
                double phi_ang = atan2(row - center_y, col - center_x);
                double theta = (omega_ang == 0) ? 0 : theta_max * (r / r_max) * cos(omega_ang * phi_ang);

                cx_double J_twist = is_ang ? J * cx_double(0, sin(theta)) : J;

                cx_double z_jr = (phi(j_right * D, 0) - phi(j_right * D + 1, 0)) * norm_factor;
                cx_double z_jd = (phi(j_down * D, 0) - phi(j_down * D + 1, 0)) * norm_factor;

                Hzz_phi(i * D, 0) = J_twist * (z_jr + z_jd) * phi(i * D, 0);
                Hzz_phi(i * D + 1, 0) = J_twist * (-z_jr - z_jd) * phi(i * D + 1, 0);
            }
        }
        return Hzz_phi;
    }

    void compute_nonzero_indices_spiral_twist(double J, double ht, int rows, int cols, int D, double omega_ang, umat& locations, cx_vec& values, uint& nz) {
        int N = rows * cols;
        const int NNZ = 2 * N; // 2N for sigma^x

        nz = 0;

        for (int i = 0; i < N; i++) {
            locations(0, nz) = i * D;     locations(1, nz) = i * D + 1; values(nz) = -ht; nz++;
            locations(0, nz) = i * D + 1; locations(1, nz) = i * D;     values(nz) = -ht; nz++;
        }
    }

    sp_cx_mat hamiltonian_cl10_90_spiral_twist(double J, double ht, double omega_ang) {
        int N_local = rows * cols;
        const int NNZ = 2 * N_local;
        umat locations(2, NNZ);
        cx_vec values(NNZ);
        uint nz;
        compute_nonzero_indices_spiral_twist(J, ht, rows, cols, D, omega_ang, locations, values, nz);
        return sp_cx_mat(locations.submat(0, 0, 1, nz - 1), values.subvec(0, nz - 1), N_local * D, N_local * D);
    }

    double compute_avg_stabilizer(const cx_mat& phi) {
        double stab_sum = 0.0;
        int count = 0;

        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < cols; col++) {
                int i = row * cols + col;
                int j_right = row * cols + ((col + 1) % cols);
                int j_down = ((row + 1) % rows) * cols + col;

                double zi = abs(phi(i * D, 0) - phi(i * D + 1, 0));
                double zj_right = abs(phi(j_right * D, 0) - phi(j_right * D + 1, 0));
                double zj_down = abs(phi(j_down * D, 0) - phi(j_down * D + 1, 0));

                stab_sum += zi * zj_right + zi * zj_down;
                count += 2;
            }
        }
        return stab_sum / count;
    }
};

int main() {
    SpiralVM vm(30, 30);
    vm.initialize_state("neel");
    vm.run_floquet(5000, "neel");
    return 0;
}
