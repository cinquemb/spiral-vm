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

using namespace arma;
using namespace std;

class SpiralVM {
public:
    static constexpr int D = 2;
    const int rows, cols;         // Lattice dimensions
    const int N;                  // Number of sites
    double J, h0, h1, omega, T;  // Hamiltonian / Floquet parameters
    bool is_ang;                 // Spiral angle flag

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
    std::vector<LogicalQubit> logical_qubits; // Add this as a member field

public:
    SpiralVM(int r, int c) : rows(r), cols(c), N(r * c),
        J(0.3), h0(0), h1(2.5 / 4.4),
        omega(20 * 2 * datum::pi), T(2 * datum::pi / omega),
        is_ang(true), steps(200),
        sx_gain(1900.0), current_period(0),
        rng(777), dist(0.0, 1.0),
        fidelities(5001, 0.0), fidelity_window(32, 0.0)
    {
        // Initialize state vectors to zero size 2*N x 1
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

    // Full Floquet evolution running N_max periods, saved to filename
    void run_floquet(int N_max, const string& initial_state) {
        ofstream fout;
        stringstream fname;
        fname << "dtc_floquet_with_" << initial_state << "_state_" << (is_ang ? "_spiral_" : "_no_spiral_") << (int)sx_gain << "_omega_" << omega << "_T_" << T << ".txt";
        fout.open(fname.str());
        fout << "Period,Fidelity,Stabilizer,Energy,sx_energy,zz_energy,Delta_F,ht_eff_end,sx_avg\n";

        double delta_F = 0.0;

        // Initial energy and stabilizer calculations
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

            // Logging fidelity, stabilizers, energy etc.
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

    // RK4 integration over one Floquet period n updating state phi and fidelity delta_F
    void step_period(int n, double &delta_F) {
        double dt = T / steps;
        double delta_F_target = 1.0;

        cx_mat phi_new = phi;

        // Coefficients for feedback, modulation, etc.
        double k_bA = 0.0; // feedback gain parameter (adjust as needed)
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

        // Update fidelity and delta_F
        fidelities[n + 1] = abs(inner_product_cl10(phi_in, phi));
        delta_F = 1.0 - fidelities[n + 1];
        current_period++;
    }

    double omega_ang_end(int n) {
        // Approximate omega_ang modulation at period end n
        // Adjust with actual formula if desired
        return 0.0;
    }

    double h_effective_end(int n) {
        // Approximate effective h1 at period end n
        return 0.0;
    }

    double sx_avg(int n) {
        // Compute sx average at period n (rough estimation)
        double sx_sum = 0.0;
        for (int i = 0; i < N; i++) {
            sx_sum += 2.0 * real(phi(i * D, 0) * conj(phi(i * D + 1, 0)));
        }
        return sx_sum / N;
    }

    //----------------------------------------------------------------
    // Add a logical qubit at (x,y) on lattice, return its ID
    uint32_t add_qubit(uint32_t x, uint32_t y) {
        logical_qubits.push_back({x, y});
        return logical_qubits.size() - 1;
    }

    //----------------------------------------------------------------
    // Run N Floquet periods, updating the quantum state.
    // This calls step_period() N times internally.
    void run_periods(uint32_t N_periods) {
        double dummy_deltaF = 0.0;
        for (uint32_t i = 0; i < N_periods; i++) {
            step_period(current_period, dummy_deltaF);
            current_period++;
        }
    }

    //----------------------------------------------------------------
    // Apply logical X gate: global pi pulse on even periods only.
    void apply_global_pi_pulse_on_even_cycles() {
        if (current_period % 2 == 0) {
            std::cout << "Applying logical X (global pi pulse) at period " << current_period << "\n";
            global_pi_pulse();
        } else {
            std::cout << "Skipping logical X on odd period " << current_period << "\n";
        }
    }

    // Helper that flips all qubits' sigma_x basis amplitudes (pi rotation around X)
    // This implements global pi pulse (logical X) on full lattice
    void global_pi_pulse() {
        for (int i = 0; i < N; i++) {
            // Swap amplitudes of |0> and |1> (indices i*D, i*D+1)
            cx_double temp = phi(i * D, 0);
            phi(i * D, 0) = phi(i * D + 1, 0);
            phi(i * D + 1, 0) = temp;
        }
    }

    //----------------------------------------------------------------
    // Measure logical qubit even cycle population (|0⟩_L) for qubit with ID qid
    // Here we define it as the population sum over physical qubits on one sublattice (even or odd)
    // For demonstration, assume logical qubit center vicinity defines the measurement region
    double measure_even_population(uint32_t qid) {
        if (qid >= logical_qubits.size()) {
            std::cerr << "Invalid logical qubit ID\n";
            return 0.0;
        }
        const LogicalQubit& q = logical_qubits[qid];

        // Define a measurement radius: count qubits within Manhattan radius R of center
        const int R = 3; // example radius, can adjust as needed

        double pop_sum = 0.0;
        int count = 0;

        for (int row = 0; row < rows; ++row) {
            for (int col = 0; col < cols; ++col) {
                int dist = abs((int)row - (int)q.center_y) + abs((int)col - (int)q.center_x);
                if (dist <= R && ((row + col) % 2 == 0)) { // even sublattice check (|0⟩ neighbors)
                    int i = row * cols + col;
                    pop_sum += std::norm(phi(i * D, 0)); // probability amplitude squared of |0⟩ component
                    count++;
                }
            }
        }

        if (count == 0) return 0.0;
        return pop_sum / count; // average population on even sublattice near logical qubit
    }

private:
    //helper functions

    cx_mat mat_vec_mult_cl10(const sp_cx_mat& H, const cx_mat& phi) {
        return H * phi;
    }

    double inner_product_cl10(const cx_mat& phi1, const cx_mat& phi2) {
        cx_double prod = as_scalar(phi1.t() * phi2);
        return real(prod);
    }

    void compute_nonzero_indices(double J, double ht, int rows, int cols, int D,
                                umat& locations, cx_vec& values, uint& nz) {
        const int N = rows * cols;
        const int NNZ = 6 * N + 2 * rows * cols + 2 * rows * cols; // 5400
        nz = 0;

        if (locations.n_cols < NNZ || values.n_elem < NNZ) {
            cerr << "Error: Arrays too small. Required: " << NNZ
                 << ", Got: " << locations.n_cols << "\n";
            return;
        }

        for (int i = 0; i < N; i++) {
            locations(0, nz) = i * D;     locations(1, nz) = i * D + 1; values(nz) = -ht; nz++;
            locations(0, nz) = i * D + 1; locations(1, nz) = i * D;     values(nz) = -ht; nz++;
        }

        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < cols; col++) {
                int i = row * cols + col;
                int j_right = row * cols + ((col + 1) % cols);
                int j_down = ((row + 1) % rows) * cols + col;

                if (j_right > i || (j_right < i && col == cols - 1)) {
                    locations(0, nz) = i * D;     locations(1, nz) = j_right * D;     values(nz) = -J; nz++;
                    locations(0, nz) = i * D + 1; locations(1, nz) = j_right * D + 1; values(nz) = -J; nz++;
                }

                if (j_down > i || (j_down < i && row == rows - 1)) {
                    locations(0, nz) = i * D;     locations(1, nz) = j_down * D;     values(nz) = -J; nz++;
                    locations(0, nz) = i * D + 1; locations(1, nz) = j_down * D + 1; values(nz) = -J; nz++;
                }
            }
        }
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

    void compute_nonzero_indices_spiral_twist(double J, double ht, int rows, int cols, int D, double omega_ang,
                                              umat& locations, cx_vec& values, uint& nz) {
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
