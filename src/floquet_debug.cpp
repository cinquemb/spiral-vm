#include <armadillo>
#include <complex>
#include <iostream>
#include <fstream>
using namespace arma;
using namespace std;

const int N = 900;
const int D = 2;
const int ROWS = 30;
const int COLS = 30;

/*

Requirements: Armadillo installed (e.g., via libarmadillo-dev on Ubuntu) and linked with LAPACK/BLAS.
Compile: g++ -O2 floquet_2d.cpp -o floquet_2d -larmadillo `pkg-config lapack --libs` `pkg-config blas --libs`
Compile: g++ -O2 floquet_2d_simple.cpp -o floquet_2d_simple -larmadillo `pkg-config lapack --libs` `pkg-config blas --libs`
Run: ./floquet_2d

sudo apt-get install libopenblas-openmp-dev libarmadillo-dev libblas-dev liblapack-dev gfortran

*/

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
    double theta_max = theta_max_base * (1.0 - cos(M_PI * period / 2.0)) / 2.0; // Peaks at 2T, 6T, etc.

    double center_x = COLS / 2.0, center_y = ROWS / 2.0;
    double r_max = sqrt(center_x * center_x + center_y * center_y);
    double norm_factor = sqrt(ROWS * COLS);

    for (int row = 0; row < ROWS; row++) {
        for (int col = 0; col < COLS; col++) {
            int i = row * COLS + col;
            int j_right = row * COLS + ((col + 1) % COLS);
            int j_down = ((row + 1) % ROWS) * COLS + col;

            double r = sqrt(pow(row - center_y, 2) + pow(col - center_x, 2));
            double phi_ang = atan2(row - center_y, col - center_x);
            double theta = (omega_ang == 0) ? 0 : theta_max * (r / r_max) * cos(omega_ang * phi_ang);
            cx_double J_twist;
            if (is_ang == true)
                J_twist = J * cx_double(0, sin(theta));
            else
                J_twist = J;

            cx_double z_i = (phi(i * D, 0) - phi(i * D + 1, 0)) * norm_factor;
            cx_double z_jr = (phi(j_right * D, 0) - phi(j_right * D + 1, 0)) * norm_factor;
            cx_double z_jd = (phi(j_down * D, 0) - phi(j_down * D + 1, 0)) * norm_factor;

            if (is_ang) {
                energy += imag(J_twist * z_i * z_jr); // Spiral: imaginary part
                energy += imag(J_twist * z_i * z_jd);
            } else {
                energy += real(J_twist * z_i * z_jr); // No spiral: real part
                energy += real(J_twist * z_i * z_jd);
            }
        }
    }
    return energy;
}

cx_mat compute_zz_energy_vector(const cx_mat& phi, double J, double omega_ang, double period, bool is_ang = false) {
    const int N = ROWS * COLS;
    cx_mat Hzz_phi = zeros<cx_mat>(N * D, 1);
    double theta_max_base = M_PI / 512;
    double theta_max = theta_max_base * (1.0 - cos(M_PI * period / 2.0)) / 2.0; // Peaks at 2T, 6T, etc.
    double center_x = COLS / 2.0, center_y = ROWS / 2.0;
    double r_max = sqrt(center_x * center_x + center_y * center_y);
    double norm_factor = sqrt(N);

    for (int row = 0; row < ROWS; row++) {
        for (int col = 0; col < COLS; col++) {
            int i = row * COLS + col;
            int j_right = row * COLS + ((col + 1) % COLS);
            int j_down = ((row + 1) % ROWS) * COLS + col;

            double r = sqrt(pow(row - center_y, 2) + pow(col - center_x, 2));
            double phi_ang = atan2(row - center_y, col - center_x);
            double theta = (omega_ang == 0) ? 0 : theta_max * (r / r_max) * cos(omega_ang * phi_ang);
            cx_double J_twist = (is_ang) ? J * cx_double(0, sin(theta)) : J;

            // Fixed: Use raw amplitudes, not norms
            cx_double z_jr = (phi(j_right * D, 0) - phi(j_right * D + 1, 0)) * norm_factor;
            cx_double z_jd = (phi(j_down * D, 0) - phi(j_down * D + 1, 0)) * norm_factor;

            Hzz_phi(i * D, 0) = J_twist * (z_jr + z_jd) * phi(i * D, 0);     // sigma^z_i = +1
            Hzz_phi(i * D + 1, 0) = J_twist * (-z_jr - z_jd) * phi(i * D + 1, 0); // sigma^z_i = -1
        }
    }
    return Hzz_phi;
}

void compute_nonzero_indices_spiral_twist(double J, double ht, int rows, int cols, int D, double omega_ang,
                                          umat& locations, cx_vec& values, uint& nz) {
    const int N = rows * cols;
    const int NNZ = 2 * N; // 2N for sigma^x, 4N for sigma^z sigma^z (4 bonds per site)

    nz = 0;

    // Transverse field terms (-ht sigma^x)
    for (int i = 0; i < N; i++) {
        locations(0, nz) = i * D;     locations(1, nz) = i * D + 1; values(nz) = -ht; nz++;
        locations(0, nz) = i * D + 1; locations(1, nz) = i * D;     values(nz) = -ht; nz++;
    }
}

sp_cx_mat hamiltonian_cl10_90(double J, double ht) {
    int NNZ = 2 * N + 2 * ROWS * COLS + 2 * ROWS * COLS;
    umat locations(2, NNZ);
    cx_vec values(NNZ);
    uint nz;
    compute_nonzero_indices(J, ht, ROWS, COLS, D, locations, values, nz);
    return sp_cx_mat(locations.submat(0, 0, 1, nz - 1), values.subvec(0, nz - 1), N * D, N * D);
}

sp_cx_mat hamiltonian_cl10_90_spiral_twist(double J, double ht, double omega_ang) {
    int N = ROWS * COLS;
    const int NNZ = 2 * N; // 2N for sigma^x, 2N for sigma^z sigma^z
    umat locations(2, NNZ);
    cx_vec values(NNZ);
    uint nz;
    compute_nonzero_indices_spiral_twist(J, ht, ROWS, COLS, D, omega_ang, locations, values, nz);
    return sp_cx_mat(locations.submat(0, 0, 1, nz - 1), values.subvec(0, nz - 1), N * D, N * D);
}

double compute_avg_stabilizer(const cx_mat& phi) {
    double stab_sum = 0.0;
    int count = 0;
    for (int row = 0; row < ROWS; row++) {
        for (int col = 0; col < COLS; col++) {
            int i = row * COLS + col;
            int j_right = row * COLS + ((col + 1) % COLS);
            int j_down = ((row + 1) % ROWS) * COLS + col;
            double zi = abs(phi(i*D,0) - phi(i*D+1,0));
            double zj_right = abs(phi(j_right*D,0) - phi(j_right*D+1,0));
            double zj_down = abs(phi(j_down*D,0) - phi(j_down*D+1,0));
            stab_sum += zi * zj_right + zi * zj_down;
            count += 2;
        }
    }
    return stab_sum / count; // ~1 for perfect order
}

void floquet_cl10_90(double J, double h0, double h1, double omega, double T, int steps, int N_max, double sx_gain, bool is_ang) {
    double dt = T / steps;

    string initial_state = "neel";
    cx_mat phi = zeros<cx_mat>(N * D, 1);
    mt19937 rng(777);
    uniform_real_distribution<double> dist(0.0, 1.0);
    for (int row = 0; row < ROWS; row++) {
        for (int col = 0; col < COLS; col++) {
            int i = row * COLS + col;
            if (initial_state == "neel") {
                if ((row + col) % 2 == 0) {
                    phi(i * D, 0) = 1.0; // |0>
                } else {
                    phi(i * D + 1, 0) = 1.0; // |1>
                }
            } else if (initial_state == "polarized") {
                phi(i * D, 0) = 1.0; // |↑⟩
            } else if (initial_state == "disordered") {
                if (dist(rng) < 0.5) {
                    phi(i * D, 0) = 1.0; // |↑⟩
                } else {
                    phi(i * D + 1, 0) = 1.0; // |↓⟩
                }
            }
        }
    }


    double norm0 = sqrt(inner_product_cl10(phi, phi)); // ~46.1111
    phi /= norm0;
    cx_mat phi_in = phi;
    double h1_base = h1;

    //h1 energy gain
    //double sx_gain = 1900;//1300.0; max limit

    double h1_limit = 10000.00;

    //h1 feedback gain
    double h1_f_gain  = 0.0;
    double k_bA = 0.0; //3 is too much, 0.4 is too little; 0.8 - 1.6 is ok;
    
    //angualr freq
    double k_oA = 1.0;//angualr amp
    double alpha_ang = 1.0;    //rotatinal amp

    double omega_ang_base = 0.0;
    
    double P1 = T;
    double P2 = T;
    double beta_ang = 0.0;
    double k_angular = 0.0;
    double P = T;
    double omega_omega_angA = omega;
    double omega_omega_angB = 2*omega_omega_angA;

    if (is_ang == false) {
        alpha_ang = 0;
        k_oA = 0;
        omega_omega_angA = 0;
        omega_omega_angB = 0;
    }

    const int W = 32;
    vector<double> fidelity_window(W, 0.0);
    double delta_F_target = 1.0;
    vec fidelities(N_max + 1, fill::zeros);
    fidelities(0) = 1.0;


    //spiral
    double omega_ang_init = omega_ang_base;
    sp_cx_mat H_init = hamiltonian_cl10_90_spiral_twist(J, 0, 0);
    cx_mat Hphi_init = mat_vec_mult_cl10(H_init, phi);
    double energy_init = real(inner_product_cl10(phi, Hphi_init)) + compute_zz_energy(phi, J, 0, 0); // Full energy
    double zz_energy_init = compute_zz_energy(phi, J, 0, 0);
    /*no spiral energy
    sp_cx_mat H_init = hamiltonian_cl10_90_spiral_twist(J, 0, 2.0); // J-term only
    cx_mat Hphi_init = mat_vec_mult_cl10(H_init, phi);
    double energy_init = real(inner_product_cl10(phi, Hphi_init));

    */

    cout << "Initial norm = " << norm0 << "\n";
    cout << "phi[0] initial = " << phi(0, 0) << "\n";
    cout << "Initial energy = " << energy_init << "\n";


    std::ostringstream os;
    //os << "dtc_floquet_with" << ((is_ang == false) ? "_no_spiral_" : "_spiral_") << (int)sx_gain << ".txt";
    os << "dtc_floquet_with_" << initial_state << "_state_" << ((is_ang == false) ? "_no_spiral_" : "_spiral_") << (int)sx_gain << "_omega_" << omega << "_omegaang_" << omega_ang_base << ".txt";
    std::string fname = os.str();

    double sx_energy_init = 0.0;
    for (int i = 0; i < N; i++) {
        sx_energy_init -= h1 * 2.0 * real(phi(i * D, 0) * conj(phi(i * D + 1, 0)));
    }


    ofstream fout(fname);
    fout << "Period,Fidelity,Stabilizer,Energy,sx_energy,zz_energy,Delta_F,ht_eff_end,sx_avg\n";
    fout << "0,1.0," << compute_avg_stabilizer(phi) << "," << energy_init << ","<< sx_energy_init<<","<<zz_energy_init <<"," << h1 << ",0,0\n";
    
    double delta_F = 0.0; // Initial fidelity
    double sx_prev = 0.0;
    for (int n = 0; n < N_max; n++) {
        double t_n = n * T; // Base time for this period

        cx_mat phi_new = phi;
        for (int k = 0; k < steps; k++) {
            double t = t_n + k * dt;
            double angular_freq_quasi = alpha_ang * (sin(omega_omega_angA * M_PI * t / P1) + sin(omega_omega_angB * M_PI * t / P2));
            double angular_freq_feedback = beta_ang * (delta_F_target - delta_F);
            double omega_ang_mod = omega_ang_base + k_oA * sin(omega_omega_angA * M_PI * t / P1);
            double omega_ang = omega_ang_mod + angular_freq_quasi + k_angular * angular_freq_feedback;
            double h1_feedback = h1_f_gain * (k_bA * (delta_F_target - delta_F));
            double ht_eff = h0 + h1 * cos(omega * (t / 2) + M_PI / 4.0) + h1_feedback;

            if (ht_eff / J > h1_limit) ht_eff = h1_limit * J;
            else if (ht_eff / J < -h1_limit) ht_eff = -h1_limit * J;

            sp_cx_mat H_sx = hamiltonian_cl10_90_spiral_twist(J, ht_eff, omega_ang);
            if (is_ang == false) {
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

        double angular_freq_quasi_end = alpha_ang * (sin(omega_omega_angA * M_PI * (n+1) * T / P1) + sin(omega_omega_angB * M_PI * (n+1) * T / P2));
        double angular_freq_feedback_end = beta_ang * (delta_F_target - delta_F);
        double omega_ang_mod_end = omega_ang_base + k_oA * sin(omega_omega_angA * M_PI * (n+1) + T / P);
        double omega_ang_end = omega_ang_mod_end + angular_freq_quasi_end + k_angular * angular_freq_feedback_end;
        double h1_feedback_end = h1_f_gain * (k_bA * (delta_F_target - delta_F));
        double ht_eff_end = h0 + h1 * cos(omega * ((n+1) * T / 2) + M_PI / 4.0) + h1_feedback_end;
        if (ht_eff_end/J > h1_limit) ht_eff_end = h1_limit * J;
        else if (ht_eff_end/J < -h1_limit) ht_eff_end = -h1_limit * J;


        sp_cx_mat H_energy;
        if (is_ang == false) {
            H_energy = hamiltonian_cl10_90_spiral_twist(J, ht_eff_end, 0); // Energy at h1
        } else {
            H_energy = hamiltonian_cl10_90_spiral_twist(J, ht_eff_end, omega_ang_end); // Energy at h1
        } 


        fidelities(n + 1) = abs(inner_product_cl10(phi_in, phi));
        fidelity_window[n % W] = fidelities(n + 1);


        double stabilizer = compute_avg_stabilizer(phi);
        cx_mat Hphi = mat_vec_mult_cl10(H_energy, phi);

        double sx_energy = 0.0;
        for (int i = 0; i < N; i++) {
            sx_energy -= ht_eff_end * 2.0 * real(phi(i * D, 0) * conj(phi(i * D + 1, 0)));
        }
        
        double zz_energy;
        if (is_ang == false) {
            zz_energy = compute_zz_energy(phi, J, 0, 0);
        } else {
            zz_energy = compute_zz_energy(phi, J, omega_ang_end, (n+1), true);
        }
        double energy = sx_energy + zz_energy;

        double sx_sum = 0.0;
        for (int i = 0; i < N; i++) {
            sx_sum += 2.0 * real(phi(i * D, 0) * conj(phi(i * D + 1, 0)));
        }
        double sx_avg = sx_sum / N;


        fout << n + 1 << "," << fidelities(n + 1) << "," << stabilizer << "," << energy << "," << sx_energy << "," << zz_energy << "," << delta_F << "," << ht_eff_end << "," << sx_avg << "\n";
        cout << "Period " << n + 1 << ": Fidelity = " << fidelities(n + 1) << ", Stabilizer = " << stabilizer 
             << ", Energy = " << energy << ", Delta_F = " << delta_F << "\n";
        cout << "sx_energy = " << sx_energy << ", zz_energy = " << zz_energy << ", sx_avg = " << sx_avg << ", ht_eff_end = " << ht_eff_end  << "\n";
        if (n % 2 == 1 && fidelities(n + 1) > 0.9) {
            cout << "Period-doubling at " << n + 1 << "T, Fidelity = " << fidelities(n + 1) << "\n";
        }

        /*
        if (fidelities(n + 1) < 0.9) {
            break;//break early if less than threshold
        }*/
    }
    fout.close();
}

int main() {
    //double J = 0.3, h0 = 0, h1 =2.5/4.4, omega = 20.0, T =  datum::pi;//optimized
    double J = 0.3, h0 = 0, h1 =2.5/4.4, omega = 20 * 2 * datum::pi; double T =  2 * datum::pi / omega;
    //double J = 4.4, h0 = 0, h1 =2.5, omega = 2.0, T =  datum::pi;//unoptimized

    vector<double> sx_gains = {0, 500, 1000, 1500, 1900, 2500, 3000};
    //vector<double> omegas = {20, 50, 100, 150, 200, 250, 300, 350, 400};
    vector<double> omegas = {20, 50, 100, 150, 200, 250, 300, 350, 400};
    int steps = 200, N_max = 5000;   double sx_gain = 1900;
    bool is_ang = true;
    /*
    for(int i =0; i<sx_gains.size();i++){
        floquet_cl10_90(J, h0, h1, omega, T, steps, N_max, sx_gains[i]);
    }*/
    floquet_cl10_90(J, h0, h1, omega, T, steps, N_max, sx_gain, is_ang);

    /*
    for(int i =17000; i<=30000;i+=50){
        floquet_cl10_90(J, h0, h1, i, T, steps, N_max, sx_gain, is_ang);
    }*/

    return 0;
}