// spiral_compiler.hpp

#ifndef SPIRAL_COMPILER_HPP
#define SPIRAL_COMPILER_HPP

#include <vector>
#include <cstdint>
#include <armadillo>

struct LogicalQubit {
    uint32_t center_x, center_y;  // spiral center
    double base_phase;
};

struct Gate {
    enum Type { X, Z, CZ, T, PHASE, MEASURE };
    Type type;
    uint32_t target;
    uint32_t control;  // only for CZ
    double angle;

    Gate(Type t, uint32_t tgt, double ang = 0.0, uint32_t ctrl = UINT32_MAX)
        : type(t), target(tgt), control(ctrl), angle(ang) {}
};

class SpiralVM {
public:
    SpiralVM(int rows, int cols);  // constructor

    static constexpr int D = 2;
    const int R = 1; // Physical neighborhood radius around each logical qubit's center
    const int rows, cols;         // Lattice dimensions
    const int N;                  // Number of sites
    double J, h0, h1, omega, T;  // Hamiltonian / Floquet parameters
    bool is_ang;                 // Spiral angle flag
    bool overlap_enabled = false; // Overlap mode toggle
    arma::cx_mat state;   // ALWAYS up-to-date state (used by all measurements)

    // Logical qubit management
    uint32_t add_qubit(uint32_t x, uint32_t y);

    // Gate scheduling and application
    void apply_gate(const Gate& g, double period_time);
    void compile_and_run(const std::vector<Gate>& program);

    // Functional control
    void run_floquet(int N_max, const std::string& initial_state);
    void run_periods(uint32_t n);
    void apply_global_pi_pulse_on_even_cycles();
    double measure_even_population(uint32_t qid);
    void apply_phase_shift(double angle);
    void apply_phase_kick_between(uint32_t qid1, uint32_t qid2, double strength, double duration_fraction);
    void apply_phase_kick_between_full(uint32_t qid1, uint32_t qid2, double strength, double duration_fraction);
    double logical_zz_correlation(uint32_t qid1, uint32_t qid2);
    double get_logical_phase(uint32_t qid);
    double measure_logical_Z(uint32_t qid) const;
    void ramp_omega_ang(double start, double end, double duration_seconds);
    void global_pi_pulse();

    // Initialization and simulation control
    void initialize_state(const std::string& initial_state = "neel");
    double omega_ang_end(int n);
    double h_effective_end(int n);
    double sx_avg(int n);
    void print_overlap_stats();

private:
    double omega_ang_base;

    arma::cx_mat phi;                  // Quantum state vector (2*N x 1)
    arma::cx_mat phi_in;               // Initial state for fidelity measurement
    int steps;                   // RK4 steps per period
    int current_period;          // Tracks Floquet periods elapsed
    double sx_gain;              // h1 gain parameter

    std::vector<LogicalQubit> logical_qubits;  // Logical qubit list

    // Mapping physical qubits to logical qubits (allows many logicals per physical)
    std::vector<std::vector<uint32_t>> phys_to_logicals;

    // Random engine for state initialization
    std::mt19937 rng;
    std::uniform_real_distribution<double> dist;

    // Fidelity tracking
    std::vector<double> fidelities;
    std::vector<double> fidelity_window;

    // Private methods: physics calculations, Floquet step
    void step_period(int n, double& delta_F);


    // Helpers for energy, Hamiltonian, etc.

    // e.g.:
    arma::cx_mat mat_vec_mult_cl10(const arma::sp_cx_mat& H, const arma::cx_mat& phi);
    double inner_product_cl10(const arma::cx_mat& phi1, const arma::cx_mat& phi2);
    double compute_zz_energy(const arma::cx_mat& phi, double J, double omega_ang, double period, bool is_ang = false);
    arma::cx_mat compute_zz_energy_vector(const arma::cx_mat& phi, double J, double omega_ang, double period, bool is_ang = false);
    void compute_nonzero_indices_spiral_twist(double J, double ht, int rows, int cols, int D, double omega_ang, arma::umat& locations, arma::cx_vec& values, uint& nz);
    arma::sp_cx_mat hamiltonian_cl10_90_spiral_twist(double J, double ht, double omega_ang);
    double compute_avg_stabilizer(const arma::cx_mat& phi);
};

#endif // SPIRAL_COMPILER_HPP
