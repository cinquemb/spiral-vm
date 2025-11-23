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

    // Logical qubit management
    uint32_t add_qubit(int x, int y);

    // Gate scheduling and application
    void apply_gate(const Gate& g, double period_time);
    void compile_and_run(const std::vector<Gate>& program);

    // Functional control
    void run_periods(uint32_t n);
    void apply_global_pi_pulse_on_even_cycles();
    double measure_even_population(uint32_t qid);
    void apply_phase_shift(double angle);
    void apply_phase_kick_between(uint32_t qid1, uint32_t qid2, double strength, double duration_fraction);
    void apply_phase_kick_between_full(uint32_t qid1, uint32_t qid2, double strength, double duration_fraction);
    double logical_zz_correlation(uint32_t qid1, uint32_t qid2);
    double get_logical_phase(uint32_t qid);
    void ramp_omega_ang(double start, double end, double duration_seconds);

    // Initialization and simulation control
    void initialize_state(const std::string& initial_state = "neel");

private:
    int rows_, cols_;
    int N_;
    double J_, h0_, h1_, omega_, T_;
    double omega_ang_base_;
    double sx_gain_;

    arma::cx_mat phi_;
    arma::cx_mat phi_in_;
    int steps_;
    int current_period_;

    std::vector<LogicalQubit> logical_qubits_;

    // Random engine for state initialization
    std::mt19937 rng_;
    std::uniform_real_distribution<double> dist_;

    // Fidelity tracking
    std::vector<double> fidelities_;
    std::vector<double> fidelity_window_;

    // Private methods: physics calculations, Floquet step
    void step_period(int n, double& delta_F);

    void global_pi_pulse();

    // Helpers for energy, Hamiltonian, etc.

    // e.g.:
    arma::cx_mat mat_vec_mult_cl10(const arma::sp_cx_mat& H, const arma::cx_mat& phi);
    double inner_product_cl10(const arma::cx_mat& phi1, const arma::cx_mat& phi2);
    double compute_zz_energy(const cx_mat& phi, double J, double omega_ang, double period, bool is_ang = false);
    cx_mat compute_zz_energy_vector(const cx_mat& phi, double J, double omega_ang, double period, bool is_ang = false);
    void compute_nonzero_indices_spiral_twist(double J, double ht, int rows, int cols, int D, double omega_ang, umat& locations, cx_vec& values, uint& nz);
    arma::sp_cx_mat hamiltonian_cl10_90_spiral_twist(double J, double ht, double omega_ang);
    double compute_avg_stabilizer(const cx_mat& phi);
};

#endif // SPIRAL_COMPILER_HPP
