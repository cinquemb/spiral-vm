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
    // Copy your helper functions here or keep in separate utility files

    // e.g.:
    arma::cx_mat mat_vec_mult_cl10(const arma::sp_cx_mat& H, const arma::cx_mat& phi);
    double inner_product_cl10(const arma::cx_mat& phi1, const arma::cx_mat& phi2);
    arma::sp_cx_mat hamiltonian_cl10_90_spiral_twist(double J, double ht, double omega_ang);

    // etc...
};

#endif // SPIRAL_COMPILER_HPP
