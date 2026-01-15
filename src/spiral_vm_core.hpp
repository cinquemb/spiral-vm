#ifndef SPIRAL_COMPILER_HPP
#define SPIRAL_COMPILER_HPP

#include <vector>
#include <cstdint>
#include <armadillo>
#include <complex>
#include <random>
#include <functional>
#include <string>
#include <limits>



// Minimal footprint Waveform/Tone types (kept public for test scripting)
struct Tone {
    double amp;      
    double freq;     
    double phase;    
    double I_component, Q_component;  // Vector components
    double envelope_start; 
    double envelope_end;   
    int logical_id;   

    // Updated Constructor
    Tone(double a=0, double f=0, double p=0, double s=0, double e=1.0, int q=-1)
        : amp(a), 
          freq(f), 
          phase(p),
          I_component(a), // Default: I is the amplitude, Q is 0 (Pure X-axis)
          Q_component(0.0),
          envelope_start(s), 
          envelope_end(e),
          logical_id(q) {}
          
    // Optional helper to set IQ directly
    void set_iq(double i, double q) {
        I_component = i;
        Q_component = q;
        amp = std::sqrt(i*i + q*q); // Keep amp synced for legacy code
    }
};


struct Waveform {
    std::vector<Tone> tones;
    // Evaluate waveform at time t (seconds)
    // Optional lowpass anti-aliasing applied externally by engine
    double eval(double t) const {
        double s = 0.0;
        for (const auto &tn : tones) {
            // basic envelope: linear ramp between start and end (wrap around handled by caller)
            // the envelope here is naive and expects caller to choose appropriate durations
            s += tn.amp * std::cos(tn.freq * t + tn.phase);
        }
        return s;
    }
};

struct LogicalQubit {
    uint32_t center_x, center_y;  // spiral center
    double base_phase;
    // allocator info
    int waveform_id = -1;
};

struct Gate {
    enum Type { X, Z, CZ, H, T, RZ, RX, RY, CNOT, PHASE, MEASURE };
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
    double LOGICAL_X_AMPLITUDE = 125;///1396999.7245;
    double PHASE_RAMP_MAGNITUDE = 0.5;


    static constexpr int D = 2;
    const int R; // Physical neighborhood radius around each logical qubit's center
    const int rows, cols;         // Lattice dimensions
    const int N;                  // Number of sites
    double J, h0, h1, omega, T;  // Hamiltonian / Floquet parameters
    bool is_ang;                 // Spiral angle flag
    bool overlap_enabled = false; // Overlap mode toggle
    arma::cx_mat state;   // ALWAYS up-to-date state (used by all measurements)

    // Waveform engine tunables (public for easy experimentation)
    double freq_base = 2.0;       // base addressable angular frequency (rad/s)
    double freq_spacing = 0.002;  // spacing between logical qubit carriers (rad/s)
    double lowpass_cutoff = 0.5;  // low-pass factor (0..1) relative to Nyquist (coarse)
    double max_tone_amp = 100000000000000000.0;    // safety clamp on tone amplitudes

    // Logical qubit management
    uint32_t add_qubit(uint32_t x, uint32_t y);

    // Gate scheduling and application
    void apply_gate(const Gate& g, double period_time = 0.0);
    void compile_and_run(const std::vector<Gate>& program);

    // Functional control
    void run_floquet(int N_max, const std::string& initial_state);
    void run_periods(uint32_t n);
    void apply_global_pi_pulse_on_even_cycles();
    void global_phase_ramp(double slope, int steps);
    void logical_phase_ramp(int target_qid, double slope, int steps);
    double measure_even_population(uint32_t qid);
    void apply_phase_shift(double angle);
    void apply_phase_kick_between(uint32_t qid1, uint32_t qid2, double strength, double duration_fraction);
    void apply_phase_kick_between_full(uint32_t qid1, uint32_t qid2, double strength, double duration_fraction);
    double logical_zz_correlation(uint32_t qid1, uint32_t qid2) const;
    double get_logical_phase(uint32_t qid);
    double measure_logical_Z(uint32_t qid) const;
    double measure_logical_global_Z(uint32_t qid) const;
    double measure_logical_X(uint32_t qid) const;
    double measure_logical_Y(uint32_t qid) const;
    void ramp_omega_ang(double start, double end, double duration_seconds);
    void global_pi_pulse();
    void logical_hadamard(uint32_t qid);
    void apply_T_gate(uint32_t qid, int steps);
    double get_logical_phase_frame_corrected(uint32_t qid);
    int find_waveform_index_for_qubit(uint32_t qid);
    void logical_x_pulse(uint32_t qid, double duration_periods);
    void logical_cz(uint32_t control, uint32_t target);
    void logical_z_rotation(uint32_t qid, double angle);
    void logical_controlled_phase(uint32_t control, uint32_t target,
                                        double max_angle, double duration_periods);

    // Printing Wavforms
    void dump_waveforms(const std::string& format = "csv",  // "csv" or "json"
        const std::string& prefix = "waveform_",
                    int period = -1) const;               // -1 = current/latest
    void dump_h_eff(const std::string& fname_base, int period = -1) const;
    size_t find_carrier_tone(int wid, uint32_t qid) const;

    // Helper to sample a waveform over one period
    void sample_waveform(const Waveform& w, double t_start, double dt, arma::vec& times, arma::cx_vec& iq, arma::vec& amps, arma::vec& phases) const;
    void compile_to_physical_waveform();  // Compiles all logical waveforms into a single global physical waveform (multi-tone broadcast)
    void dump_frequency_mapping(const std::string& fname = "frequency_to_logical.json") const;

    std::pair<double, double> reconstruct_logical_amp_phase_from_csv(const std::string& fname, uint32_t qid, bool from_file);

    // Initialization and simulation control
    void initialize_state(const std::string& initial_state = "neel");
    double omega_ang_end(int n) const;
    double h_effective_end(int n);
    double sx_avg(int n);
    int get_period();
    void print_overlap_stats();

    // Add these 2 methods
    void virtual_phase_gate(uint32_t qid, double angle);  // Virtual T/Z: NO decoherence
    double measure_logical_Z_frame_corrected(uint32_t qid) const;  // Applies virtual compensation
    int get_total_logical_qubits();
    size_t find_carrier_tone(int wid, uint32_t qid);

private:
    // internal state
    double omega_ang_base;
    double drive_phase = 0.0;

    bool auto_compile_enabled = true;  // member variable


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

    // Waveform engine (bank + mapping)
    std::vector<Waveform> waveforms;       // global waveform bank
    Waveform physical_waveform;               // index 0 only: merged physical;
    std::vector<int> drive_index;          // size N, drive_index[i] = waveform ID


    std::vector<double> virtual_frame_phase;  // VIRTUAL Z: infinite T2 bookkeeping


    // Frequency allocation bookkeeping
    std::vector<double> allocated_carriers;

    // Private methods: physics calculations, Floquet step
    void step_period(int n, double& delta_F);

    // Waveform helpers
    int allocate_waveform_for_qubit(uint32_t qid);
    Waveform make_default_logical_waveform(uint32_t qid);
    Waveform make_cz_waveform(uint32_t qid1, uint32_t qid2, double strength, double duration_fraction);
    void lowpass_filter_waveform(Waveform &w, double cutoff_factor);
    double eval_waveform_with_envelope(const Waveform &w, double t, double local_period_fraction=0.0) const;

    // Hamiltonian helpers
    arma::cx_mat mat_vec_mult_cl10(const arma::sp_cx_mat& H, const arma::cx_mat& phi);
    double inner_product_cl10(const arma::cx_mat& phi1, const arma::cx_mat& phi2);
    double compute_zz_energy(const arma::cx_mat& phi, double J, double omega_ang, double period, bool is_ang = false);
    double compute_zz_energy_edgeaware(const arma::cx_mat& phi_inp, double J, double omega_ang, double period, bool is_ang = false);
    arma::cx_mat compute_zz_energy_vector(const arma::cx_mat& phi, double J, double omega_ang, double period, bool is_ang = false);
    arma::cx_mat compute_zz_energy_vector_edgeaware(const arma::cx_mat& phi_inp, double J, double omega_ang, double period, bool is_ang = false);
    void compute_nonzero_indices_spiral_twist(double J, double ht, int rows, int cols, int D, double omega_ang, arma::umat& locations, arma::cx_vec& values, uint& nz);
    arma::sp_cx_mat hamiltonian_cl10_90_spiral_twist(double J, double ht, double omega_ang);
    arma::sp_cx_mat hamiltonian_cl10_90_spiral_twist_inhomogeneous(double J, const std::vector<double> &local_hx, double omega_ang);

    int get_right_neighbor(int row, int col) const;
    int get_down_neighbor(int row, int col) const;

    double compute_avg_stabilizer(const arma::cx_mat& phi);
    double current_orbit_phase(uint32_t qid) const;

    // Safety helpers
    double clamp_tone_amp(double a) const;
};

#endif // SPIRAL_COMPILER_HPP
