// spiral_compiler.hpp — the first header of quantum Linux
struct LogicalQubit {
    uint32_t center_x, center_y;  // spiral center
    double   base_phase;
};

struct Gate {
    enum Type { X, Z, CZ, T, PHASE, MEASURE };
    Type     type;
    uint32_t target;
    uint32_t control;  // for CZ
    double   angle;
};

class SpiralVM {
    std::vector<LogicalQubit> qubits;
    double J, h0, h1, omega;
    double omega_ang_base = 126.0;   // your magic number
    double sx_gain = 1900.0;

public:
    void add_qubit(int x, int y);
    void apply_gate(const Gate& g, double period_time);
    void compile_and_run(const std::vector<Gate>& program);
};