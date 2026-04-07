#include "q2ns-spiral.hpp"
#include <sstream>
#include <iomanip>

namespace q2ns {

QStateSpiral::QStateSpiral(std::size_t grid_size, bool use_phi_direct)
    : vm_(static_cast<int>(grid_size), static_cast<int>(grid_size)),
      use_phi_direct_(use_phi_direct)
{
    vm_.is_ang = true;
    vm_.use_phi_direct = use_phi_direct;
    vm_.initialize_state("neel");
    vm_.auto_compile_enabled = false;
}

std::size_t QStateSpiral::NumQubits() const {
    return static_cast<std::size_t>(vm_.get_total_logical_qubits());
}

uint32_t QStateSpiral::AddQubit(uint32_t x, uint32_t y) {
    return vm_.add_qubit(x, y);
}

int64_t QStateSpiral::AssignStreams(int64_t stream) {
    constexpr uint64_t SPIRAL_SALT = 0x53504952414C564DULL;
    return AssignStreamsGlobal<SPIRAL_SALT>(stream, [this](uint64_t seed) {
        vm_.rng.seed(seed);
    });
}

void QStateSpiral::Apply(const QGate& g, const std::vector<Index>& targets)
{
    if (targets.empty()) return;

    // Q2NS uses QGateKind, while SpiralVM uses its own Gate struct.
    // We map Q2NS gates to SpiralVM's logical operations where possible.

    if (g.Kind() == QGateKind::Custom) {
        std::cerr << "[SpiralVM] Custom unitary gates not supported (waveform backend)\n";
        return;
    }

    // Most Q2NS gates are single-qubit or two-qubit with fixed targets
    if (targets.size() == 1) {
        uint32_t q = static_cast<uint32_t>(targets[0]);

        switch (g.Kind()) {
            case QGateKind::I:
                break;                                      // no-op

            case QGateKind::H:
                vm_.logical_hadamard(q);                    // uses your tuned Hadamard
                break;

            case QGateKind::X:
                vm_.logical_x_pulse(q, 1.0);
                break;

            case QGateKind::Y:
                vm_.logical_y_pulse(q, 1.0);
                break;

            case QGateKind::Z:
                vm_.logical_z_rotation(q, M_PI);            // full Z
                break;

            case QGateKind::S:
                vm_.logical_phase_ramp(q, M_PI/2.0, 2);
                break;

            case QGateKind::SDG:
                vm_.logical_phase_ramp(q, -M_PI/2.0, 2);
                break;

            case QGateKind::T:   // T = Z(π/4)
                vm_.logical_phase_ramp(q, M_PI/4.0, 2);
                break;

            default:
                std::cerr << "[SpiralVM] Unsupported single-qubit gate: " 
                          << static_cast<int>(g.Kind()) << "\n";
        }
    }
    else if (targets.size() >= 2) {
        // Two-qubit gates – assume first is control, second is target (common convention)
        uint32_t ctrl = static_cast<uint32_t>(targets[0]);
        uint32_t tgt  = static_cast<uint32_t>(targets[1]);

        switch (g.Kind()) {
            case QGateKind::CZ:
                vm_.logical_cz(ctrl, tgt);
                break;

            case QGateKind::CNOT:
                // Standard decomposition: H(tgt) - CZ - H(tgt)
                vm_.logical_hadamard(tgt);
                vm_.logical_cz(ctrl, tgt);
                vm_.logical_hadamard(tgt);
                break;

            case QGateKind::SWAP:
                // Rough 3-CNOT decomposition
                vm_.logical_hadamard(tgt); vm_.logical_cz(ctrl, tgt); vm_.logical_hadamard(tgt);
                vm_.logical_hadamard(ctrl); vm_.logical_cz(tgt, ctrl); vm_.logical_hadamard(ctrl);
                vm_.logical_hadamard(tgt); vm_.logical_cz(ctrl, tgt); vm_.logical_hadamard(tgt);
                break;

            default:
                std::cerr << "[SpiralVM] Unsupported two-qubit gate: " 
                          << static_cast<int>(g.Kind()) << "\n";
        }
    }

    // Always run a stabilization period after gate application (as in your apply_gate)
    vm_.run_periods(1);
}

std::shared_ptr<QState> QStateSpiral::MergeDisjoint(const QState& other) const {
    auto copy = std::make_shared<QStateSpiral>(*this);

    const auto* other_spiral = dynamic_cast<const QStateSpiral*>(&other);
    if (!other_spiral) {
        std::cerr << "[SpiralVM] MergeDisjoint: other is not QStateSpiral\n";
        return copy;
    }

    uint32_t current_max_x = 0;
    for (const auto& q : vm_.logical_qubits) {
        if (q.center_x > current_max_x) current_max_x = q.center_x;
    }

    uint32_t offset_x = current_max_x + 6;

    for (uint32_t i = 0; i < other_spiral->NumQubits(); ++i) {
        copy->AddQubit(offset_x + (i % 10)*4, 12 + (i / 10)*4);
    }

    return copy;
}

// Improved collapse: rotate toward the measured basis before projecting
void QStateSpiral::ApproximateCollapse(uint32_t qid, Basis basis, int outcome)
{
    double angle = 0.0;
    double axis  = 0.0;

    if (basis == Basis::X)      { angle = M_PI; axis = 0.0; }      // X rotation
    else if (basis == Basis::Y) { angle = M_PI; axis = M_PI/2.0; } // Y rotation
    else /* Z */                { angle = (outcome == 1 ? M_PI : 0.0); axis = 0.0; }

    vm_.apply_local_rotation(qid, angle, axis);
    vm_.run_periods(2);
}

MeasureResult QStateSpiral::Measure(Index target, Basis basis)
{
    uint32_t q = static_cast<uint32_t>(target);
    double expectation = 0.0;

    switch (basis) {
        case Basis::Z:
            expectation = vm_.measure_logical_Z(q);
            break;
        case Basis::X:
            expectation = vm_.measure_logical_X(q);
            break;
        case Basis::Y:
            expectation = vm_.measure_logical_Y(q);
            break;
        default:
            std::cerr << "[SpiralVM] Unsupported measurement basis\n";
            expectation = vm_.measure_logical_Z(q);
    }

    double p0 = (1.0 + expectation) / 2.0;
    int outcome = (vm_.dist(vm_.rng) < p0) ? 0 : 1;

    ApproximateCollapse(q, basis, outcome);

    MeasureResult res;
    res.outcome = outcome;
    res.measured = std::make_shared<QStateSpiral>(*this);
    res.survivor = std::make_shared<QStateSpiral>(*this);   // TODO: proper removal

    return res;
}

void QStateSpiral::Print(std::ostream& os) const {
    PrintHeader(os, "SpiralVM (Waveform / Phase-Kick)");

    os << "  Logical qubits : " << NumQubits() << "\n";
    os << "  Use phi direct : " << (use_phi_direct_ ? "yes" : "no") << "\n";
    os << "  Current period : " << vm_.get_period() << "\n";
    os << "  Underlying state: full lattice phi (" << vm_.N << " sites)\n";
    os << "  Crosstalk regime: negligible below 121 logical : 1 physical\n";
    os << "                    spectral crowding starts at 122:1\n\n";

    os << "  Qubit  |   <Z>     |   <X>     |   <Y>     |   Phase (×2π)\n";
    os << "  -------+-----------+-----------+-----------+---------------\n";

    for (std::size_t q = 0; q < NumQubits(); ++q) {
        uint32_t qid = static_cast<uint32_t>(q);
        double z = vm_.measure_logical_Z(qid);
        double x = vm_.measure_logical_X(qid);
        double y = vm_.measure_logical_Y(qid);
        double phi_val = vm_.get_logical_phase(qid);

        os << "  " << std::setw(6) << q
           << " | " << std::fixed << std::setprecision(4) << std::setw(9) << z
           << " | " << std::setw(9) << x
           << " | " << std::setw(9) << y
           << " | " << std::setw(12) << (phi_val / (2.0 * M_PI)) << "\n";
    }
}

std::shared_ptr<QState> CreateQStateSpiral(std::size_t grid_size, bool use_phi_direct)
{
    return std::make_shared<QStateSpiral>(grid_size, use_phi_direct);
}

} // namespace q2ns