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
    vm_.initialize_state("neel");   // or "neel" — change as needed
    vm_.auto_compile_enabled = false;    // usually better to control manually
}

std::size_t QStateSpiral::NumQubits() const {
    return static_cast<std::size_t>(vm_.get_total_logical_qubits());
}

uint32_t QStateSpiral::AddQubit(uint32_t x, uint32_t y) {
    return vm_.add_qubit(x, y);
}

int64_t QStateSpiral::AssignStreams(int64_t stream) {
    // Unique salt so SpiralVM does not collide with other backends
    constexpr uint64_t SPIRAL_SALT = 0x53504952414C564DULL; // "SPIRALVM"

    return AssignStreamsGlobal<SPIRAL_SALT>(stream, [this](uint64_t seed) {
        vm_.rng.seed(seed);
        // You can reseed other internal generators here if needed
    });
}

void QStateSpiral::Apply(const QGate& g, const std::vector<Index>& targets)
{
    if (targets.empty()) return;

    // Helper lambda to apply a built-in gate
    auto apply_builtin = [&](QGateKind kind) {
        switch (kind) {
            case QGateKind::I:  /* no-op */ break;
            case QGateKind::H:
                for (Index q : targets) vm_.logical_hadamard(static_cast<uint32_t>(q), 200, 1.0);
                break;
            case QGateKind::X:
                for (Index q : targets) vm_.logical_x_pulse(static_cast<uint32_t>(q), 1.0);
                break;
            case QGateKind::Y:
                for (Index q : targets) vm_.logical_y_pulse(static_cast<uint32_t>(q), 1.0);
                break;
            case QGateKind::Z:
            case QGateKind::S:
                for (Index q : targets) vm_.logical_phase_ramp(static_cast<uint32_t>(q), M_PI/2.0, 2);
                break;
            case QGateKind::SDG:
                for (Index q : targets) vm_.logical_phase_ramp(static_cast<uint32_t>(q), -M_PI/2.0, 2);
                break;
            case QGateKind::CNOT:
            case QGateKind::CZ:
            case QGateKind::SWAP:
                if (targets.size() >= 2) {
                    uint32_t ctrl = static_cast<uint32_t>(targets[0]);
                    uint32_t tgt  = static_cast<uint32_t>(targets[1]);
                    if (kind == QGateKind::CZ || kind == QGateKind::CNOT) {
                        vm_.logical_cz(ctrl, tgt);           // or logical_controlled_phase
                        if (kind == QGateKind::CNOT) {
                            vm_.logical_hadamard(tgt);
                            vm_.logical_cz(ctrl, tgt);
                            vm_.logical_hadamard(tgt);
                        }
                    } else if (kind == QGateKind::SWAP) {
                        // Rough SWAP via 3 CNOTs
                        vm_.logical_hadamard(tgt); vm_.logical_cz(ctrl, tgt); vm_.logical_hadamard(tgt);
                        vm_.logical_hadamard(ctrl); vm_.logical_cz(tgt, ctrl); vm_.logical_hadamard(ctrl);
                        vm_.logical_hadamard(tgt); vm_.logical_cz(ctrl, tgt); vm_.logical_hadamard(tgt);
                    }
                }
                break;
            case QGateKind::Custom:
                std::cerr << "[SpiralVM] Custom unitary gates not supported yet (matrix-based)\n";
                break;
            default:
                std::cerr << "[SpiralVM] Unknown gate kind\n";
        }
    };

    if (g.Kind() == QGateKind::Custom) {
        std::cerr << "[SpiralVM] Custom gates not implemented (waveform backend)\n";
    } else {
        apply_builtin(g.Kind());
    }

    // Stabilization step after every gate (tune this value)
    vm_.run_periods(1);
}

std::shared_ptr<QState> QStateSpiral::MergeDisjoint(const QState& other) const {
    // Improved MergeDisjoint taking advantage of the high logical density support.
    // Since crosstalk is negligible below 121 logical qubits per physical site,
    // we can pack the new logical qubits quite densely without immediate problems.

    auto copy = std::make_shared<QStateSpiral>(*this);

    const auto* other_spiral = dynamic_cast<const QStateSpiral*>(&other);
    if (!other_spiral) {
        std::cerr << "[SpiralVM] MergeDisjoint: other is not QStateSpiral\n";
        return copy;
    }

    // Place new qubits with moderate spacing (can be made denser if needed)
    uint32_t current_max_x = 0;
    for (const auto& q : vm_.logical_qubits) {
        if (q.center_x > current_max_x) current_max_x = q.center_x;
    }

    uint32_t offset_x = current_max_x + 4;   // relatively tight packing is acceptable

    for (uint32_t i = 0; i < other_spiral->NumQubits(); ++i) {
        copy->AddQubit(offset_x + (i % 8)*3, 10 + (i / 8)*3);  // 2D packing
    }

    return copy;
}

void QStateSpiral::ApproximateCollapse(uint32_t qid, int outcome) {
    if (outcome == 1) {
        vm_.apply_local_rotation(qid, M_PI, 0.0);   // bias toward |1>
    }
    vm_.run_periods(2);
}

MeasureResult QStateSpiral::Measure(Index target, Basis basis) {
    if (basis != Basis::Z) {
        std::cerr << "[SpiralVM] Only Z-basis measurement is currently supported\n";
    }

    uint32_t q = static_cast<uint32_t>(target);
    double z_exp = vm_.measure_logical_Z(q);

    // Probabilistic outcome based on expectation value
    double p0 = (1.0 + z_exp) / 2.0;
    int outcome = (vm_.dist(vm_.rng) < p0) ? 0 : 1;

    ApproximateCollapse(q, outcome);

    MeasureResult res;
    res.outcome = outcome;

    // For many networking protocols a shallow copy is acceptable
    res.measured = std::make_shared<QStateSpiral>(*this);
    res.survivor = std::make_shared<QStateSpiral>(*this);   // TODO: proper qubit removal + index shift

    return res;
}

void QStateSpiral::Print(std::ostream& os) const {
    PrintHeader(os, "SpiralVM (Waveform / Phase-Kick)");

    os << "  Logical qubits : " << NumQubits() << "\n";
    os << "  Use phi direct : " << (use_phi_direct_ ? "yes" : "no") << "\n";
    os << "  Current period : " << vm_.get_period() << "\n";
    os << "  Underlying state: full lattice phi (" << vm_.N << " sites)\n\n";

    os << "  Qubit  |   <Z>     |   Phase (wrapped)\n";
    os << "  -------+-----------+-------------------\n";

    for (std::size_t q = 0; q < NumQubits(); ++q) {
        double z = vm_.measure_logical_Z(static_cast<uint32_t>(q));
        double phi = vm_.get_logical_phase(static_cast<uint32_t>(q));
        os << "  " << std::setw(6) << q
           << " | " << std::fixed << std::setprecision(6) << std::setw(9) << z
           << " | " << std::setw(12) << (phi / (2.0 * M_PI)) << " * 2π\n";
    }

    // Optional: show some pairwise correlations
    if (NumQubits() >= 2) {
        os << "\n  ZZ correlations (selected):\n";
        for (std::size_t i = 0; i < std::min<std::size_t>(NumQubits(), 4); ++i) {
            for (std::size_t j = i+1; j < std::min<std::size_t>(NumQubits(), 5); ++j) {
                double zz = vm_.logical_zz_correlation(static_cast<uint32_t>(i),
                                                       static_cast<uint32_t>(j));
                os << "    <Z" << i << " Z" << j << "> = " << std::fixed << std::setprecision(4) << zz << "\n";
            }
        }
    }
}

// =============================================================================
// Factory for QStateRegistry
// =============================================================================

std::shared_ptr<QState> CreateQStateSpiral(std::size_t grid_size, bool use_phi_direct)
{
    return std::make_shared<QStateSpiral>(grid_size, use_phi_direct);
}

} // namespace q2ns