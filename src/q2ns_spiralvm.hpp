/*-----------------------------------------------------------------------------
 * Q2NS SpiralVM Backend
 * Backend-agnostic wrapper around your waveform-based SpiralVM simulator.
 *---------------------------------------------------------------------------*/

#pragma once

#include "q2ns-qstate.h"
#include "q2ns-qgate.h"
#include "spiral_vm_core.hpp"   // Your header

#include <memory>
#include <iostream>
#include <iomanip>

namespace q2ns {

/**
 * @brief Q2NS backend using SpiralVM (analog waveform / Floquet / phase-kick model)
 *
 * This backend maps Q2NS gate descriptors to your logical_* and phase-ramp functions.
 * It is particularly suited for protocols that rely on Bell states, GHZ, phase estimation,
 * and controlled-phase kicks (like your Shor demo).
 *
 * Note: The underlying quantum state is the full lattice wavefunction `phi` in SpiralVM.
 * Logical qubits are encoded as localized neighborhoods on the physical lattice.
 *
 * Limitations (by design):
 *  - Measurements are approximate (expectation → probabilistic outcome + bias collapse)
 *  - MergeDisjoint is a simple copy + spatial relocation (true disjoint tensor product is hard)
 *  - Only Z-basis measurement is supported for now
 *  - Not mathematically equivalent to QStateKet due to lattice encoding and continuous evolution
 */
class QStateSpiral : public QState {
public:
    /**
     * @param grid_size   Square lattice side length (e.g. 30)
     * @param use_phi_direct  true = fast direct phi manipulation (recommended for most protocols)
     */
    explicit QStateSpiral(std::size_t grid_size = 30, bool use_phi_direct = true);

    ~QStateSpiral() override = default;

    std::size_t NumQubits() const override;

    void Apply(const QGate& g, const std::vector<Index>& targets) override;

    std::shared_ptr<QState> MergeDisjoint(const QState& other) const override;

    MeasureResult Measure(Index target, Basis basis = Basis::Z) override;

    void Print(std::ostream& os) const override;

    int64_t AssignStreams(int64_t stream) override;

    // Convenience: add a logical qubit at a physical location
    uint32_t AddQubit(uint32_t x, uint32_t y);

    // Expose VM for debugging / custom operations (e.g. your Shor-style Ua)
    SpiralVM& GetVM() { return vm_; }
    const SpiralVM& GetVM() const { return vm_; }

private:
    SpiralVM vm_;
    bool use_phi_direct_;

    // Simple collapse helper used by Measure
    void ApproximateCollapse(uint32_t qid, int outcome);
};

/**
 * @brief Factory function to create a SpiralVM backend state.
 *
 * Used by QStateRegistry when Spiral backend is selected.
 */
std::shared_ptr<QState> CreateQStateSpiral(std::size_t grid_size = 30,
                                           bool use_phi_direct = true);

} // namespace q2ns