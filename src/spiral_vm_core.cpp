// spiral_vm_core.cpp
// Implementation for SpiralVM Option 1 (waveform multiplexing + frequency allocator + shaped envelopes)
// Requires armadillo and C++17
#include "spiral_vm_core.hpp"
#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <cassert>

using namespace arma;
using namespace std;

// ---------- Utility small helpers ----------
static inline double sq(double x){ return x*x; }

// ---------- Constructor ----------
SpiralVM::SpiralVM(int r, int c)
: rows(r), cols(c), N(r*c),
  J(0.3), h0(0.0), h1(2.5/4.4),
  omega(65536*2*M_PI), T(2*M_PI/(65536*2*M_PI)), // default values; T overwritten if omega set
  is_ang(true),
  state(2*r*c, 1, fill::zeros),
  phi(2*r*c, 1, fill::zeros),
  phi_in(2*r*c, 1, fill::zeros),
  steps(200),
  current_period(0),
  sx_gain(1900.0),
  omega_ang_base(0.0),
  rng(777),
  dist(0.0, 1.0),
  fidelities(5001, 0.0),
  fidelity_window(32, 0.0),
  waveforms(),
  drive_index(N, -1),
  virtual_frame_phase(32, 0.0),
  R((int)(sqrt((double)r*(double)c) * 0.5)), // Physical neighborhood radius around each logical qubit's center, this is the optimal
  allocated_carriers()
{
    // Keep T consistent with omega if user sets later
    T = 2.0 * M_PI / omega;

    // seed RNG
    rng.seed(777);
    // create a default global waveform so unassigned drive_index maps to 0
    Waveform global;
    global.tones.push_back(Tone(h1, omega/2.0, M_PI/4.0)); // base modulation tone (example)
    waveforms.push_back(global);
    // assign all physical sites to global waveform by default
    for (int i = 0; i < N; ++i) drive_index[i] = 0;
}

// clamp tone amplitude to avoid exploding fields
double SpiralVM::clamp_tone_amp(double a) const {
    if (std::isnan(a) || std::isinf(a)) return 0.0;
    if (a > max_tone_amp) return max_tone_amp;
    if (a < -max_tone_amp) return -max_tone_amp;
    return a;
}

// ---------- Initialize state ----------
void SpiralVM::initialize_state(const string& initial_state) {
    // uses phi member
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            int i = row * cols + col;
            if (initial_state == "neel") {
                if ((row + col) % 2 == 0) {
                    phi(i * D, 0) = 1.0;
                    phi(i * D + 1, 0) = 0.0;
                } else {
                    phi(i * D, 0) = 0.0;
                    phi(i * D + 1, 0) = 1.0;
                }
            } else if (initial_state == "polarized") {
                phi(i * D, 0) = 1.0;
                phi(i * D + 1, 0) = 0.0;
            } else if (initial_state == "disordered") {
                if (dist(rng) < 0.5) {
                    phi(i * D, 0) = 1.0;
                    phi(i * D + 1, 0) = 0.0;
                } else {
                    phi(i * D, 0) = 0.0;
                    phi(i * D + 1, 0) = 1.0;
                }
            } else {
                // default polarized
                phi(i * D, 0) = 1.0;
                phi(i * D + 1, 0) = 0.0;
            }
        }
    }

    double norm0 = sqrt(real(inner_product_cl10(phi, phi)));
    if (initial_state != "disordered") {
        phi /= norm0;
        phi_in = phi;
    } else {
        phi_in = phi / norm0;
    }
    fidelities.assign(fidelities.size(), 0.0);
    fidelities[0] = 1.0;
    current_period = 0;
    state = phi;
    cout << "[SpiralVM] Initialized (" << initial_state << "), norm0=" << norm0 << "\n";
}

// ---------- Add logical qubit + allocate waveform ----------------
int SpiralVM::allocate_waveform_for_qubit(uint32_t qid) {
    double base = freq_base;
    double spacing = freq_spacing;

    double target = base + static_cast<double>(qid) * spacing;
    bool ok = false;
    int attempt = 0;
    const int max_attempts = 1000;
    while (!ok && attempt < max_attempts) {
        ok = true;
        for (double c : allocated_carriers) {
            if (fabs(c - target) < 0.5 * spacing) { ok = false; break; }
        }
        if (!ok) target += spacing * 0.5 * (1.0 + (attempt%2?1:-1));
        attempt++;
    }
    if (attempt >= max_attempts) target = base + allocated_carriers.size() * spacing * 1.2;

    allocated_carriers.push_back(target);

    Waveform w = make_default_logical_waveform(qid);
    if (!w.tones.empty()) w.tones[0].freq = target;

    lowpass_filter_waveform(w, lowpass_cutoff);

    int wid = waveforms.size();
    waveforms.push_back(w);
    return wid;
}


// Add this implementation to spiral_vm_core.cpp (place it near allocate_waveform_for_qubit or other setup functions)

// Compilation to single physical (global multi-tone) waveform
void SpiralVM::compile_to_physical_waveform() {
    if (waveforms.empty()) { std::cerr << "[SpiralVM] No waveforms to compile\n"; return; }

    // Use the global physical_waveform member - DON'T touch logical waveforms
    physical_waveform.tones.clear();

    constexpr double freq_eps  = 1e-9;
    constexpr double phase_eps = 1e-9;
    using Key = std::tuple<double,double,int>;

    auto cmp = [freq_eps, phase_eps](const Key& a, const Key& b) {
        if (std::abs(std::get<0>(a) - std::get<0>(b)) > freq_eps)
            return std::get<0>(a) < std::get<0>(b);
        if (std::abs(std::get<1>(a) - std::get<1>(b)) > phase_eps)
            return std::get<1>(a) < std::get<1>(b);
        return std::get<2>(a) < std::get<2>(b);
    };

    std::map<Key, double, decltype(cmp)> tone_map(cmp);

    // Merge ALL waveforms → single physical (logical waveforms preserved)
    for (const auto& wf : waveforms) {
        for (const auto& tn : wf.tones) {
            tone_map[{tn.freq, tn.phase, tn.logical_id}] += tn.amp;
        }
    }

    for (const auto& entry : tone_map) {
        double freq  = std::get<0>(entry.first);
        double phase = std::get<1>(entry.first);
        int    qid   = std::get<2>(entry.first);
        double amp   = clamp_tone_amp(entry.second);
        physical_waveform.tones.emplace_back(amp, freq, phase, 0.0, 1.0, qid);
    }

    std::fill(drive_index.begin(), drive_index.end(), -2);  // flag: use physical_waveform

    cout << "[SpiralVM] Compiled " << waveforms.size() 
         << " logical → 1 physical waveform (" << physical_waveform.tones.size() << " tones)\n";
}



//For Which frequency addresses which logical qubit
void SpiralVM::dump_frequency_mapping(const std::string& fname) const {
    if (allocated_carriers.empty()) return;

    std::ofstream fout(fname);
    fout << "{\n  \"mapping\": [\n";
    for (size_t qid = 0; qid < allocated_carriers.size(); ++qid) {
        if (qid > 0) fout << ",\n";
        fout << "    {\"qid\": " << qid << ", \"carrier_freq\": " << allocated_carriers[qid] << "}";
    }
    fout << "\n  ]\n}\n";
    fout.close();
    cout << "[SpiralVM] Dumped frequency → logical qubit mapping to " << fname << "\n";
}


Waveform SpiralVM::make_default_logical_waveform(uint32_t qid) {
    Waveform w;

    // --- 1. Single dominant tone per qubit ---
    double carrier = freq_base + qid * freq_spacing;  
    double amp = clamp_tone_amp(h1);        // main amplitude
    double phase = (qid % 2 ? M_PI/3 : -M_PI/3); // staggered initial phase per qubit
    w.tones.push_back(Tone(amp, carrier, phase, 0.0, 1.0, qid));

    // --- 2. Tiny orthogonal “helper” tone for crosstalk cancellation ---
    double helper_freq = carrier + 0.1*freq_spacing; 
    double helper_amp = clamp_tone_amp(0.02 * h1); // tiny amplitude
    double helper_phase = M_PI/4 + 0.1*qid;        // slight phase shift
    w.tones.push_back(Tone(helper_amp, helper_freq, helper_phase, 0.0, 1.0, qid));

    // --- 3. Low-pass / spectral cleanup ---
    lowpass_filter_waveform(w, lowpass_cutoff);

    return w;
}


// create a CZ-style waveform that briefly introduces correlated phase near both qubits
Waveform SpiralVM::make_cz_waveform(uint32_t qid1, uint32_t qid2,
                                    double zz_strength,      // desired max ZZ coupling (rad/period)
                                    double gate_duration_frac) {
    Waveform w;

    double ω1 = freq_base + qid1 * freq_spacing;
    double ω2 = freq_base + qid2 * freq_spacing;
    double δ = std::abs(ω1 - ω2);           // detuning / beat freq

    // Main carriers — keep them on (or slightly boosted) to maintain dressing / AC Stark
    double carrier_amp = clamp_tone_amp(0.8 * h1);   // nominal drive level
    w.tones.emplace_back(carrier_amp, ω1, 0.0, 0.0, 1.0, qid1);
    w.tones.emplace_back(carrier_amp, ω2, 0.0, 0.0, 1.0, qid2);

    // Interaction tone — most physical choice is usually at |ω1 - ω2| (two-photon process)
    // Amplitude tuned so that integrated ZZ coupling ≈ zz_strength over gate time
    double interaction_amp = clamp_tone_amp(zz_strength / (gate_duration_frac * M_PI)); // rough scaling
    double interaction_freq = δ;  // or (ω1 + ω2)/2 for Raman-like, but δ is more common for ZZ

    w.tones.emplace_back(interaction_amp, interaction_freq, 0.0, 0.0, gate_duration_frac, -1);

    // Optional: small AC Stark correction tone at average frequency if needed
    // w.tones.emplace_back(0.05 * carrier_amp, (ω1 + ω2)/2, M_PI, 0.0, 1.0, -1);

    lowpass_filter_waveform(w, lowpass_cutoff);

    return w;
}

// crude low-pass filter: reduce amplitude of high-frequency tones relative to cutoff factor
void SpiralVM::lowpass_filter_waveform(Waveform &w, double cutoff_factor) {
    // cutoff_factor in (0,1]: relative attenuation for components beyond cutoff location.
    // We implement a soft attenuation proportional to freq: high freqs reduced.
    double maxfreq = 0.0;
    for (auto &t : w.tones) maxfreq = max(maxfreq, fabs(t.freq));
    if (maxfreq < 1e-12) return;
    double cutoff = cutoff_factor * maxfreq;
    for (auto &t : w.tones) {
        double f = fabs(t.freq);
        if (f > cutoff) {
            double factor = cutoff / f;
            t.amp *= factor;
        }
        // safety clamp
        t.amp = clamp_tone_amp(t.amp);
    }
}

// evaluate waveform with optional local envelope-handling (period fraction ignored for now)
double SpiralVM::eval_waveform_with_envelope(const Waveform &w, double t, double local_phase) const {
    double I_total = 0.0, Q_total = 0.0;
    
    for (const auto &tn : w.tones) {
        double arg = tn.freq * t + tn.phase + local_phase;
        // IQ MODULATION → full XY control!
        I_total += tn.I_component * cos(arg) - tn.Q_component * sin(arg);
        Q_total += tn.I_component * sin(arg) + tn.Q_component * cos(arg);
    }
    
    return sqrt(I_total*I_total + Q_total*Q_total);  // Vector Rabi drive
}


// ---------- Add logical qubit ----------
uint32_t SpiralVM::add_qubit(uint32_t x, uint32_t y) {
    LogicalQubit q;
    q.center_x = x;
    q.center_y = y;
    q.base_phase = 0.0;
    logical_qubits.push_back(q);
    uint32_t qid = logical_qubits.size() - 1;

    // allocate waveform and attach to physical sites in radius R
    int wf_id = allocate_waveform_for_qubit(qid);
    logical_qubits[qid].waveform_id = wf_id;

    for (int row = (int)q.center_y - R; row <= (int)q.center_y + R; ++row) {
        if (row < 0 || row >= rows) continue;
        for (int col = (int)q.center_x - R; col <= (int)q.center_x + R; ++col) {
            if (col < 0 || col >= cols) continue;
            int phys_idx = row * cols + col;
            drive_index[phys_idx] = wf_id;
            // phys_to_logicals bookkeeping
            if ((int)phys_to_logicals.size() < N) phys_to_logicals.resize(N);
            phys_to_logicals[phys_idx].push_back(qid);
        }
    }

    cout << "[SpiralVM] Added logical qubit " << qid << " @("<<x<<","<<y<<"), wf="<<wf_id<<"\n";
    return qid;
}



// ---------- Sample one waveform over one full Floquet period ----------
// Sample one waveform over one full Floquet period
void SpiralVM::sample_waveform(const Waveform& w, double t_start, double dt,
                               arma::vec& times,
                               arma::cx_vec& iq,
                               arma::vec& amps,
                               arma::vec& phases) const {
    int n_samples = static_cast<int>(std::ceil(T / dt));
    times.set_size(n_samples);
    iq.set_size(n_samples);
    amps.set_size(n_samples);
    phases.set_size(n_samples);

    for (int k = 0; k < n_samples; ++k) {
        double t = t_start + k * dt;
        double real_part = 0.0;
        double imag_part = 0.0;

        // Sum contributions from all tones
        for (const auto& tn : w.tones) {
            double arg = tn.freq * t + tn.phase;
            real_part += tn.amp * std::cos(arg);
            imag_part += tn.amp * std::sin(arg);
        }

        times(k) = t;
        iq(k) = cx_double(real_part, imag_part);

        double amp = std::sqrt(real_part*real_part + imag_part*imag_part);
        double phase = (amp > 1e-12) ? std::atan2(imag_part, real_part) : 0.0;

        amps(k)   = amp;
        phases(k) = phase;
    }
}

// ---------- Main dump function ----------
void SpiralVM::dump_waveforms(const std::string& format,
                              const std::string& prefix,
                              int period) const {

    // Only dump the physical global (id 0)
    size_t wid = 0;
    std::string fname_base = prefix + "global_period" + (period < 0 ? "_latest" : std::to_string(period));

    if (format == "csv") {
        //for (size_t wid = 0; wid < waveforms.size(); ++wid) {
            std::string fname = fname_base + "_wf" + std::to_string(wid) + ".csv";
            std::ofstream fout(fname);
            if (!fout) {
                std::cerr << "Failed to open " << fname << "\n";
                return;//continue;
            }

            // Metadata header
            fout << "# SpiralVM Waveform Dump\n";
            fout << "# Waveform ID: " << wid << "\n";
            fout << "# Period: " << (period < 0 ? current_period : period) << "\n";
            fout << "# T (s): " << T << "\n";
            fout << "# Sample dt (s): " << T/1000.0 << " (1000 pts per period)\n";
            fout << "# time_s,I,Q,amp,phase_rad\n";

            arma::vec times;
            arma::cx_vec iq;
            arma::vec amps;
            arma::vec phases;

            double dt = T / 1000.0;  // 1000 samples per period — adjust as needed
            sample_waveform(waveforms[wid], 0.0, dt, times, iq, amps, phases);

            for (size_t k = 0; k < times.size(); ++k) {
                fout << times(k) << "," << iq(k).real() << "," << iq(k).imag() << "," << amps(k) << "," << phases(k) << "\n";
            }
            fout.close();
            //std::cout << "Dumped CSV: " << fname << "\n";
        //}
    } else if (format == "json") {
        // Similar but JSON structure — more verbose, good for metadata-heavy use
        std::string fname = fname_base + ".json";
        std::ofstream fout(fname);
        if (!fout) {
            std::cerr << "Failed to open " << fname << "\n";
            return;
        }

        fout << "{\n";
        fout << "  \"spiral_vm_version\": \"0.1\",\n";
        fout << "  \"period\": " << (period < 0 ? current_period : period) << ",\n";
        fout << "  \"T\": " << T << ",\n";
        fout << "  \"physical_waveforms\": [\n";

        //for (size_t wid = 0; wid < waveforms.size(); ++wid) {
            if (wid > 0) fout << ",\n";
            fout << "    {\n";
            fout << "      \"id\": " << wid << ",\n";
            fout << "      \"tones\": [\n";
            for (size_t t = 0; t < physical_waveform.tones.size(); ++t) {
                if (t > 0) fout << ",\n";
                const auto& tn = physical_waveform.tones[t];
                fout << "        {\"amp\": " << tn.amp
                     << ", \"freq\": " << tn.freq
                     << ", \"phase\": " << tn.phase << "}";
            }
            fout << "\n      ],\n";

            // Sampled IQ
            arma::vec times;
            arma::cx_vec iq;
            arma::vec amps;
            arma::vec phases;
            double dt = T / 1000.0;
            sample_waveform(physical_waveform, 0.0, dt, times, iq, amps, phases);
            fout << "      \"samples\": [\n";
            for (size_t k = 0; k < times.size(); ++k) {
                if (k > 0) fout << ",\n";
                fout << "        {\"t\": " << times(k)
                     << ", \"I\": " << iq(k).real()
                     << ", \"Q\": " << iq(k).imag() << "}";
            }
            fout << "\n      ]\n";
            fout << "    }";
        //}
        fout << "\n  ]\n}\n";
        fout.close();
        //std::cout << "Dumped JSON: " << fname << "\n";
    } else {
        std::cerr << "Unsupported format: " << format << " (use csv or json)\n";
    }
}

void SpiralVM::apply_gate(const Gate& g, double period_time) {
    // period_time is typically in units of periods (e.g. 1.0 = one full Floquet period)
    // We use it to scale duration / strength of control pulses

    switch (g.type) {

        case Gate::X: {
            // Goal: logical X on target qubit → should become a π-pulse on that logical qubit only
            // Current: global_pi_pulse() flips EVERY physical site → we want to avoid this
            logical_x_pulse(g.target, period_time);
            break;
        }

        case Gate::Z: {
            // Logical Z rotation by angle g.angle on the target qubit
            // We already have the machinery → just phase-ramp the logical qubit's carrier tone
            logical_z_rotation(g.target, g.angle);   // Phase π/2
            // Note: if angle is large, you may want to split into multiple smaller ramps
            break;
        }

        case Gate::PHASE: {
            // Same as Z — single-qubit phase gate
            logical_phase_ramp(g.target, g.angle, 1);
            break;
        }

        case Gate::H: {
            // Logical Hadamard on target
            // Standard decomposition: Z^(π/2) → X → Z^(-π/2)
            // We already have logical_hadamard doing exactly this — but it currently calls global_pi_pulse()
            // → we replace the global X with logical_x_pulse
            logical_hadamard(g.target);   // ← this function will be updated below
            break;
        }

        case Gate::T: {
            // T gate = Z^(π/4) on |1⟩ (phase gate by π/4)
            logical_phase_ramp(g.target, M_PI / 4.0, 1);  // single period, or increase steps if needed
            break;
        }

        case Gate::CZ: {
            // Control-Z between control and target
            // Current: global pi pulse + conditional phase kick → both problematic
            // Desired: controlled phase between the two logical qubits without global flip
            //
            logical_cz(g.control, g.target);
            break;
        }

        case Gate::CNOT: {
            // Standard decomposition: H target → CZ → H target
            // We keep the structure but rely on updated logical_hadamard and logical_controlled_phase
            logical_hadamard(g.target);
            logical_controlled_phase(g.control, g.target, M_PI, period_time * 0.8); // slightly shorter to avoid over-rotation
            logical_hadamard(g.target);
            break;
        }

        case Gate::RX: {
            // Variable X rotation: scale pulse strength by angle
            double x_strength = LOGICAL_X_AMPLITUDE * (g.angle / M_PI);
            // Temporarily override for partial flip
            logical_x_pulse(g.target, 1.0 * (g.angle / M_PI));  // Scale duration
            break;
        }

        case Gate::RY:
        case Gate::RZ: {
            // Arbitrary single-qubit rotations — to be implemented properly later
            // Placeholder strategy:
            //   RX(θ) → X^(θ/π) pulse (variable strength or duration π-pulse)
            //   RZ(θ) → phase ramp by θ
            //   RY(θ) → composite of RX + RZ or different axis modulation
            //
            // For now: minimal implementation (only RZ is trivial)
            if (g.type == Gate::RZ) {
                logical_phase_ramp(g.target, g.angle, 1);
            } else {
                // TODO: implement logical RX / RY via amplitude/phase modulated pulses on the logical carrier
                //       Example path: temporarily increase tone amplitude for duration corresponding to rotation angle
                std::cout << "[SpiralVM] " << g.type << "(θ=" << g.angle << ") on q" << g.target
                          << " not yet implemented at logical level\n";
            }
            break;
        }

        case Gate::MEASURE: {
            // Measurement in Z basis — in waveform simulation this is tricky
            // Options:
            //   1. Classical readout simulation: estimate <Z> from waveform phase accumulation or neighborhood populations
            //   2. Full projective measurement (collapse) — hard in waveform picture without feedback
            //
            // Pragmatic choice for now: just log / return the logical Z value (non-destructive)
            double z = measure_logical_Z(g.target);
            std::cout << "[SpiralVM] Measured logical Z on q" << g.target << " = " << z << "\n";
            // TODO (future): implement actual destructive measurement + state preparation if in circuit mode
            break;
        }

        default:
            std::cout << "[SpiralVM] Gate type not recognized or not implemented: " << g.type << "\n";
            break;
    }

    // Optional: force a small evolution step after every gate so waveform changes take effect
    run_periods(1);   // ← you may want this depending on your timing model
}

void SpiralVM::virtual_phase_gate(uint32_t qid, double angle) {
    if (qid >= virtual_frame_phase.size()) 
        virtual_frame_phase.resize(qid + 1, 0.0);
    virtual_frame_phase[qid] += angle;
    std::cout << "[SpiralVM] Virtual phase gate on q" << qid << ": φ_frame = " << virtual_frame_phase[qid] << "\n";
    // NO hardware interaction, NO decoherence!
}

// Controlled phase gate — the most important two-qubit gate
// Implements a hardware-like entangling gate by temporarily applying
// an interaction waveform (carriers + beat/interaction tone) to both qubit neighborhoods.
// No mid-circuit measurement or feedforward — entanglement arises from the drive-mediated interaction.
void SpiralVM::logical_controlled_phase(uint32_t control, uint32_t target,
                                        double max_angle, double duration_periods) {
    if (control >= logical_qubits.size() || target >= logical_qubits.size()) {
        std::cout << "[SpiralVM] Invalid qubit IDs for controlled phase: control=" 
                  << control << ", target=" << target << "\n";
        return;
    }

    // Convert max_angle (desired conditional phase when both in |11>) into waveform strength
    // Rough heuristic: strength ≈ max_angle / (duration * some calibration factor)
    // Tune the 2.0 factor based on your sim (run test gates and measure acquired phase)
    double zz_strength = LOGICAL_X_AMPLITUDE * 2.5;//; * cos(2.0 * current_period * M_PI/(512*T)) / 4;// * max_angle / (duration_periods * 2.0);  // rad per period scaling

    // Use the full waveform-based CZ implementation
    apply_phase_kick_between_full(control, target, zz_strength, duration_periods);

    std::cout << "[SpiralVM] Applied controlled phase (CZ-like) between control=" 
              << control << " and target=" << target 
              << " with max_angle=" << max_angle 
              << " over " << duration_periods << " periods\n";
}

void SpiralVM::logical_cz(uint32_t control, uint32_t target) {
    logical_controlled_phase(control, target, M_PI, 0.8);  // ~π conditional phase, slightly shorter duration
}

double SpiralVM::get_logical_phase_frame_corrected(uint32_t qid) {
    double raw = get_logical_phase(qid);
    // subtract the global drive phase
    return raw - drive_phase;
}

int SpiralVM::find_waveform_index_for_qubit(uint32_t qid) {
    if (qid >= logical_qubits.size()) return -1;

    // search all waveforms
    for (size_t w = 0; w < waveforms.size(); ++w) {
        for (const auto &tn : waveforms[w].tones) {
            if (tn.logical_id == qid) return static_cast<int>(w);
        }
    }

    // not found
    return -1;
}


// ---------- Floquet runner / logging ----------
void SpiralVM::run_floquet(int N_max, const string& initial_state) {
    // simple wrapper that writes CSV file similar to your earlier code
    ofstream fout;
    stringstream fname;
    fname << "dtc_floquet_with_" << initial_state << "_state_" << (is_ang ? "_spiral_" : "_no_spiral_") << (int)sx_gain << "_omega_" << (int)omega << ".txt";
    fout.open(fname.str());
    fout << "Period,Fidelity,Stabilizer,Energy,sx_energy,zz_energy,Delta_F,ht_eff_end,sx_avg\n";

    double delta_F = 0.0;
    // initial energies:
    sp_cx_mat H_init = hamiltonian_cl10_90_spiral_twist(J, 0, 0);
    cx_mat Hphi_init = mat_vec_mult_cl10(H_init, phi);
    double energy_init = real(inner_product_cl10(phi, Hphi_init)) + compute_zz_energy_edgeaware(phi, J, 0, 0);
    double zz_energy_init = compute_zz_energy_edgeaware(phi, J, 0, 0);
    double sx_energy_init = 0.0;
    for (int i = 0; i < N; i++) {
        sx_energy_init -= h1 * 2.0 * real(phi(i * D, 0) * conj(phi(i * D + 1, 0)));
    }
    fout << "0,1.0," << compute_avg_stabilizer(phi) << "," << energy_init << "," << sx_energy_init << "," << zz_energy_init << "," << h1 << ",0,0\n";

    for (int n = 0; n < N_max; n++) {
        step_period(n, delta_F);

        double sx_energy = 0.0;
        for (int i = 0; i < N; i++) {
            sx_energy -= 2.0 * h1 * real(phi(i * D, 0) * conj(phi(i * D + 1, 0)));
        }
        double zz_energy = is_ang ? compute_zz_energy_edgeaware(phi, J, omega_ang_end(n), n + 1, true) : compute_zz_energy_edgeaware(phi, J, 0, 0);
        double energy = sx_energy + zz_energy;

        fout << n + 1 << "," << fidelities[n + 1] << "," << compute_avg_stabilizer(phi) << "," << energy << "," << sx_energy << "," << zz_energy << "," << delta_F << "," << h_effective_end(n) << "," << sx_avg(n) << "\n";
        cout << "Period " << n + 1 << ": Fidelity = " << fidelities[n + 1] << ", Energy = " << energy << "\n";
    }
    fout.close();
}

void SpiralVM::step_period(int n, double &delta_F) {
    double dt = T / steps;
    cx_mat phi_new = phi;

    // ---- Frequency scattering parameters ----
    double f_min = 1.0e6;    // minimal frequency (Hz)
    double f_max = 1.0e9;    // maximal frequency (Hz)
    std::uniform_real_distribution<double> dist(f_min, f_max);

    std::vector<double> qubit_freqs(N, 0.0);
    for (int i = 0; i < N; ++i) {
        qubit_freqs[i] = dist(rng);   // assign a random frequency per qubit
    }

    // ---- RK4 integration per site ----
    for (int k = 0; k < steps; k++) {
        double t = n * T + k * dt;

        // local field per qubit
        std::vector<double> local_hx(N, 0.0);
        for (int i = 0; i < N; ++i) {
            int wid = drive_index[i];

            // compute physical coordinates
            int x = i % cols;
            int y = i / cols;

            if (wid == -2) {
                // use physical_waveform (flag set by compile_to_physical_waveform)
                double val = eval_waveform_with_envelope(physical_waveform, t, drive_phase);
                local_hx[i] = h0 + val;
            } else if (wid < 0 || wid >= (int)waveforms.size()) {
                // simple oscillation using scattered frequency
                local_hx[i] = h0 + h1 * cos(qubit_freqs[i] * t + M_PI/4.0 + drive_phase);
            } else {
                // logical waveform
                double val = eval_waveform_with_envelope(waveforms[wid], t, drive_phase);
                local_hx[i] = h0 + val;
            }

        }

        // angular modulation
        double angular_freq_quasi = is_ang ? (sin(omega * M_PI * t / T) + sin(2*omega * M_PI * t / T)) : 0.0;
        double omega_ang = omega_ang_base + angular_freq_quasi;

        // Hamiltonian
        sp_cx_mat H_sx = hamiltonian_cl10_90_spiral_twist_inhomogeneous(J, local_hx, omega_ang);
        cx_mat Hzz_phi = compute_zz_energy_vector_edgeaware(phi_new, J, omega_ang, n, is_ang);

        // RK4 steps
        cx_mat k1 = mat_vec_mult_cl10(H_sx, phi_new) + Hzz_phi;
        cx_mat phi_temp = phi_new + (-cx_double(0,1) * dt/2.0) * k1;

        cx_mat k2 = mat_vec_mult_cl10(H_sx, phi_temp) + compute_zz_energy_vector_edgeaware(phi_temp, J, omega_ang, n, is_ang);
        phi_temp = phi_new + (-cx_double(0,1) * dt/2.0) * k2;
        cx_mat k3 = mat_vec_mult_cl10(H_sx, phi_temp) + compute_zz_energy_vector_edgeaware(phi_temp, J, omega_ang, n, is_ang);
        phi_temp = phi_new + (-cx_double(0,1) * dt) * k3;
        cx_mat k4 = mat_vec_mult_cl10(H_sx, phi_temp) + compute_zz_energy_vector_edgeaware(phi_temp, J, omega_ang, n, is_ang);

        phi_new += (-cx_double(0,1) * dt/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);

        double normrk4 = sqrt(real(inner_product_cl10(phi_new, phi_new)));
        if (normrk4 > 0) phi_new /= normrk4;
    }

    phi = phi_new;
    state = phi;

    // fidelity
    double fid = fabs(inner_product_cl10(phi_in, phi));
    if ((int)fidelities.size() <= current_period+1) fidelities.resize(current_period+2, 0.0);
    fidelities[current_period+1] = fid;
    //std::cout << "period: " << current_period+1 << " fidelity: " << fid << std::endl;
    //std::cout << std::flush;
    delta_F = 1.0 - fid;


    // WITH:
    if (auto_compile_enabled) {
        compile_to_physical_waveform();                   // always merge
    }

    std::string fname_base = "global_" + std::to_string(current_period);
    dump_waveforms("csv", fname_base, -1); // always dump the physical one
    std::string fname_base_h = "h_eff_" + std::to_string(current_period);
    dump_h_eff(fname_base_h, -1);

    current_period++;
}

int SpiralVM::get_total_logical_qubits(){
    return logical_qubits.size();
}

// ---------- run N periods ----------------
void SpiralVM::run_periods(uint32_t N_periods) {
    double dummy_deltaF = 0.0;
    for (uint32_t i = 0; i < N_periods; i++) {
        step_period(current_period, dummy_deltaF);
    }
}

// ---------- logical measurements and gates ----------


// logical Hadamard on single qubit
void SpiralVM::logical_hadamard(uint32_t qid) {
    int wid = logical_qubits[qid].waveform_id;
    
    // INLINE SEARCH - no new function needed
    size_t carrier_idx = 0;
    for (size_t i = 0; i < waveforms[wid].tones.size(); ++i) {
        if (waveforms[wid].tones[i].logical_id == static_cast<int>(qid)) {
            carrier_idx = i;
            break;
        }
    }
    
    double h_amp = h1 / sqrt(2.0);
    waveforms[wid].tones[carrier_idx].set_iq(h_amp, h_amp);
    
    compile_to_physical_waveform();
    run_periods(1.0);
    
    waveforms[wid].tones[carrier_idx].set_iq(h1, 0.0);
    compile_to_physical_waveform();
}



double SpiralVM::current_orbit_phase(uint32_t qid) const {
    if (qid >= logical_qubits.size()) return 0.0;

    int wid = logical_qubits[qid].waveform_id;
    if (wid < 0 || wid >= (int)waveforms.size()) return 0.0;

    // Find the main carrier tone for this qubit
    for (const auto& tn : waveforms[wid].tones) {
        if (tn.logical_id == static_cast<int>(qid) &&
            std::abs(tn.freq - allocated_carriers[qid]) < 1e-6) {

            // Base phase of the tone
            double tone_phase = tn.phase;

            // Add accumulated phase from frequency over time
            double accumulated = tn.freq * (current_period * T);

            // Subtract global drive phase (so it's relative to the global modulation frame)
            double relative_phase = tone_phase + accumulated - drive_phase;

            // Normalize to [0, 2π)
            return fmod(relative_phase, 2.0 * M_PI);
        }
    }

    // Fallback: use global drive phase only
    return fmod(-drive_phase, 2.0 * M_PI);
}


// Applies a logical X (pi rotation) to the specified qubit by temporarily boosting the amplitude
// of its main carrier tone to induce a strong transverse field pulse over a calibrated duration.
// This evolves the state naturally via run_periods without touching phi directly.
void SpiralVM::logical_x_pulse(uint32_t qid, double duration_periods) {
    if (qid >= logical_qubits.size()) {
        std::cout << "[SpiralVM] Invalid qubit id for logical X: " << qid << "\n";
        return;
    }

    int wid = logical_qubits[qid].waveform_id;
    if (wid < 0 || wid >= (int)waveforms.size()) {
        std::cout << "[SpiralVM] No waveform found for qubit " << qid << "\n";
        return;
    }

    std::cout << "[DEBUG] allocated_carriers["<<qid<<"]=" << allocated_carriers[qid] << "\n";
    for (size_t i=0; i<waveforms[wid].tones.size(); i++) {
        std::cout << "[DEBUG] tone["<<i<<"] freq=" << waveforms[wid].tones[i].freq 
                  << " logical_id=" << waveforms[wid].tones[i].logical_id << "\n";
    }


    // Find the main carrier tone for this qubit (assuming first tone is the dominant carrier)
    size_t carrier_idx = -1;
    for (size_t i = 0; i < waveforms[wid].tones.size(); ++i) {
        if (waveforms[wid].tones[i].logical_id == (int)qid &&
            std::abs(waveforms[wid].tones[i].freq - allocated_carriers[qid]) < 1e-6) {
            carrier_idx = i;
            break;
        }
    }
    if (carrier_idx == static_cast<size_t>(-1)) {
        std::cout << "[SpiralVM] No carrier tone found for qubit " << qid << "\n";
        return;
    }

    // SAVE ALL STATE
    std::vector<int> saved_drive_index = drive_index;
    double original_amp = waveforms[wid].tones[carrier_idx].amp;
    

    // BOOST + ORBIT ALIGN (1 extra line!)
    double quasi = is_ang ? (sin(omega * current_period) + sin(2*omega * current_period)) : 0.0;
    waveforms[wid].tones[carrier_idx].amp *= LOGICAL_X_AMPLITUDE * cos(2.0 * current_period * M_PI/(512*T));
    waveforms[wid].tones[carrier_idx].phase += M_PI * (quasi > 0 ? 1.0 : -1.0);
    //waveforms[wid].tones[carrier_idx].amp = clamp_tone_amp(waveforms[wid].tones[carrier_idx].amp);
    std::cout << "[DEBUG] waveforms[wid].tones[carrier_idx].phase=" << waveforms[wid].tones[carrier_idx].phase << "\n";

    compile_to_physical_waveform();
    
    std::cout << "[DEBUG] physical[0].amp=" << physical_waveform.tones[0].amp << "\n";
    std::cout << "[DEBUG] physical[0].phase=" << physical_waveform.tones[0].phase << "\n";
    
    // PULSE WITH NO SIDE EFFECTS
    bool saved_auto = auto_compile_enabled;
    auto_compile_enabled = false;
    std::vector<int> pulse_drive_index = drive_index;  // capture -2's
    
    uint32_t steps = std::ceil(duration_periods);
    for(uint32_t i = 0; i < steps; i++) {
        double dummy_F = 0.0;
        drive_index = pulse_drive_index;  // FORCE -2 every step!
        step_period(current_period + i, dummy_F);  
    }
    
    // CLEANUP
    drive_index = saved_drive_index;
    waveforms[wid].tones[carrier_idx].amp = original_amp;
    auto_compile_enabled = saved_auto;
    compile_to_physical_waveform();

}

double SpiralVM::measure_logical_Z(uint32_t qid) const {
    if (qid >= logical_qubits.size()) return 0.0;
    const LogicalQubit &q = logical_qubits[qid];

    double staggered = 0.0;
    int count = 0;

    for (int row = q.center_y - R; row <= q.center_y + R; ++row) {
        if (row < 0 || row >= rows) continue;
        for (int col = q.center_x - R; col <= q.center_x + R; ++col) {
            if (col < 0 || col >= cols) continue;
            int i = row * cols + col;
            double sz = N * (std::norm(phi(i*D + 0,0)) - std::norm(phi(i*D + 1,0)));  // Add N factor
            staggered += ((row + col) % 2 == 0) ? sz : -sz;
            ++count;
        }
    }

    return (count > 0) ? staggered / (double)count : 0.0;
}

double SpiralVM::measure_logical_Z_frame_corrected(uint32_t qid) const {
    double physical_z = measure_logical_Z(qid);
    
    // Virtual frame compensation: rotate measured value
    if (qid < virtual_frame_phase.size()) {
        double frame = virtual_frame_phase[qid];
        return physical_z * std::cos(frame) + measure_logical_X(qid) * std::sin(frame);
        // Perfect frame rotation in measurement basis
    }
    return physical_z;
}


double SpiralVM::measure_logical_X(uint32_t qid) const {
    if (qid >= logical_qubits.size()) return 0.0;
    const LogicalQubit &q = logical_qubits[qid];
    double re = 0.0;
    int count = 0;
    
    for (int row = q.center_y - R; row <= q.center_y + R; ++row) {
        if (row < 0 || row >= rows) continue;
        for (int col = q.center_x - R; col <= q.center_x + R; ++col) {
            if (col < 0 || col >= cols) continue;
            int i = row * cols + col;
            
            // **CORRECT ⟨X⟩ = Re[ψ₀* ψ₁ + ψ₁* ψ₀]**
            cx_double psi0 = phi(i*D + 0, 0);
            cx_double psi1 = phi(i*D + 1, 0);
            re += std::real( std::conj(psi0)*psi1 + std::conj(psi1)*psi0 );
            
            // Staggered sign for Néel state
            if ((row + col) % 2 == 1) re = -re;
            ++count;
        }
    }
    return (count > 0) ? re / count : 0.0;
}


double SpiralVM::measure_logical_Y(uint32_t qid) const {
    if (qid >= logical_qubits.size()) return 0.0;
    const LogicalQubit &q = logical_qubits[qid];
    double sum = 0.0;
    int count = 0;
    for (int row = q.center_y - R; row <= q.center_y + R; ++row) {
        if (row < 0 || row >= rows) continue;
        for (int col = q.center_x - R; col <= q.center_x + R; ++col) {
            if (col < 0 || col >= cols) continue;
            int i = row * cols + col;
            // ⟨Y⟩ ≈ 2 Im(⟨0|1⟩) per site
            double sy = 2.0 * imag(phi(i*D + 0, 0) * conj(phi(i*D + 1, 0)));
            sum += ((row + col) % 2 == 0) ? sy : -sy;
            ++count;
        }
    }
    return (count > 0) ? sum / count : 0.0;
}

void SpiralVM::logical_phase_ramp(int target_qid, double slope, int steps) {
    for (int k = 0; k < steps; ++k) {
        for (auto &w : waveforms) {
            for (auto &tn : w.tones) {
                if (tn.logical_id == target_qid) {
                    tn.phase += slope;
                }
            }
        }

        compile_to_physical_waveform();  // rebuild merged drive
        run_periods(1);                  // let it act on φ
    }
}


double SpiralVM::get_logical_phase(uint32_t qid) {
    if (qid >= logical_qubits.size()) return 0.0;
    const LogicalQubit &q = logical_qubits[qid];

    // pick the center site of the qubit
    int row = q.center_y;
    int col = q.center_x;
    if (row < 0 || row >= rows || col < 0 || col >= cols) return 0.0;

    int i = row*cols + col;
    // D offset: pick the "physical" component corresponding to the |1> amplitude
    std::complex<double> amp_center = phi(i*D + 1,0);

    return std::arg(amp_center);
}

void SpiralVM::logical_z_rotation(uint32_t qid, double angle) {
    if (qid >= logical_qubits.size()) return;
    
    int wid = logical_qubits[qid].waveform_id;
    if (wid < 0 || wid >= (int)waveforms.size()) return;
    
    // Find carrier tone (like your X pulse)
    size_t carrier_idx = -1;
    for (size_t i = 0; i < waveforms[wid].tones.size(); ++i) {
        if (waveforms[wid].tones[i].logical_id == (int)qid &&
            std::abs(waveforms[wid].tones[i].freq - allocated_carriers[qid]) < 1e-6) {
            carrier_idx = i;
            break;
        }
    }
    if (carrier_idx == static_cast<size_t>(-1)) return;

    // SAVE STATE
    std::vector<int> saved_drive_index = drive_index;
    double original_freq = waveforms[wid].tones[carrier_idx].freq;
    double original_phase = waveforms[wid].tones[carrier_idx].phase;
    
    // **Z FIELD = CARRIER DETUNING** (orthogonal to X drive)
    waveforms[wid].tones[carrier_idx].amp /=  0.389 * M_PI * angle * sin(angle);//cos(current_period * angle / T) * sin(current_period * angle / T);;// * angle;  // detuning!
    waveforms[wid].tones[carrier_idx].phase -= M_PI/2 * LOGICAL_X_AMPLITUDE;
    
    compile_to_physical_waveform();
    
    // **EVOLVE under pure Z Hamiltonian**
    bool saved_auto = auto_compile_enabled;
    auto_compile_enabled = false;
    std::vector<int> pulse_drive_index = drive_index;
    
    run_periods(1);  // 1 period Z evolution
    
    // RESTORE
    drive_index = saved_drive_index;
    waveforms[wid].tones[carrier_idx].freq = original_freq;
    waveforms[wid].tones[carrier_idx].phase = original_phase;
    auto_compile_enabled = saved_auto;
    compile_to_physical_waveform();
}



double SpiralVM::logical_zz_correlation(uint32_t qid1, uint32_t qid2) const {
    if (qid1 >= logical_qubits.size() || qid2 >= logical_qubits.size()) return 0.0;

    const LogicalQubit &qa = logical_qubits[qid1];
    const LogicalQubit &qb = logical_qubits[qid2];

    double sum_zz = 0.0;
    int count = 0;

    for (int ra = qa.center_y - R; ra <= qa.center_y + R; ++ra) {
        if (ra < 0 || ra >= rows) continue;
        for (int ca = qa.center_x - R; ca <= qa.center_x + R; ++ca) {
            if (ca < 0 || ca >= cols) continue;
            int ia = ra * cols + ca;

            double za = norm(phi(ia*D + 0, 0)) - norm(phi(ia*D + 1, 0));  // raw P0 - P1 per site

            for (int rb = qb.center_y - R; rb <= qb.center_y + R; ++rb) {
                if (rb < 0 || rb >= rows) continue;
                for (int cb = qb.center_x - R; cb <= qb.center_x + R; ++cb) {
                    if (cb < 0 || cb >= cols) continue;
                    int ib = rb * cols + cb;

                    double zb = norm(phi(ib*D + 0, 0)) - norm(phi(ib*D + 1, 0));

                    /**/
                    // Apply the same staggering as measure_logical_Z
                    double sign_a = ((ra + ca) % 2 == 0) ? 1.0 : -1.0;
                    double sign_b = ((rb + cb) % 2 == 0) ? 1.0 : -1.0;
                    /**/

                    // WITH THIS (raw correlations):for neel init?
                    //double sign_a = 1.0;
                    //double sign_b = 1.0;

                    sum_zz += (za * sign_a) * (zb * sign_b);
                    count++;
                }
            }
        }
    }

    return count > 0 ? sum_zz / count : 0.0;
}

// Preferred / physical CZ implementation
void SpiralVM::apply_phase_kick_between_full(uint32_t qid1, uint32_t qid2,
                                             double zz_strength, double duration_fraction) {
    if (qid1 >= logical_qubits.size() || qid2 >= logical_qubits.size()) return;

    // Generate CZ waveform with desired integrated ZZ
    Waveform w = make_cz_waveform(qid1, qid2, zz_strength, duration_fraction);

    int temp_wid = waveforms.size();
    waveforms.push_back(std::move(w));

    // Temporarily assign to both logical neighborhoods
    std::vector<int> old_drive_indices;
    for (int q : {qid1, qid2}) {
        auto& qb = logical_qubits[q];
        for (int r = qb.center_y - R; r <= qb.center_y + R; ++r) {
            if (r < 0 || r >= rows) continue;
            for (int c = qb.center_x - R; c <= qb.center_x + R; ++c) {
                if (c < 0 || c >= cols) continue;
                int idx = r * cols + c;
                old_drive_indices.push_back(drive_index[idx]);  // save
                drive_index[idx] = temp_wid;
            }
        }
    }

    // Let the interaction act
    run_periods(1);   // or more if duration_fraction > 1

    // Restore original drive assignments
    size_t restore_idx = 0;
    for (int q : {qid1, qid2}) {
        auto& qb = logical_qubits[q];
        for (int r = qb.center_y - R; r <= qb.center_y + R; ++r) {
            if (r < 0 || r >= rows) continue;
            for (int c = qb.center_x - R; c <= qb.center_x + R; ++c) {
                if (c < 0 || c >= cols) continue;
                int idx = r * cols + c;
                drive_index[idx] = old_drive_indices[restore_idx++];
            }
        }
    }

    // Optional: remove temp waveform to save memory
    waveforms.pop_back();
}

// ---------- Hamiltonian builder: keep old signature intact and add new inhomogeneous builder ----------
void SpiralVM::compute_nonzero_indices_spiral_twist(double J, double ht, int rows, int cols, int D, double omega_ang, arma::umat& locations, arma::cx_vec& values, uint& nz) {
    int N_local = rows * cols;
    const int NNZ = 2 * N_local;
    locations.set_size(2, NNZ);
    values.set_size(NNZ);
    nz = 0;
    for (int i = 0; i < N_local; i++) {
        locations(0, nz) = i * D;     locations(1, nz) = i * D + 1; values(nz) = -ht; nz++;
        locations(0, nz) = i * D + 1; locations(1, nz) = i * D;     values(nz) = -ht; nz++;
    }
}

// original (uniform) Hamiltonian
arma::sp_cx_mat SpiralVM::hamiltonian_cl10_90_spiral_twist(double J, double ht, double omega_ang) {
    int N_local = rows * cols;
    const int NNZ = 2 * N_local;
    arma::umat locations(2, NNZ);
    arma::cx_vec values(NNZ);
    uint nz;
    compute_nonzero_indices_spiral_twist(J, ht, rows, cols, D, omega_ang, locations, values, nz);
    return arma::sp_cx_mat(locations.submat(0,0,1,nz-1), values.subvec(0,nz-1), N_local*D, N_local*D);
}

// new: inhomogeneous transverse field Hamiltonian
arma::sp_cx_mat SpiralVM::hamiltonian_cl10_90_spiral_twist_inhomogeneous(double J, const std::vector<double> &local_hx, double omega_ang) {
    int N_local = rows * cols;
    // We construct sigma_x terms with per-site coefficients.
    // Each physical site i contributes -hx_i * (|0><1| + |1><0|)
    // plus we keep ZZ terms out of this routine (handled separately via compute_zz_energy_vector)
    // Build sparse matrix entries for all |0><1| and |1><0|
    int NNZ = 2 * N_local;
    arma::umat locations(2, NNZ);
    arma::cx_vec values(NNZ);
    uint nz = 0;
    for (int i = 0; i < N_local; ++i) {
        double hx = (i < (int)local_hx.size()) ? local_hx[i] : (h0 + h1 * cos(omega * (0.0/2.0) + M_PI/4.0));
        // off-diagonal entries:
        locations(0, nz) = i*D; locations(1, nz) = i*D + 1; values(nz) = -hx; nz++;
        locations(0, nz) = i*D + 1; locations(1, nz) = i*D; values(nz) = -hx; nz++;
    }
    return arma::sp_cx_mat(locations.submat(0,0,1,nz-1), values.subvec(0,nz-1), N_local*D, N_local*D);
}

// ---------- Edge-aware neighbor lookup ----------
inline int SpiralVM::get_right_neighbor(int row, int col) const {
    if (col + 1 < cols) return row * cols + (col + 1);
    if (col - 1 >= 0) return row * cols + (col - 1);
    return row * cols + col; // fallback to self
}

inline int SpiralVM::get_down_neighbor(int row, int col) const {
    if (row + 1 < rows) return (row + 1) * cols + col;
    if (row - 1 >= 0) return (row - 1) * cols + col;
    return row * cols + col; // fallback to self
}

// ---------- Edge-aware ZZ energy ----------
double SpiralVM::compute_zz_energy_edgeaware(const arma::cx_mat& phi_inp, double J, double omega_ang, double period, bool is_angf) {
    double energy = 0.0;
    double norm_factor = sqrt((double)rows * cols);

    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            int i = row * cols + col;
            int j_right = get_right_neighbor(row, col);
            int j_down  = get_down_neighbor(row, col);

            cx_double zi  = (phi_inp(i*D,0) - phi_inp(i*D + 1,0)) * norm_factor;
            cx_double zjr = (phi_inp(j_right*D,0) - phi_inp(j_right*D + 1,0)) * norm_factor;
            cx_double zjd = (phi_inp(j_down*D,0)  - phi_inp(j_down*D + 1,0)) * norm_factor;

            cx_double J_twist = is_angf ? J * cx_double(0,1) : cx_double(J,0);

            if (is_angf) {
                energy += imag(J_twist * zi * zjr);
                energy += imag(J_twist * zi * zjd);
            } else {
                energy += real(J_twist * zi * zjr);
                energy += real(J_twist * zi * zjd);
            }
        }
    }
    return energy;
}

// ---------- Edge-aware ZZ energy vector ----------
arma::cx_mat SpiralVM::compute_zz_energy_vector_edgeaware(const arma::cx_mat& phi_inp, double J, double omega_ang, double period, bool is_angf) {
    arma::cx_mat Hzz_phi = arma::zeros<arma::cx_mat>(N*D, 1);
    double norm_factor = sqrt((double)rows * cols);

    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            int i = row * cols + col;
            int j_right = get_right_neighbor(row, col);
            int j_down  = get_down_neighbor(row, col);

            cx_double zi  = (phi_inp(i*D,0) - phi_inp(i*D + 1,0)) * norm_factor;
            cx_double zjr = (phi_inp(j_right*D,0) - phi_inp(j_right*D + 1,0)) * norm_factor;
            cx_double zjd = (phi_inp(j_down*D,0)  - phi_inp(j_down*D + 1,0)) * norm_factor;

            cx_double J_twist = is_angf ? J * cx_double(0,1) : cx_double(J,0);

            if (is_angf) {
                Hzz_phi(i*D + 0,0) += imag(J_twist * (zjr + zjd));
                Hzz_phi(i*D + 1,0) -= imag(J_twist * (zjr + zjd));
            } else {
                Hzz_phi(i*D + 0,0) += real(J_twist * (zjr + zjd));
                Hzz_phi(i*D + 1,0) -= real(J_twist * (zjr + zjd));
            }
        }
    }
    return Hzz_phi;
}


// ---------- basic linear algebra helpers ----------
arma::cx_mat SpiralVM::mat_vec_mult_cl10(const arma::sp_cx_mat& H, const arma::cx_mat& phi_inp) {
    return H * phi_inp;
}

double SpiralVM::inner_product_cl10(const arma::cx_mat& phi1, const arma::cx_mat& phi2) {
    cx_double prod = as_scalar(phi1.t() * phi2);
    return real(prod);
}

double SpiralVM::sx_avg(int n) {
    double sx_sum = 0.0;
    for (int i = 0; i < N; ++i) {
        sx_sum += 2.0 * real(phi(i*D,0) * conj(phi(i*D + 1,0)));
    }
    return sx_sum / (double)N;
}

double SpiralVM::omega_ang_end(int n) const {
    double t_end = n * T;
    double quasi = is_ang ? (sin(omega * M_PI * t_end / T) + sin(2*omega * M_PI * t_end / T)) : 0.0;
    return omega_ang_base + quasi;
}

double SpiralVM::h_effective_end(int n) {
    double t_end = n * T;
    double h1_osc = h1 * cos(omega * (t_end / 2.0) + M_PI/4.0);
    double h_eff = h0 + h1_osc;
    double h1_limit = 10000.0;
    if (h_eff / J > h1_limit) h_eff = h1_limit * J;
    else if (h_eff / J < -h1_limit) h_eff = -h1_limit * J;
    return h_eff;
}

int SpiralVM::get_period() {
    return current_period + 1;
}

double SpiralVM::compute_avg_stabilizer(const arma::cx_mat& phi_inp) {
    double stab_sum = 0.0;
    int count = 0;
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            int i = row * cols + col;
            int j_right = row * cols + ((col + 1) % cols);
            int j_down = ((row + 1) % rows) * cols + col;
            double zi = abs(phi_inp(i*D,0) - phi_inp(i*D + 1,0));
            double zj_right = abs(phi_inp(j_right*D,0) - phi_inp(j_right*D + 1,0));
            double zj_down = abs(phi_inp(j_down*D,0) - phi_inp(j_down*D + 1,0));
            stab_sum += zi * zj_right + zi * zj_down;
            count += 2;
        }
    }
    return (count>0) ? stab_sum / (double)count : 0.0;
}

// Dumps effective Hamiltonian components to CSV for the current period.
// Includes: period, site_idx, local_hx (transverse field per site), omega_ang,
// and global tone details (amp,freq,phase per tone).
// Call this in step_period after computing local_hx and omega_ang.
void SpiralVM::dump_h_eff(const std::string& fname_base, int period) const {
    std::string fname = fname_base + "_h_eff_period" + (period < 0 ? "_latest" : std::to_string(period)) + ".csv";
    std::ofstream fout(fname);
    if (!fout) {
        std::cerr << "[SpiralVM] Failed to open " << fname << " for h_eff dump\n";
        return;
    }

    int curr_period = (period < 0) ? current_period : period;
    // Header
    fout << "# SpiralVM h_eff Dump\n";
    fout << "# Period: " << curr_period << "\n";
    fout << "# T (s): " << T << "\n";
    fout << "# omega_ang: " << omega_ang_end(curr_period) << "\n";  // End-of-period angular mod
    fout << "# J: " << J << "\n";
    fout << "# Columns: site_idx,local_hx,drive_wid\n";

    // Per-site local_hx and drive index (computed at t_end for simplicity)
    double t_end = curr_period * T;
    for (int i = 0; i < N; ++i) {
        int wid = drive_index[i];
        double local_hx_val = h0;
        if (wid >= 0 && wid < (int)waveforms.size()) {
            local_hx_val += eval_waveform_with_envelope(waveforms[wid], t_end, drive_phase);
        } else {
            // Fallback (from your step_period)
            local_hx_val += h1 * cos(omega * (t_end / 2.0) + M_PI/4.0 + drive_phase);
        }
        fout << i << "," << local_hx_val << "," << wid << "\n";
    }

    // Section for global/merged tones (if compiled)
    fout << "# Tones (global if compiled): tone_idx,amp,freq,phase,logical_id\n";
    if (!waveforms.empty()) {
        const Waveform& global_w = waveforms[0];  // Assuming compiled to 0
        for (size_t t = 0; t < global_w.tones.size(); ++t) {
            const auto& tn = global_w.tones[t];
            fout << t << "," << tn.amp << "," << tn.freq << "," << tn.phase << "," << tn.logical_id << "\n";
        }
    }

    fout.close();
    std::cout << "[SpiralVM] Dumped h_eff to " << fname << "\n";
}

// Reconstructs approximate amplitude and phase for a logical qubit from a dumped waveform CSV.
// Assumes CSV from dump_waveforms("csv", ...) with columns: time_s,I,Q,amp,phase_rad
// Uses FFT on the complex IQ signal to find the component at the qubit's carrier freq.
// Returns a pair: <amplitude, phase> at the carrier.
// For in-memory: overload with bool from_file = true; if false, sample current global waveform in-memory.
std::pair<double, double> SpiralVM::reconstruct_logical_amp_phase_from_csv(const std::string& fname, uint32_t qid, bool from_file = true) {
    if (qid >= logical_qubits.size() || qid >= allocated_carriers.size()) {
        std::cout << "[SpiralVM] Invalid qid " << qid << " for reconstruction\n";
        return {0.0, 0.0};
    }
    double carrier_freq = allocated_carriers[qid];

    arma::vec times;
    arma::cx_vec iq;

    if (from_file) {
        // Load from CSV
        std::ifstream fin(fname);
        if (!fin) {
            std::cerr << "[SpiralVM] Failed to open " << fname << " for reconstruction\n";
            return {0.0, 0.0};
        }

        std::string line;
        std::vector<double> t_vec, i_vec, q_vec;
        while (std::getline(fin, line)) {
            if (line.empty() || line[0] == '#') continue;  // Skip headers/comments
            std::stringstream ss(line);
            std::string token;
            std::vector<std::string> cols;
            while (std::getline(ss, token, ',')) cols.push_back(token);
            if (cols.size() < 5) continue;  // time,I,Q,amp,phase
            t_vec.push_back(std::stod(cols[0]));
            i_vec.push_back(std::stod(cols[1]));
            q_vec.push_back(std::stod(cols[2]));
        }
        fin.close();

        times = arma::vec(t_vec);
        iq = arma::cx_vec(arma::vec(i_vec), arma::vec(q_vec));
    } else {
        // In-memory: sample current global waveform (id 0)
        double dt = T / 1000.0;  // Match dump resolution
        arma::vec amps, phases;  // Unused
        sample_waveform(waveforms[0], 0.0, dt, times, iq, amps, phases);
    }

    if (times.empty() || iq.empty()) {
        std::cout << "[SpiralVM] No data for reconstruction of q" << qid << "\n";
        return {0.0, 0.0};
    }

    // FFT of IQ signal
    arma::cx_vec fft_iq = arma::fft(iq);
    arma::vec freqs = arma::linspace<arma::vec>(0, 1.0 / (times(1) - times(0)), iq.n_elem);

    // Find index closest to carrier_freq
    arma::uword idx = arma::index_min(arma::abs(freqs - carrier_freq));
    std::complex<double> component = fft_iq(idx) / static_cast<double>(iq.n_elem);  // Normalize

    double amp = std::abs(component);
    double phase = std::arg(component);

    std::cout << "[SpiralVM] Reconstructed for q" << qid << ": amp=" << amp << ", phase=" << phase << "\n";
    return {amp, phase};
}