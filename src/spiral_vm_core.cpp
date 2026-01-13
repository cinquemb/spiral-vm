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

    Waveform physical_global = waveforms[0];

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

    for (const auto& tn : physical_global.tones) {
        tone_map[{tn.freq, tn.phase, tn.logical_id}] += tn.amp;
    }

    for (size_t wid = 1; wid < waveforms.size(); ++wid)
        for (const auto& tn : waveforms[wid].tones)
            tone_map[{tn.freq, tn.phase, tn.logical_id}] += tn.amp;

    physical_global.tones.clear();
    for (const auto& entry : tone_map) {
        double freq  = std::get<0>(entry.first);
        double phase = std::get<1>(entry.first);
        int    qid   = std::get<2>(entry.first);
        double amp   = clamp_tone_amp(entry.second);
        physical_global.tones.emplace_back(amp, freq, phase, 0.0, 1.0, qid);
    }

    waveforms.clear();
    waveforms.push_back(physical_global);
    std::fill(drive_index.begin(), drive_index.end(), 0);

    cout << "[SpiralVM] Compiled to single physical global waveform with "
         << physical_global.tones.size()
         << " merged tones (logical IDs preserved)\n";
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
Waveform SpiralVM::make_cz_waveform(uint32_t qid1, uint32_t qid2, double strength, double duration_fraction) {
    // create combined waveform: carriers for each qubit plus a short burst tone at beat frequency
    Waveform w;
    double c1 = freq_base + qid1 * freq_spacing;
    double c2 = freq_base + qid2 * freq_spacing;
    double amp = clamp_tone_amp(strength * 0.7);
    w.tones.push_back(Tone(amp, c1, 0.0, 0.0, 1.0, qid1));
    w.tones.push_back(Tone(amp, c2, 0.0, 0.0, 1.0, qid2));
    double beat = fabs(c1 - c2);
    if (beat < 1e-12) beat = freq_spacing * 0.5;
    w.tones.push_back(Tone(clamp_tone_amp(0.5*strength), beat + max(c1,c2), 0.0, 0.0, 1.0, -1));  // shared
    w.tones.push_back(Tone(clamp_tone_amp(0.05*strength), 0.5*freq_base, 0.0, 0.0, 1.0, -1));  // shared
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
    // local_period_fraction isn't used in this simple engine, but kept for future gating.
    double s = 0.0;
    for (const auto &tn : w.tones) {
        s += tn.amp * cos(tn.freq * t + tn.phase + local_phase);
    }
    return s;
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
        fout << "  \"waveforms\": [\n";

        //for (size_t wid = 0; wid < waveforms.size(); ++wid) {
            if (wid > 0) fout << ",\n";
            fout << "    {\n";
            fout << "      \"id\": " << wid << ",\n";
            fout << "      \"tones\": [\n";
            for (size_t t = 0; t < waveforms[wid].tones.size(); ++t) {
                if (t > 0) fout << ",\n";
                const auto& tn = waveforms[wid].tones[t];
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
            sample_waveform(waveforms[wid], 0.0, dt, times, iq, amps, phases);
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

// ---------- Gate application (simple immediate implementation) ----------
void SpiralVM::apply_gate(const Gate& g, double period_time) {
    switch (g.type) {  // ← ADD THIS
        case Gate::X:
            global_pi_pulse();
            break;
        case Gate::Z:
            apply_phase_shift(g.angle);
            break;
        case Gate::CZ: {
            global_pi_pulse();
            uint32_t t = g.target;
            uint32_t c = g.control;
            apply_phase_kick_between(t, c, 0.25, period_time*T);  // strength & duration
            break;
        }
        case Gate::MEASURE:
            cout << "[SpiralVM] Measure gate requested for " << g.target << "\n";
            break;
        case Gate::PHASE:
            apply_phase_shift(g.angle);
            break;
        case Gate::H:
            logical_hadamard(g.target);
            break;
        case Gate::CNOT: {
            // Standard: H target, CZ, H target
            logical_hadamard(g.target);
            apply_phase_kick_between(g.control, g.target, M_PI, 0.2 * period_time);  // full π for CZ
            logical_hadamard(g.target);
            break;
        }
        case Gate::RX: case Gate::RY: case Gate::RZ:  // arbitrary rotation
            // Placeholder: implement via phase ramp + pi-pulse sequences later
            cout << "[SpiralVM] Rotation gate " << g.type << " not fully implemented\n";
            break;
        default:
            cout << "[SpiralVM] Gate type not implemented in apply_gate()\n";
            break;
    }
}

// compile and run simple program - immediate mode
void SpiralVM::compile_and_run(const std::vector<Gate>& program) {
    for (const auto &g : program) {
        apply_gate(g, 1); // default run for ~1 period
    }
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

            if (wid < 0 || wid >= (int)waveforms.size()) {
                // simple oscillation using scattered frequency
                local_hx[i] = h0 + h1 * cos(qubit_freqs[i] * t + M_PI/4.0 + drive_phase);
            } else {
                // evaluate waveform if assigned
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
    delta_F = 1.0 - fid;

    std::string fname_base = "waveform_" + std::to_string(current_period);
    dump_waveforms("csv", fname_base, -1);

    current_period++;
}


// ---------- run N periods ----------------
void SpiralVM::run_periods(uint32_t N_periods) {
    double dummy_deltaF = 0.0;
    for (uint32_t i = 0; i < N_periods; i++) {
        step_period(current_period, dummy_deltaF);
    }
}

// ---------- logical measurements and gates ----------
void SpiralVM::apply_global_pi_pulse_on_even_cycles() {
    if (current_period % 2 == 0) {
        cout << "[SpiralVM] Applying logical X (global pi pulse) at period " << current_period << "\n";
        global_pi_pulse();
    } else {
        cout << "[SpiralVM] Skipping logical X on odd period " << current_period << "\n";
    }
}

void SpiralVM::global_pi_pulse() {
    // swap |0> and |1> amplitudes for every site
    for (int i = 0; i < N; ++i) {
        std::swap(phi(i*D + 0, 0), phi(i*D + 1, 0));
    }
    state = phi;
    std::string fname_base = "logical_x_waveform_period_"+ std::to_string(current_period);
        dump_waveforms("csv", fname_base, current_period);
}

// logical Hadamard on single qubit
void SpiralVM::logical_hadamard(uint32_t qid) {
    if (qid >= logical_qubits.size()) return;

    // Z(π/2) only on target neighborhood
    logical_phase_ramp(qid, M_PI/2.0, 1);  // ramp over 1 periods

    // Global X
    global_pi_pulse();

    // Z(-π/2) on target
    logical_phase_ramp(qid, -M_PI/2.0, 1);
}

double SpiralVM::measure_even_population(uint32_t qid) {
    if (qid >= logical_qubits.size()) return 0.0;
    const LogicalQubit &q = logical_qubits[qid];
    double pop_sum = 0.0;
    int count = 0;
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            int dist = abs(row - (int)q.center_y) + abs(col - (int)q.center_x);
            if (dist <= R && ((row + col) % 2 == 0)) {
                int i = row*cols + col;
                pop_sum += std::norm(phi(i*D,0));
                ++count;
            }
        }
    }
    if (count==0) return 0.0;
    return pop_sum / count;
}

double SpiralVM::measure_logical_global_Z(uint32_t qid) const {
    if (qid >= logical_qubits.size()) return 0.0;
    double staggered = 0.0;
    int count = 0;
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            int i = r*cols + c;
            double sz = N * (std::norm(phi(i*D + 0,0)) - std::norm(phi(i*D + 1,0)));  // Add N factor
            staggered += ((r+c)%2==0) ? sz : -sz;
            ++count;
        }
    }
    return staggered / (double)count;
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

void SpiralVM::logical_phase_ramp(int target_qid, double slope, int steps) {
    for (int k = 0; k < steps; ++k) {
        for (auto &w : waveforms) {
            for (auto &tn : w.tones) {
                if (tn.logical_id == target_qid) {
                    tn.phase += slope;
                }
            }
        }
        run_periods(1);
    }
}


double SpiralVM::get_logical_phase(uint32_t qid) {
    if (qid >= logical_qubits.size()) return 0.0;
    const LogicalQubit &q = logical_qubits[qid];
    complex<double> phase_sum(0.0, 0.0);
    int count = 0;
    for (int row = q.center_y - R; row <= q.center_y + R; ++row) {
        if (row < 0 || row >= rows) continue;
        for (int col = q.center_x - R; col <= q.center_x + R; ++col) {
            if (col < 0 || col >= cols) continue;
            int i = row*cols + col;
            complex<double> amp1 = phi(i*D + 1,0);
            if (abs(amp1) > 1e-12) {
                phase_sum += amp1 / abs(amp1);
                ++count;
            }
        }
    }
    if (count==0) return 0.0;
    complex<double> avg_phase = phase_sum / (double)count;
    return arg(avg_phase);
}


double SpiralVM::logical_zz_correlation(uint32_t qid1, uint32_t qid2) {
    if (qid1 >= logical_qubits.size() || qid2 >= logical_qubits.size()) return 0.0;
    double z1 = measure_logical_Z(qid1);
    double z2 = measure_logical_Z(qid2);
    return z1 * z2;
}

void SpiralVM::apply_phase_shift(double angle) {
    drive_phase += angle;
}

void SpiralVM::global_phase_ramp(double slope, int steps) {
    for (int k = 0; k < steps; ++k) {
        apply_phase_shift(slope);  // small increment each step
        run_periods(1);            // let dynamics respond
    }
}


void SpiralVM::apply_phase_kick_between(uint32_t qid1, uint32_t qid2,
                                        double strength, double duration_fraction) {
    if (qid1 >= logical_qubits.size() || qid2 >= logical_qubits.size()) return;

    const LogicalQubit &q1 = logical_qubits[qid1];
    const LogicalQubit &q2 = logical_qubits[qid2];

    // Estimate logical Z for each qubit (coarse, local)
    auto logical_Z_estimate = [&](const LogicalQubit &q) {
        double sum = 0.0;
        int count = 0;
        for (int row = q.center_y - R; row <= q.center_y + R; ++row) {
            if (row < 0 || row >= rows) continue;
            for (int col = q.center_x - R; col <= q.center_x + R; ++col) {
                if (col < 0 || col >= cols) continue;
                int i = row * cols + col;
                double sz = N * (std::norm(phi(i*D + 0,0)) - std::norm(phi(i*D + 1,0)));  // Add N factor
                sum += ((row + col) % 2 == 0) ? sz : -sz;
                ++count;
            }
        }
        return (count > 0) ? sum / (double)count : 0.0;
    };

    double Z1 = logical_Z_estimate(q1);
    double Z2 = logical_Z_estimate(q2);

    // Conditional phase angle ~ Z1 * Z2
    double theta = strength * duration_fraction * Z1 * Z2;
    cx_double phase = std::exp(cx_double(0, theta));

    // Apply only to |1> amplitudes in BOTH neighborhoods
    for (int row = q1.center_y - R; row <= q1.center_y + R; ++row) {
        if (row < 0 || row >= rows) continue;
        for (int col = q1.center_x - R; col <= q1.center_x + R; ++col) {
            if (col < 0 || col >= cols) continue;
            int i = row * cols + col;
            phi(i*D + 1,0) *= phase;
        }
    }

    for (int row = q2.center_y - R; row <= q2.center_y + R; ++row) {
        if (row < 0 || row >= rows) continue;
        for (int col = q2.center_x - R; col <= q2.center_x + R; ++col) {
            if (col < 0 || col >= cols) continue;
            int i = row * cols + col;
            phi(i*D + 1,0) *= phase;
        }
    }

    // Renormalize
    double norm = std::sqrt(real(inner_product_cl10(phi, phi)));
    if (norm > 0) phi /= norm;
}


void SpiralVM::apply_phase_kick_between_full(uint32_t qid1, uint32_t qid2, double strength, double duration_fraction) {
    // more detailed: temporarily inject a CZ waveform (uses make_cz_waveform) and run for one period
    if (qid1 >= logical_qubits.size() || qid2 >= logical_qubits.size()) return;
    Waveform w = make_cz_waveform(qid1,qid2,strength,duration_fraction);
    int wid = waveforms.size();
    waveforms.push_back(w);
    // attach to neighborhoods
    for (int row = logical_qubits[qid1].center_y - R; row <= logical_qubits[qid1].center_y + R; ++row) {
        if (row < 0 || row >= rows) continue;
        for (int col = logical_qubits[qid1].center_x - R; col <= logical_qubits[qid1].center_x + R; ++col) {
            if (col < 0 || col >= cols) continue;
            drive_index[row*cols + col] = wid;
        }
    }
    for (int row = logical_qubits[qid2].center_y - R; row <= logical_qubits[qid2].center_y + R; ++row) {
        if (row < 0 || row >= rows) continue;
        for (int col = logical_qubits[qid2].center_x - R; col <= logical_qubits[qid2].center_x + R; ++col) {
            if (col < 0 || col >= cols) continue;
            drive_index[row*cols + col] = wid;
        }
    }
    // run one or a few periods
    run_periods(1);
    // restore allocations
    for (int row = logical_qubits[qid1].center_y - R; row <= logical_qubits[qid1].center_y + R; ++row) {
        if (row < 0 || row >= rows) continue;
        for (int col = logical_qubits[qid1].center_x - R; col <= logical_qubits[qid1].center_x + R; ++col) {
            if (col < 0 || col >= cols) continue;
            drive_index[row*cols + col] = logical_qubits[qid1].waveform_id;
        }
    }
    for (int row = logical_qubits[qid2].center_y - R; row <= logical_qubits[qid2].center_y + R; ++row) {
        if (row < 0 || row >= rows) continue;
        for (int col = logical_qubits[qid2].center_x - R; col <= logical_qubits[qid2].center_x + R; ++col) {
            if (col < 0 || col >= cols) continue;
            drive_index[row*cols + col] = logical_qubits[qid2].waveform_id;
        }
    }
    // optionally remove the waveform
    // waveforms.pop_back();
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

// ---------- compute zz energy and vector ----------
double SpiralVM::compute_zz_energy(const arma::cx_mat& phi_inp, double J, double omega_ang, double period, bool is_angf) {
    // Keep your original energy computation semantics, but simplified numerics
    double energy = 0.0;
    double theta_max_base = M_PI / 512.0;
    double theta_max = theta_max_base * (1.0 - cos(M_PI * period / 2.0)) / 2.0;
    double center_x = cols / 2.0, center_y = rows / 2.0;
    double r_max = sqrt(center_x*center_x + center_y*center_y);
    double norm_factor = sqrt((double)rows * cols);
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            int i = row * cols + col;
            int j_right = row * cols + ((col + 1) % cols);
            int j_down = ((row + 1) % rows) * cols + col;
            double r = sqrt(sq(row - center_y) + sq(col - center_x));
            double phi_ang = atan2(row - center_y, col - center_x);
            double theta = (omega_ang == 0.0) ? 0.0 : theta_max * (r / r_max) * cos(omega_ang * phi_ang);
            cx_double J_twist = is_angf ? J * cx_double(0, sin(theta)) : cx_double(J,0);
            // approximate z projection:
            cx_double zi = (phi_inp(i*D,0) - phi_inp(i*D + 1,0)) * norm_factor;
            cx_double zjr = (phi_inp(j_right*D,0) - phi_inp(j_right*D + 1,0)) * norm_factor;
            cx_double zjd = (phi_inp(j_down*D,0) - phi_inp(j_down*D + 1,0)) * norm_factor;
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

arma::cx_mat SpiralVM::compute_zz_energy_vector(const arma::cx_mat& phi_inp, double J, double omega_ang, double period, bool is_angf) {
    arma::cx_mat Hzz_phi = arma::zeros<arma::cx_mat>(N*D, 1);
    double theta_max_base = M_PI / 512.0;
    double theta_max = theta_max_base * (1.0 - cos(M_PI * period / 2.0)) / 2.0;
    double center_x = cols / 2.0, center_y = rows / 2.0;
    double r_max = sqrt(center_x*center_x + center_y*center_y);
    double norm_factor = sqrt((double)N);
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            int i = row * cols + col;
            int j_right = row * cols + ((col + 1) % cols);
            int j_down = ((row + 1) % rows) * cols + col;
            double r = sqrt(sq(row - center_y) + sq(col - center_x));
            double phi_ang = atan2(row - center_y, col - center_x);
            double theta = (omega_ang == 0.0) ? 0.0 : theta_max * (r / r_max) * cos(omega_ang * phi_ang);
            cx_double J_twist = is_angf ? J * cx_double(0, sin(theta)) : cx_double(J,0);
            cx_double zjr = (phi_inp(j_right*D,0) - phi_inp(j_right*D + 1,0)) * norm_factor;
            cx_double zjd = (phi_inp(j_down*D,0) - phi_inp(j_down*D + 1,0)) * norm_factor;
            Hzz_phi(i*D,0) = J_twist * (zjr + zjd) * phi_inp(i*D,0);
            Hzz_phi(i*D + 1,0) = J_twist * (-zjr - zjd) * phi_inp(i*D + 1,0);
        }
    }
    return Hzz_phi;
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

double SpiralVM::omega_ang_end(int n) {
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

void SpiralVM::print_overlap_stats() {
    if (!overlap_enabled) {
        cout << "[SpiralVM] Overlap mode disabled.\n";
        return;
    }
    for (int i = 0; i < N; ++i) {
        int count = phys_to_logicals[i].size();
        if (count > 1) {
            cout << "Physical qubit " << i << " belongs to " << count << " logical qubits.\n";
        }
    }
}