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
    // Assign a carrier frequency avoiding previously allocated carriers.
    // Very simple linear spacing allocator with basic collision avoidance.
    double base = freq_base;
    double spacing = freq_spacing;

    // Attempt to allocate at base + qid*spacing, but if too close to existing carrier, shift.
    double target = base + static_cast<double>(qid) * spacing;
    bool ok = false;
    int attempt = 0;
    const int max_attempts = 1000;
    while (!ok && attempt < max_attempts) {
        ok = true;
        for (double c : allocated_carriers) {
            if (fabs(c - target) < 0.5 * spacing) { ok = false; break; }
        }
        if (!ok) target += spacing * 0.5 * (1.0 + (attempt%2?1:-1)); // jitter
        attempt++;
    }
    if (attempt >= max_attempts) {
        // fallback: place after last
        target = base + allocated_carriers.size() * spacing * 1.2;
    }
    allocated_carriers.push_back(target);

    // Create waveform
    Waveform w = make_default_logical_waveform(qid);
    // set primary tone freq to allocated carrier
    if (!w.tones.empty()) w.tones[0].freq = target;

    // anti-alias filter
    lowpass_filter_waveform(w, lowpass_cutoff);

    int wid = waveforms.size();
    waveforms.push_back(w);
    return wid;
}

Waveform SpiralVM::make_default_logical_waveform(uint32_t qid) {
    // default: carrier + two weak sidebands + small static offset
    Waveform w;
    double carrier = freq_base + qid * freq_spacing;
    // main carrier (amplitude scaled to h1)
    Tone main(clamp_tone_amp(h1), carrier, 0.0, 0.0, 1.0);
    w.tones.push_back(main);
    // small symmetric sidebands (addressability / beat-note)
    double side_delta = 0.05 * freq_spacing + 0.01 * (1.0 + (qid%3));
    w.tones.push_back(Tone(clamp_tone_amp(0.12*h1), carrier + side_delta, M_PI/4.0));
    w.tones.push_back(Tone(clamp_tone_amp(0.12*h1), carrier - side_delta, -M_PI/4.0));
    // optional low-frequency envelope (slow amplitude modulation)
    w.tones.push_back(Tone(clamp_tone_amp(0.02*h1), 0.5 * freq_base, 0.0));
    return w;
}

// create a CZ-style waveform that briefly introduces correlated phase near both qubits
Waveform SpiralVM::make_cz_waveform(uint32_t qid1, uint32_t qid2, double strength, double duration_fraction) {
    // create combined waveform: carriers for each qubit plus a short burst tone at beat frequency
    Waveform w;
    double c1 = freq_base + qid1 * freq_spacing;
    double c2 = freq_base + qid2 * freq_spacing;
    // amplitude scaling
    double amp = clamp_tone_amp( strength * 0.7 );
    w.tones.push_back(Tone(amp, c1, 0.0));
    w.tones.push_back(Tone(amp, c2, 0.0));
    // add a short high-frequency burst (emulated by a tone with envelope timing handled at eval)
    double beat = fabs(c1 - c2);
    if (beat < 1e-12) beat = freq_spacing * 0.5;
    w.tones.push_back(Tone(clamp_tone_amp(0.5*strength), beat + max(c1,c2), 0.0));
    // optionally add a slow gaussian-like envelope via additional low-frequency tone
    w.tones.push_back(Tone(clamp_tone_amp(0.05*strength), 0.5*freq_base, 0.0));
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
double SpiralVM::eval_waveform_with_envelope(const Waveform &w, double t, double local_period_fraction) const {
    // local_period_fraction isn't used in this simple engine, but kept for future gating.
    double s = 0.0;
    for (const auto &tn : w.tones) {
        s += tn.amp * cos(tn.freq * t + tn.phase);
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

// ---------- Gate application (simple immediate implementation) ----------
void SpiralVM::apply_gate(const Gate& g, double period_time) {
    if (g.type == Gate::X) {
        apply_global_pi_pulse_on_even_cycles();
    } else if (g.type == Gate::Z) {
        apply_phase_shift(g.angle);
    } else if (g.type == Gate::CZ) {
        uint32_t t = g.target;
        uint32_t c = g.control;
        if (t < logical_qubits.size() && c < logical_qubits.size()) {
            // compile CZ waveform and insert it briefly
            Waveform w = make_cz_waveform(t,c, g.angle > 0.0 ? g.angle : 0.12, 0.02);
            int wid = waveforms.size();
            waveforms.push_back(w);
            // attach to neighborhoods (temporary injection)
            for (int row = (int)logical_qubits[t].center_y - R; row <= (int)logical_qubits[t].center_y + R; ++row) {
                if (row < 0 || row >= rows) continue;
                for (int col = (int)logical_qubits[t].center_x - R; col <= (int)logical_qubits[t].center_x + R; ++col) {
                    if (col < 0 || col >= cols) continue;
                    int phys_idx = row*cols + col;
                    drive_index[phys_idx] = wid;
                }
            }
            for (int row = (int)logical_qubits[c].center_y - R; row <= (int)logical_qubits[c].center_y + R; ++row) {
                if (row < 0 || row >= rows) continue;
                for (int col = (int)logical_qubits[c].center_x - R; col <= (int)logical_qubits[c].center_x + R; ++col) {
                    if (col < 0 || col >= cols) continue;
                    int phys_idx = row*cols + col;
                    drive_index[phys_idx] = wid;
                }
            }
            // run a few periods with this waveform active
            run_periods(max<uint32_t>(1, (uint32_t)round(period_time / T)));
            // restore waveforms for neighborhoods to their original allocated wf id
            for (int row = (int)logical_qubits[t].center_y - R; row <= (int)logical_qubits[t].center_y + R; ++row) {
                if (row < 0 || row >= rows) continue;
                for (int col = (int)logical_qubits[t].center_x - R; col <= (int)logical_qubits[t].center_x + R; ++col) {
                    if (col < 0 || col >= cols) continue;
                    int phys_idx = row*cols + col;
                    drive_index[phys_idx] = logical_qubits[t].waveform_id;
                }
            }
            for (int row = (int)logical_qubits[c].center_y - R; row <= (int)logical_qubits[c].center_y + R; ++row) {
                if (row < 0 || row >= rows) continue;
                for (int col = (int)logical_qubits[c].center_x - R; col <= (int)logical_qubits[c].center_x + R; ++col) {
                    if (col < 0 || col >= cols) continue;
                    int phys_idx = row*cols + col;
                    drive_index[phys_idx] = logical_qubits[c].waveform_id;
                }
            }
            // pop waveform (or keep if you want persistent)
            // waveforms.pop_back();
        }
    } else if (g.type == Gate::MEASURE) {
        cout << "[SpiralVM] Measure gate requested for " << g.target << "\n";
    } else if (g.type == Gate::PHASE) {
        apply_phase_shift(g.angle);
    } else {
        cout << "[SpiralVM] Gate type not implemented in apply_gate()\n";
    }
}

// compile and run simple program - immediate mode
void SpiralVM::compile_and_run(const std::vector<Gate>& program) {
    for (const auto &g : program) {
        apply_gate(g, T); // default run for ~1 period
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
    double energy_init = real(inner_product_cl10(phi, Hphi_init)) + compute_zz_energy(phi, J, 0, 0);
    double zz_energy_init = compute_zz_energy(phi, J, 0, 0);
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
        double zz_energy = is_ang ? compute_zz_energy(phi, J, omega_ang_end(n), n + 1, true) : compute_zz_energy(phi, J, 0, 0);
        double energy = sx_energy + zz_energy;

        fout << n + 1 << "," << fidelities[n + 1] << "," << compute_avg_stabilizer(phi) << "," << energy << "," << sx_energy << "," << zz_energy << "," << delta_F << "," << h_effective_end(n) << "," << sx_avg(n) << "\n";
        cout << "Period " << n + 1 << ": Fidelity = " << fidelities[n + 1] << ", Energy = " << energy << "\n";
    }
    fout.close();
}

// ---------- RK4 integration step (modified to evaluate per-site waveform) ----------
void SpiralVM::step_period(int n, double &delta_F) {
    double dt = T / steps;
    double delta_F_target = 1.0;

    cx_mat phi_new = phi;

    // precompute local hx per site for this period/time sample
    std::vector<double> local_hx(N, 0.0);

    for (int k = 0; k < steps; k++) {
        double t = n * T + k * dt;

        // fill local_hx by evaluating the waveform bank for each physical site
        for (int i = 0; i < N; ++i) {
            int wid = drive_index[i];
            if (wid < 0 || wid >= (int)waveforms.size()) {
                local_hx[i] = h0 + h1 * cos(omega * (t/2.0) + M_PI/4.0);
            } else {
                // local eval
                double val = eval_waveform_with_envelope(waveforms[wid], t, 0.0);
                local_hx[i] = h0 + val;
            }
        }

        // compute angular modulation as before
        double angular_freq_quasi = is_ang ? (sin(omega * M_PI * t / T) + sin(2*omega * M_PI * t / T)) : 0.0;
        double omega_ang = omega_ang_base + angular_freq_quasi;

        // Build inhomogeneous Hamiltonian for this time step
        sp_cx_mat H_sx = hamiltonian_cl10_90_spiral_twist_inhomogeneous(J, local_hx, omega_ang);

        // compute Hzz contribution vector
        cx_mat Hzz_phi = compute_zz_energy_vector(phi_new, J, omega_ang, n, is_ang);

        // RK4 k1..k4 using H_sx and Hzz_phi
        cx_mat k1 = mat_vec_mult_cl10(H_sx, phi_new) + Hzz_phi;
        cx_mat phi_temp = phi_new + (-cx_double(0,1) * dt/2.0) * k1;

        // recompute local_hx at mid-step for better accuracy (optional here)
        // We keep local_hx constant within step for performance; that's already a decent approximation.

        cx_mat k2 = mat_vec_mult_cl10(H_sx, phi_temp) + compute_zz_energy_vector(phi_temp, J, omega_ang, n, is_ang);
        phi_temp = phi_new + (-cx_double(0,1) * dt/2.0) * k2;
        cx_mat k3 = mat_vec_mult_cl10(H_sx, phi_temp) + compute_zz_energy_vector(phi_temp, J, omega_ang, n, is_ang);
        phi_temp = phi_new + (-cx_double(0,1) * dt) * k3;
        cx_mat k4 = mat_vec_mult_cl10(H_sx, phi_temp) + compute_zz_energy_vector(phi_temp, J, omega_ang, n, is_ang);

        phi_new += (-cx_double(0,1) * dt/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);

        double normrk4 = sqrt(real(inner_product_cl10(phi_new, phi_new)));
        if (normrk4 > 0) phi_new /= normrk4;
    }

    phi = phi_new;
    state = phi;
    // fidelity tracking
    double fid = fabs(inner_product_cl10(phi_in, phi));
    // push into vector (ensure size)
    if ((int)fidelities.size() <= current_period+1) fidelities.resize(current_period+2, 0.0);
    fidelities[current_period+1] = fid;
    delta_F = 1.0 - fid;
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

double SpiralVM::measure_logical_Z(uint32_t qid) const {
    if (qid >= logical_qubits.size()) return 0.0;
    // measure global staggered magnetization as before
    double staggered = 0.0;
    int count = 0;
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            int i = r*cols + c;
            double sz = std::norm(phi(i*D + 0,0)) - std::norm(phi(i*D + 1,0));
            staggered += ((r+c)%2==0) ? sz : -sz;
            ++count;
        }
    }
    return staggered / (double)count;
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
    const LogicalQubit &q1 = logical_qubits[qid1];
    const LogicalQubit &q2 = logical_qubits[qid2];
    double sum_correlation = 0.0;
    int count = 0;
    for (int row1 = q1.center_y - R; row1 <= q1.center_y + R; ++row1) {
        if (row1 < 0 || row1 >= rows) continue;
        for (int col1 = q1.center_x - R; col1 <= q1.center_x + R; ++col1) {
            if (col1 < 0 || col1 >= cols) continue;
            int i1 = row1*cols + col1;
            double z1 = std::norm(phi(i1*D + 0,0)) - std::norm(phi(i1*D + 1,0));
            for (int row2 = q2.center_y - R; row2 <= q2.center_y + R; ++row2) {
                if (row2 < 0 || row2 >= rows) continue;
                for (int col2 = q2.center_x - R; col2 <= q2.center_x + R; ++col2) {
                    if (col2 < 0 || col2 >= cols) continue;
                    int i2 = row2*cols + col2;
                    double z2 = std::norm(phi(i2*D + 0,0)) - std::norm(phi(i2*D + 1,0));
                    sum_correlation += z1 * z2;
                    ++count;
                }
            }
        }
    }
    if (count==0) return 0.0;
    return sum_correlation / (double)count;
}

void SpiralVM::apply_phase_shift(double angle) {
    cx_double phase = std::exp(cx_double(0, angle));
    for (int i = 0; i < N*D; ++i) phi(i,0) *= phase;
}

void SpiralVM::apply_phase_kick_between(uint32_t qid1, uint32_t qid2, double strength, double duration_fraction) {
    // lightweight approximate implementation: briefly multiply local |1> amplitudes by phase on neighborhoods
    if (qid1 >= logical_qubits.size() || qid2 >= logical_qubits.size()) return;
    const LogicalQubit &q1 = logical_qubits[qid1];
    const LogicalQubit &q2 = logical_qubits[qid2];
    cx_double phase = std::exp(cx_double(0, strength * duration_fraction));
    for (int row1 = q1.center_y - R; row1 <= q1.center_y + R; ++row1) {
        if (row1 < 0 || row1 >= rows) continue;
        for (int col1 = q1.center_x - R; col1 <= q1.center_x + R; ++col1) {
            if (col1 < 0 || col1 >= cols) continue;
            int i1 = row1*cols + col1;
            phi(i1*D + 1,0) *= phase;
        }
    }
    for (int row2 = q2.center_y - R; row2 <= q2.center_y + R; ++row2) {
        if (row2 < 0 || row2 >= rows) continue;
        for (int col2 = q2.center_x - R; col2 <= q2.center_x + R; ++col2) {
            if (col2 < 0 || col2 >= cols) continue;
            int i2 = row2*cols + col2;
            phi(i2*D + 1,0) *= phase;
        }
    }
    // renormalize
    double norm = sqrt(real(inner_product_cl10(phi, phi)));
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


/*
int main() {
    SpiralVM vm(30, 30);
    vm.initialize_state("neel");
    vm.run_floquet(5000, "neel");
    return 0;
}*/
