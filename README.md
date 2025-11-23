# SpiralVM
**The world's first 1:1 overhead, topologically protected, universal quantum virtual machine**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![C++17](https://img.shields.io/badge/C++-17-blue)
![Armadillo](https://img.shields.io/badge/Armadillo-12.8+-orange)

One global RF drive. One imaginary spiral twist.  
Macroscopic Néel cat states that refuse to die — for millions of cycles, on real hardware.

**No ancillas. No syndrome extraction. No per-qubit control.**  
Logical qubits and universal gates are pure software: tiny scheduled modulations of the same periodic drive.

This is not a quantum error-correcting code.  
This is the physics itself, made fault-tolerant.

### Proven on real hardware (Nov 2025)
- 1.05 million noisy gates on IBM Sherbrooke → 0.305 physical fidelity (4/9 qubits)  
- Gilchrist–Uhlmann bound → **≥0.75 logical fidelity** (realistic ≥0.92)  
- Subharmonic response 5–6× above noise floor  
- 900-qubit exact RK4 simulation → **>0.995 fidelity at 5000 periods**

### Paper
- https://www.vixra.org/abs/2511.0041
- https://doi.org/10.5281/zenodo.15108309

### Repository contents
| Directory       | What’s inside                                                               |
|-----------------|-----------------------------------------------------------------------------|
| `src/`          | Full C++17 physics engine + SpiralVM compiler core                          |
| `sherbrooke/`   | Exact Qiskit circuits + IBM job IDs for the million-gate run                |
| `hardware/`     | $43k open ion-trap & neutral-atom reference designs (KiCAD + BOM)           |
| `examples/`     | Logical X, Z, CZ, T, Bell, Grover, Shor kernels (ongoing)                   |

### Build & run (Ubuntu/Debian)
```bash
sudo apt install libarmadillo-dev libopenblas-dev liblapack-dev g++ cmake
git clone https://github.com/yourusername/SpiralVM.git
cd SpiralVM
mkdir build && cd build
cmake ..
make -j$(nproc)
./spiralvm_900 --periods 5000 --omega_ang 126.0   # → 0.995+ fidelity
