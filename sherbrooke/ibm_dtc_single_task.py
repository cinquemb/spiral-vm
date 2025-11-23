import numpy as np
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit, ClassicalRegister, transpile
from qiskit_ibm_runtime import QiskitRuntimeService, Sampler
from qiskit.transpiler import PassManager, CouplingMap
from qiskit_ibm_runtime.transpiler.passes import FoldRzzAngle
from qiskit.transpiler.target import Target
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

from qiskit.transpiler.passes import (
    SetLayout,
    ApplyLayout,
    SabreSwap,
    LookaheadSwap,
    CheckMap,
    Optimize1qGates,
)
from qiskit.transpiler.layout import Layout


from time import sleep
import sys

# Parameters (same as your Braket code)
ROWS, COLS = 3, 3
N = ROWS * COLS
J = 0.3
h1 = 0.568
hz = 20
omega = hz * 2 * np.pi
T = 2 * np.pi / omega
steps = 5
dt = T / steps
omega_ang = 2 * omega
periods = 50
shots = 1024

k_oA = 1
alpha_ang = 1
theta_max = np.pi / 512

# Define pairs for 3x3 lattice with periodic boundaries
pairs = []
for i in range(ROWS):
    for j in range(COLS):
        q = i * COLS + j
        pairs.append((q, (i * COLS + (j + 1) % COLS)))  # Horizontal
        pairs.append((q, (((i + 1) % ROWS) * COLS + j)))  # Vertical

def theta_ij(t, pair):
    r_max = np.sqrt(8)
    xi, yi = pair[0] // COLS, pair[0] % COLS
    xj, yj = pair[1] // COLS, pair[1] % COLS
    r = np.sqrt((xi - ROWS/2)**2 + (yi - COLS/2)**2)  # Center-adjusted
    phi = np.arctan2(yj - yi, xj - xi)
    n = int(t / T)
    factor = (1 - np.cos(np.pi * n / 2)) / 2
    omega_app = 0 + k_oA * np.sin(omega * np.pi * t / T) + alpha_ang * (np.sin(omega * np.pi * t / T) + np.sin(omega_ang * np.pi * t / T))
    return theta_max * factor * (r / r_max) * np.cos(omega_app * phi)

# Initialize Qiskit circuit with classical bits for mid-circuit measurements
circ = QuantumCircuit(N, N * periods)

# Prepare Néel state
for i in range(ROWS):
    for j in range(COLS):
        if (i + j) % 2 == 1:
            circ.x(i * COLS + j)

# Copy initial circuit for fidelity baseline (no evolution, just Néel state)
initial_circ = circ.copy()

# Main simulation loop: evolve and measure at each period
for n in range(periods):
    for step in range(steps):
        t = n * T + step * dt
        h_t = h1 * np.cos(omega * (t / 2) + np.pi / 4)
        theta_x = -2 * h_t * dt

        # Apply rzz gates (Qiskit's zz rotation)
        for (i, j) in pairs:
            sin_theta = np.sin(theta_ij(t, (i, j)))
            theta_zz = J * sin_theta * dt
            circ.rzz(-theta_zz, i, j)

        # Apply rx rotations
        for q in range(N):
            circ.rx(theta_x, q)

        # Repeat rzz gates after rx
        for (i, j) in pairs:
            sin_theta = np.sin(theta_ij(t, (i, j)))
            theta_zz = J * sin_theta * dt
            circ.rzz(-theta_zz, i, j)

    # Mid-circuit measurement for fidelity estimation
    for q in range(N):
        circ.measure(q, n * N + q)

# Setup IBM Quantum service (replace with your API key and CRN)
your_api_key = ""
your_crn = "crn:v1:bluemix:public:quantum-computing:us-east:a/e596892d7b30494c83fc2e71b6a3e1a0:690d156c-0928-4625-9c79-9fd82fa31b37::"

'''
QiskitRuntimeService.save_account(
    channel="ibm_quantum_platform",
    token=your_api_key,
    instance=your_crn,
    name="qgss-2025-fun-2",
)
'''

service = QiskitRuntimeService(name="qgss-2025-fun-2")
# Choose least busy real backend (not simulator)
backend = service.least_busy(operational=True, simulator=False, use_fractional_gates=True)
print(f"Running on backend: {backend.name}")



print(type(circ.qubits))
print(circ.qubits)
print([type(q) for q in circ.qubits])


# Create a custom pass manager with Optimize1qGates eps set very low
pass_manager = generate_preset_pass_manager(backend=backend, optimization_level=1)

# Routing (SWAP insertion)
#pass_manager.append(Optimize1qGates(eps=1e-15))
#pass_manager.append(FoldRzzAngle())





# Transpile circuit for backend with optimization
#circ_transpiled = transpile(circ, backend=backend, optimization_level=0)

circ_transpiled = pass_manager.run(circ)

initial_circ_transpiled = transpile(initial_circ, backend=backend, optimization_level=0)

circ_transpiled_counts = circ_transpiled.count_ops()
circ_transpiled_depth = circ_transpiled.depth()
print("circ_transpiled:", circ_transpiled_counts)
print("circ_transpiled_depth:", circ_transpiled_depth)

#sys.exit()


# Get all operational real backends with >= 9 qubits (your circuit size)
backends = service.backends(simulator=False, operational=True)
suitable_backends = [b for b in backends if b.configuration().n_qubits >= 9]

# Collect average T1 and T2 per backend
backend_stats = []
for backend1 in suitable_backends:
    props = backend1.properties()
    t1_vals = [q[0].value for q in props.qubits if q[0].name == 'T1']
    t2_vals = [q[1].value for q in props.qubits if q[1].name == 'T2']
    avg_t1 = sum(t1_vals) / len(t1_vals)
    avg_t2 = sum(t2_vals) / len(t2_vals)
    backend_stats.append((backend1.name, avg_t1, avg_t2))

# Sort by average T2 descending (better coherence)
backend_stats.sort(key=lambda x: x[2], reverse=True)

for name, avg_t1, avg_t2 in backend_stats:
    print(f"{name}: Avg T1 = {avg_t1:.1f} µs, Avg T2 = {avg_t2:.1f} µs")


# Run main evolved circuit
sampler = Sampler(backend)
job = sampler.run([circ_transpiled], shots=shots)

# Run initial Néel state circuit for baseline
initial_job = sampler.run([initial_circ_transpiled], shots=shots)

while True:
    job_status = job.status()
    initial_status = initial_job.status()

    print(f"Job status 1: {job_status}")
    print(f"Job status 2: {initial_status}")
    if job_status.is_final() and initial_status.is_final():
        break
    sleep(10)  # wait 10 seconds before checking again


job_result = job.result()
initial_result = initial_job.result()

# Extract counts (quasi-distributions)
counts = job_result[0].quasi_dists
initial_counts = initial_result[0].quasi_dists

# Compute fidelity for each period
fidelities = []
for n in range(periods):
    period_counts = {}
    for bitstring, prob in counts.items():
        # bitstring is a string of length N*periods, extract bits for period n
        period_outcome = bitstring[n * N:(n + 1) * N]
        period_counts[period_outcome] = period_counts.get(period_outcome, 0) + prob

    fidelity = sum(prob for outcome, prob in period_counts.items() if outcome in initial_counts) / shots
    fidelities.append(fidelity)
    print(f"Fidelity at {n + 1}T: {fidelity}")

# Save fidelities to file
with open("fidelity_results.txt", "w") as f:
    for n, fidelity in enumerate(fidelities):
        f.write(f"Fidelity at {n + 1}T: {fidelity}\n")

'''
'''

'''
$ python ibm_dtc_single_task.py 
Running on backend: ibm_sherbrooke
circ_transpiled: OrderedDict([('rz', 1984032), ('sx', 991920), ('ecr', 314808), ('x', 72004), ('measure', 1800)])
circ_transpiled_depth: 1050507
ibm_sherbrooke: Avg T1 = 274.2 µs, Avg T2 = 210.0 µs
ibm_brisbane: Avg T1 = 235.2 µs, Avg T2 = 137.6 µs
ibm_torino: Avg T1 = 170.6 µs, Avg T2 = 136.4 µs
'''

'''
import numpy as np
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit, ClassicalRegister, transpile
from qiskit_ibm_runtime import QiskitRuntimeService, Sampler
from qiskit.transpiler import PassManager
from qiskit.transpiler.passes import Optimize1qGates

from time import sleep
service = QiskitRuntimeService(name="qgss-2025-fun-2")

# Get all operational real backends with >= 9 qubits (your circuit size)
backends = service.backends(simulator=False, operational=True)
suitable_backends = [b for b in backends if b.configuration().n_qubits >= 9]

# Collect average T1 and T2 per backend
backend_stats = []
for backend1 in suitable_backends:
    props = backend1.properties()
    t1_vals = [q[0].value for q in props.qubits if q[0].name == 'T1']
    t2_vals = [q[1].value for q in props.qubits if q[1].name == 'T2']
    avg_t1 = sum(t1_vals) / len(t1_vals)
    avg_t2 = sum(t2_vals) / len(t2_vals)
    backend_stats.append((backend1.name, avg_t1, avg_t2))

# Sort by average T2 descending (better coherence)
backend_stats.sort(key=lambda x: x[2], reverse=True)


for name, avg_t1, avg_t2 in backend_stats:
    print(f"{name}: Avg T1 = {avg_t1:.1f} s, Avg T2 = {avg_t2:.1f} s")


ibm_sherbrooke: Avg T1 = 274.2 s, Avg T2 = 210.0 s
ibm_brisbane: Avg T1 = 235.2 s, Avg T2 = 137.6 s
ibm_torino: Avg T1 = 170.6 s, Avg T2 = 136.4 s
'''
