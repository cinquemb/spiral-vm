from qiskit_ibm_runtime import QiskitRuntimeService
import math

# Important qubit indices (classical bit indices) to focus on
IMPORTANT_BITS = [1, 3, 5, 7]

def databin_to_counts(databin):
    if not hasattr(databin, 'c'):
        raise ValueError("DataBin does not contain expected attribute 'c' for samples.")
    bitarray = databin.c
    counts = bitarray.get_counts()
    return counts

def normalize_counts(counts: dict, shots: int) -> dict:
    return {k: v / shots for k, v in counts.items()}

def extract_important_bits_per_period(bitstring, period_idx, period_size, important_bits):
    """
    Extract important bits for a given period from the bitstring.

    Args:
        bitstring (str): Full bitstring (MSB first).
        period_idx (int): Index of the period (0-based).
        period_size (int): Number of bits per period.
        important_bits (list[int]): Indices of important bits relative to period start.

    Returns:
        str: Extracted bits concatenated as a string.
    """
    # Calculate start and end indices of the period in the bitstring
    start = period_idx * period_size
    end = start + period_size
    period_bits = bitstring[start:end]

    # Extract only important bits within this period
    extracted = ''.join(period_bits[i] for i in important_bits if i < len(period_bits))
    return extracted

def calculate_fidelities(job_id_evolved, job_id_initial, periods=50, period_size=9, shots=1024, filename="fidelity_results.txt"):
    service = QiskitRuntimeService(name="qgss-2025-fun-2")

    job_evolved = service.job(job_id_evolved)
    job_initial = service.job(job_id_initial)

    result_evolved = job_evolved.result()
    result_initial = job_initial.result()

    databin_evolved = result_evolved[0].data
    databin_initial = result_initial[0].data

    counts_evolved = databin_to_counts(databin_evolved)
    counts_initial = databin_to_counts(databin_initial)

    probs_evolved = normalize_counts(counts_evolved, shots)
    probs_initial = normalize_counts(counts_initial, shots)

    fidelities = []

    for n in range(periods):
        period_probs_evolved = {}
        period_probs_initial = {}

        for bitstring, prob in probs_evolved.items():
            # Reverse bitstring to MSB-first order (Qiskit returns LSB-first)
            bitstring_msb = bitstring[::-1]
            key = extract_important_bits_per_period(bitstring_msb, n, period_size, IMPORTANT_BITS)
            period_probs_evolved[key] = period_probs_evolved.get(key, 0) + prob

        for bitstring, prob in probs_initial.items():
            bitstring_msb = bitstring[::-1]
            key = extract_important_bits_per_period(bitstring_msb, n, period_size, IMPORTANT_BITS)
            period_probs_initial[key] = period_probs_initial.get(key, 0) + prob

        fidelity = 0.0
        for outcome, prob in period_probs_evolved.items():
            fidelity += math.sqrt(prob * period_probs_initial.get(outcome, 0))

        fidelities.append(fidelity)
        print(f"Classical fidelity on important bits at period {n + 1}T: {fidelity:.6f}")

    with open(filename, "w") as f:
        for n, fidelity in enumerate(fidelities):
            f.write(f"Fidelity at period {n + 1}T: {fidelity:.6f}\n")

    print(f"\nFidelity results saved to '{filename}'")
    return fidelities


# Example usage:
if __name__ == "__main__":
    evolved_job_id = "d1ltiqna572c7397te00"
    initial_job_id = "d1ltiqv29o4s73aq9lig"

    fidelities = calculate_fidelities(
        job_id_evolved=evolved_job_id,
        job_id_initial=initial_job_id,
        periods=50,
        period_size=9,
        shots=1024,
        filename="fidelity_results.txt"
    )
