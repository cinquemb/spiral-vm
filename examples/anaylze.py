import numpy as np

def load_log(path):
    steps = []
    Z = []
    X = []
    with open(path, "r") as f:
        for line in f:
            if "|" not in line:
                continue
            parts = line.split("|")
            try:
                step = int(parts[0].strip())
                z = float(parts[1])
                x = float(parts[2])
            except:
                continue
            steps.append(step)
            Z.append(z)
            X.append(x)
    return np.array(steps), np.array(Z), np.array(X)

def analyze(path, window=50):
    steps, Z, X = load_log(path)

    dZ = Z[1:] - Z[:-1]
    dX = X[1:] - X[:-1]

    results = []

    for i in range(len(Z) - window - 1):
        Zw = Z[i:i+window]
        Xw = X[i:i+window]
        dZw = dZ[i:i+window]
        dXw = dX[i:i+window]

        M = np.vstack([Zw, Xw]).T
        Y = np.vstack([dZw, dXw]).T

        # least squares fit
        A, _, _, _ = np.linalg.lstsq(M, Y, rcond=None)
        # A maps [Z, X] -> [dZ, dX]
        J = np.eye(2) + A  # discrete-time Jacobian

        eigvals = np.linalg.eigvals(J)
        spectral_radius = max(abs(eigvals))

        results.append((steps[i], J, eigvals, spectral_radius))

    # Find region where system is near neutral (saddle-like)
    results.sort(key=lambda r: abs(r[3] - 1.0))
    step, J, eig, rho = results[0]

    print("Closest saddle-like region around step:", step)
    print("Jacobian:")
    print(J)
    print("Eigenvalues:", eig)
    print("Spectral radius:", rho)

    # Simple linear feedback to drive toward (0,1)
    target = np.array([0.0, 1.0])
    # We want: [Z,X]_{t+1} = J*[Z,X] + u  ->  target
    # u = target - J*target
    u = target - J @ target
    print("\nSuggested constant bias (dZ, dX):", u)
    print("Interpretation:")
    print("  logical_z_rotation ≈ u[0]")
    print("  logical_x_pulse    ≈ u[1]")

if __name__ == "__main__":
    analyze("run.log", window=50)
