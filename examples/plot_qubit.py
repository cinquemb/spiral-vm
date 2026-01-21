import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
import os

def plot_qubit_data(files):
    """Plot N qubit data files - 4th column = X_norm (not Y!)"""
    dfs = []
    colors = plt.cm.tab10(np.linspace(0, 1, len(files)))
    
    # Parse all files
    for i, file in enumerate(files):
        data = []
        with open(file, 'r') as f:
            for line in f:
                parts = line.strip().split('|')
                if len(parts) >= 4:
                    step = int(parts[0].strip())
                    Z, X, X_norm = map(float, [p.strip() for p in parts[1:4]])
                    data.append([step, Z, X, X_norm])
        
        # FIXED: 4 columns = step, Z, X, X_norm
        df = pd.DataFrame(data, columns=['step', 'Z', 'X', 'X_norm'])
        dfs.append(df)
        print(f"Loaded {file}: {len(df)} steps")
    
    # Plotting
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(14, 10))
    
    # Z evolution
    for i, df in enumerate(dfs):
        ax1.plot(df['step'], df['Z'], color=colors[i], label=os.path.basename(files[i]), linewidth=2)
    ax1.set_ylabel('Z')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax1.grid(True, alpha=0.3)
    ax1.set_title('Qubit State Evolution')
    
    # X evolution  
    for i, df in enumerate(dfs):
        ax2.plot(df['step'], df['X'], color=colors[i], label=os.path.basename(files[i]), linewidth=2)
    ax2.set_ylabel('X')
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax2.grid(True, alpha=0.3)
    
    # X_norm evolution (was incorrectly labeled Y!)
    for i, df in enumerate(dfs):
        ax3.plot(df['step'], df['X_norm'], color=colors[i], label=os.path.basename(files[i]), linewidth=2)
    ax3.set_ylabel('X_norm')
    ax3.set_xlabel('Step')
    ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('qubit_evolution.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Bloch sphere projection (X, ?, Z) - using X_norm as Y for visualization
    fig2, ax4 = plt.subplots(1, 1, figsize=(10, 8), subplot_kw={'projection': '3d'})
    for i, df in enumerate(dfs):
        # X, X_norm, Z for Bloch sphere (since no true Y)
        ax4.plot(df['X'], df['X_norm'], df['Z'], color=colors[i], label=os.path.basename(files[i]), linewidth=2)
        ax4.scatter(df['X'], df['X_norm'], df['Z'], c=df['step'], cmap=f'viridis_{i}', s=1, alpha=0.6)
    ax4.set_xlabel('X')
    ax4.set_ylabel('X_norm')
    ax4.set_zlabel('Z')
    ax4.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.title('Bloch Sphere Trajectories (X, X_norm, Z)')
    plt.savefig('bloch_sphere.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Print final states
    print("\nFinal states:")
    print("-" * 60)
    print("File            | Steps |   Z      X     X_norm")
    print("-" * 60)
    for i, (file, df) in enumerate(zip(files, dfs)):
        final = df.iloc[-1]
        print(f"{os.path.basename(file):15} | {len(df):5} | {final.Z:7.4f} {final.X:7.4f} {final.X_norm:7.4f}")

# Usage
if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        plot_qubit_data(sys.argv[1:])
    else:
        # Example usage
        #plot_qubit_data(['no_q_drive.txt', '2q_drive.txt', '1q_drive.txt', "degen_q_drive.txt", "degen1_q_drive.txt", "degen2_q_drive.txt","degen3_q_drive.txt","degen4_q_drive.txt","degen5_q_drive.txt","hadamard_drive.txt", "hadamard_drive_K20.txt","opt_q_drive.txt"])
        plot_qubit_data(["degen_q_drive.txt","hadamard_drive.txt", "hadamard_drive_K20.txt","opt_q_drive.txt"])
        print("\nOr run: python plot_qubit.py *.txt")
