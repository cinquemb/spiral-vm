# quick_spiralvm.py
import spiralvm
from quick import Backend, Circuit, GateType

class SpiralVMBackend(Backend):
    def __init__(self, size=30):
        self.size = size
        self.vm = None  # Initialize in run()
    
    def run(self, circuit: Circuit):
        # Fresh VM per circuit (Quick standard)
        self.vm = SpiralVM(self.size, self.size)
        self.vm.is_ang = True
        self.vm.initialize_state("neel")
        
        # Map Quick gates → YOUR apply_gate()
        for gate in circuit.gates:
            g = self.make_gate(gate)  # Convert Quick IR → your Gate struct
            
            if gate.type in [GateType.CX, GateType.X, GateType.Z, GateType.H]:
                self.vm.apply_gate(g, 1.0)  # YOUR dispatcher handles everything
                self.vm.run_periods(1)      # Gate duration
            
            elif gate.type == GateType.measure:
                pass  # Z_L is continuous
        
        self.vm.compile_to_physical_waveform()
        self.vm.run_periods(100)  # Final measurement
        
        # Quick-compatible results format
        n_qubits = len(self.vm.logical_qubits)
        return {
            "Z": [self.vm.measure_logical_Z(i) for i in range(n_qubits)],
            "ZZ": [self.vm.logical_zz_correlation(i,i+1) for i in range(n_qubits-1)]
        }
    
    def make_gate(self, quick_gate):
        # Convert Quick IR → your Gate enum
        if quick_gate.type == GateType.CX:
            return Gate{Gate::CNOT, quick_gate.control, quick_gate.target}
        elif quick_gate.type == GateType.X:
            return Gate{Gate::X, quick_gate.target}
        elif quick_gate.type == GateType.Z:
            return Gate{Gate::Z, quick_gate.target, 0, quick_gate.params[0]}
        elif quick_gate.type == GateType.H:
            return Gate{Gate::H, quick_gate.target}
        # Add more mappings...
