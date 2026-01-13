// spiralvm_pybind.cpp
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>  // automatic std::vector support
#include "../src/spiral_vm_core.hpp"

namespace py = pybind11;

PYBIND11_MODULE(spiralvm, m) {
    m.doc() = "SpiralVM: Global-drive Floquet logical qubits";
    
    // Export SpiralVM class
    py::class_<SpiralVM>(m, "SpiralVM")
        .def(py::init<int, int>())
        .def("initialize_state", &SpiralVM::initialize_state)
        .def("add_qubit", &SpiralVM::add_qubit)
        .def("apply_gate", &SpiralVM::apply_gate)
        .def("compile_to_physical_waveform", &SpiralVM::compile_to_physical_waveform)
        .def("run_periods", &SpiralVM::run_periods)
        .def("measure_logical_Z", &SpiralVM::measure_logical_Z)
        .def("logical_zz_correlation", &SpiralVM::logical_zz_correlation)
        .def("dump_frequency_mapping", &SpiralVM::dump_frequency_mapping)
        // Add more methods as needed
        ;
    
    // Export Gate enum/struct
    py::enum_<Gate::Type>(m, "GateType")
        .value("X", Gate::X)
        .value("Z", Gate::Z)
        .value("CZ", Gate::CZ)
        .value("CNOT", Gate::CNOT)
        .value("H", Gate::H)
        .value("MEASURE", Gate::MEASURE)
        .value("PHASE", Gate::PHASE);
}

Gate make_gate_pybind(Gate::Type type, uint32_t target=0, uint32_t control=0, double angle=0.0) {
    return Gate{type, target, control, angle};
}
