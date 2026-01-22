# Find pybind11
find_package(pybind11 REQUIRED)

# Python module
pybind11_add_module(spiralvm MODULE spiralvm_pybind.cpp spiral_vm_core.cpp)
target_link_libraries(spiralvm PRIVATE armadillo lapack blas)
target_compile_options(spiralvm PRIVATE -O2)
