#include <pybind11/pybind11.h>
#include "descriptor.h"
#include <pybind11/stl.h>


#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

PYBIND11_MODULE(ext, m) {
    py::class_<SomeClass>(m, "SomeClass")
        .def(py::init<float>())
        .def("multiply", &SomeClass::multiply)
        .def("multiply_items", &SomeClass::multiply_items);
}
