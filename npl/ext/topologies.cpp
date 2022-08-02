#include <pybind11/pybind11.h>

namespace py = pybind11;

int get_bonds_count(int i, int j) {
    return i + j;
}