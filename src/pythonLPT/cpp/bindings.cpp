//
// Created by Georgia Stuart on 2/3/21.
//

#ifndef PYTHONLPT_BINDINGS_CPP
#define PYTHONLPT_BINDINGS_CPP

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include "library.h"

namespace py = pybind11;

PYBIND11_MODULE(pythonlpt_cpp, m) {
    m.def("hello", &hello);
}

#endif //PYTHONLPT_BINDINGS_CPP