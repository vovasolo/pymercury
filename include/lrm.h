#ifndef lrm_h
#define lrm_h

#include <pybind11/pybind11.h>
namespace py = pybind11;

void init_lrm(py::module& m);

#endif