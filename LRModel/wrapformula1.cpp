#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "lrformula1.h"

namespace py = pybind11;

//void init_LRFaxial(py::module_ &m) {
PYBIND11_MODULE(lrformula1, m) {
    py::class_<LRFormula1>(m, "LRFormula1")
        // Constructors
        .def(py::init<double, double, double>(),
            py::arg("x0"), py::arg("y0"), py::arg("rmax"))
//        .def(py::init<const Json &>())
        .def(py::init<std::string &>())

        // Clone
        .def("clone",
            [](const LRFormula1 &self) {
                return self.clone();
            },
            py::return_value_policy::take_ownership)

        // Domain & readiness
        .def("inDomain", &LRFormula1::inDomain,
             py::arg("x"), py::arg("y"), py::arg("z") = 0.0)
        .def("isReady", &LRFormula1::isReady)

        // Parameters
        .def("getRmax", &LRFormula1::getRmax)

        // Evaluation
        .def("eval", &LRFormula1::eval,
             py::arg("x"), py::arg("y"), py::arg("z") = 0.0)
        .def("evalAxial", &LRFormula1::evalAxial)
        .def("evalDrvX", &LRFormula1::evalDrvX,
             py::arg("x"), py::arg("y"), py::arg("z") = 0.0)
        .def("evalDrvY", &LRFormula1::evalDrvY,
             py::arg("x"), py::arg("y"), py::arg("z") = 0.0)

        // Data fitting
        .def("fitData", &LRFormula1::fitData)
        .def("addData", &LRFormula1::addData)
        .def("doFit", &LRFormula1::doFit)
        .def("clearData", &LRFormula1::clearData)

        // Spline + JSON
//        .def("getSpline",
//             &LRFormula1::getSpline,
//             py::return_value_policy::reference_internal)
        .def("type", &LRFormula1::type)
//        .def("ToJsonObject", &LRFormula1::ToJsonObject)

        // Setters
        .def("SetOrigin", &LRFormula1::SetOrigin)
        .def("SetRmin", &LRFormula1::SetRmin)
        .def("SetRmax", &LRFormula1::SetRmax)
        .def("SetParameters", &LRFormula1::SetParameters)

        // Geometry
        .def("R", &LRFormula1::R)
        .def("R2", &LRFormula1::R2)

        // Getters
        .def("GetRmin", &LRFormula1::GetRmin)
        .def("GetOriginX", &LRFormula1::GetOriginX)
        .def("GetOriginY", &LRFormula1::GetOriginY)
        .def("GetParameters", &LRFormula1::GetParameters)

//        .def("GetJsonString", &LRFormula1::GetJsonString)
//        .def("GetHist",
//             &LRFormula1::GetHist,
//             py::return_value_policy::reference_internal)

        // Ratio
//        .def("GetRatio", &LRFormula1::GetRatio)
        ;
}
