#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "lrfaxial.h"

namespace py = pybind11;

//void init_LRFaxial(py::module_ &m) {
PYBIND11_MODULE(lrfaxial, m) {
    py::class_<LRFaxial>(m, "LRFaxial")
        // Constructors
        .def(py::init<double, int>(),
             py::arg("rmax"), py::arg("nint"))
//        .def(py::init<const Json &>())
        .def(py::init<std::string &>())

        // Clone
        .def("clone",
            [](const LRFaxial &self) {
                return self.clone();
            },
            py::return_value_policy::take_ownership)

        // Domain & readiness
        .def("inDomain", &LRFaxial::inDomain,
             py::arg("x"), py::arg("y"), py::arg("z") = 0.0)
        .def("isReady", &LRFaxial::isReady)

        // Parameters
        .def("getRmax", &LRFaxial::getRmax)
        .def("getNint", &LRFaxial::getNint)
        .def("GetNodes", &LRFaxial::GetNodes)

        // Evaluation
        .def("eval", &LRFaxial::eval,
             py::arg("x"), py::arg("y"), py::arg("z") = 0.0)
        .def("evalraw", &LRFaxial::evalraw,
             py::arg("x"), py::arg("y"), py::arg("z") = 0.0)
        .def("evalAxial", &LRFaxial::evalAxial)
        .def("evalDrvX", &LRFaxial::evalDrvX,
             py::arg("x"), py::arg("y"), py::arg("z") = 0.0)
        .def("evalDrvY", &LRFaxial::evalDrvY,
             py::arg("x"), py::arg("y"), py::arg("z") = 0.0)

        // Data fitting
        .def("fitData", &LRFaxial::fitData)
        .def("addData", &LRFaxial::addData)
        .def("doFit", &LRFaxial::doFit)
        .def("clearData", &LRFaxial::clearData)

        // Spline + JSON
//        .def("getSpline",
//             &LRFaxial::getSpline,
//             py::return_value_policy::reference_internal)
        .def("type", &LRFaxial::type)
//        .def("ToJsonObject", &LRFaxial::ToJsonObject)

        // Setters
        .def("SetOrigin", &LRFaxial::SetOrigin)
        .def("SetRmin", &LRFaxial::SetRmin)
        .def("SetRmax", &LRFaxial::SetRmax)
//        .def("SetCompression", &LRFaxial::SetCompression)
//        .def("SetSpline", &LRFaxial::SetSpline)
        .def("SetFlatTop", &LRFaxial::SetFlatTop)
        .def("SetNonIncreasing", &LRFaxial::SetNonIncreasing)

        // Geometry
        .def("R", &LRFaxial::R)
        .def("R2", &LRFaxial::R2)

        .def("Rho", (double (LRFaxial::*) (double) const ) &LRFaxial::Rho)
        .def("Rho_xy", (double (LRFaxial::*) (double, double) const) &LRFaxial::Rho)
        .def("RhoDrvX", &LRFaxial::RhoDrvX)
        .def("RhoDrvY", &LRFaxial::RhoDrvY)

        // Getters
        .def("GetRmin", &LRFaxial::GetRmin)
        .def("GetOriginX", &LRFaxial::GetOriginX)
        .def("GetOriginY", &LRFaxial::GetOriginY)
//        .def("GetHist",
//             &LRFaxial::GetHist,
//             py::return_value_policy::reference_internal)

        // Ratio
//        .def("GetRatio", &LRFaxial::GetRatio)
        ;
}
