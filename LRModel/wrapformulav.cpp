#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "lrformulav.h"

namespace py = pybind11;

//void init_LRFaxial(py::module_ &m) {
PYBIND11_MODULE(lrformulav, m) {
    py::class_<LRFormulaV>(m, "LRFormulaV")
        // Constructors
        .def(py::init<double, double, double>(),
            py::arg("x0"), py::arg("y0"), py::arg("rmax"))
//        .def(py::init<const Json &>())
        .def(py::init<std::string &>())

        // Clone
        .def("clone",
            [](const LRFormulaV &self) {
                return self.clone();
            },
            py::return_value_policy::take_ownership)

        // Domain & readiness
        .def("inDomain", &LRFormulaV::inDomain,
             py::arg("x"), py::arg("y"), py::arg("z") = 0.0)
        .def("isReady", &LRFormulaV::isReady)

        // Parameters
        .def("getRmax", &LRFormulaV::getRmax)

        // Evaluation
        .def("eval", &LRFormulaV::eval,
             py::arg("x"), py::arg("y"), py::arg("z") = 0.0)
        .def("evalAxial", &LRFormulaV::evalAxial)
        .def("evalDrvX", &LRFormulaV::evalDrvX,
             py::arg("x"), py::arg("y"), py::arg("z") = 0.0)
        .def("evalDrvY", &LRFormulaV::evalDrvY,
             py::arg("x"), py::arg("y"), py::arg("z") = 0.0)

        // Data fitting
        .def("fitData", &LRFormulaV::fitData)
        .def("addData", &LRFormulaV::addData)
        .def("doFit", &LRFormulaV::doFit)
        .def("clearData", &LRFormulaV::clearData)

        // Spline + JSON
//        .def("getSpline",
//             &LRFormulaV::getSpline,
//             py::return_value_policy::reference_internal)
        .def("type", &LRFormulaV::type)
//        .def("ToJsonObject", &LRFormulaV::ToJsonObject)

        // Setters
        .def("SetOrigin", &LRFormulaV::SetOrigin)
        .def("SetRmin", &LRFormulaV::SetRmin)
        .def("SetRmax", &LRFormulaV::SetRmax)
//        .def("SetParameters", &LRFormulaV::SetParameters)

        // Geometry
        .def("R", &LRFormulaV::R)
        .def("R2", &LRFormulaV::R2)

        // Getters
        .def("GetRmin", &LRFormulaV::GetRmin)
        .def("GetOriginX", &LRFormulaV::GetOriginX)
        .def("GetOriginY", &LRFormulaV::GetOriginY)
//        .def("GetParameters", &LRFormulaV::GetParameters)

//  VFormula-related calls
    // setters
        .def("SetParameter", &LRFormulaV::SetParameter)
        .def("SetExpression", &LRFormulaV::SetExpression)

        // call this to update the virtual machine if parameter names of expression were changed
        .def("InitVF", &LRFormulaV::InitVF)

        // getters
        .def("GetParameter", &LRFormulaV::GetParameter)
        .def("GetParVector", &LRFormulaV::GetParVector)

        // VM-debugging
        .def("GetPrg", &LRFormulaV::GetPrg)
        .def("GetConstMap", &LRFormulaV::GetConstMap)
        .def("GetVarMap", &LRFormulaV::GetVarMap)
        .def("GetOperMap", &LRFormulaV::GetOperMap)
        .def("GetFuncMap", &LRFormulaV::GetFuncMap) 
        .def("GetVersion", &LRFormulaV::GetVersion)
        
        .def("GetJsonString", &LRFormulaV::GetJsonString)
//        .def("GetHist",
//             &LRFormulaV::GetHist,
//             py::return_value_policy::reference_internal)

        // Ratio
//        .def("GetRatio", &LRFormulaV::GetRatio)
        ;
}
