#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "lrformulaxy.h"

namespace py = pybind11;

//void init_LRFaxial(py::module_ &m) {
PYBIND11_MODULE(lrformulaxy, m) {
    py::class_<LRFormulaXY>(m, "LRFormulaXY")
        // Constructors
        .def(py::init<double, double, double, double>(),
            py::arg("xmin"), py::arg("xmax"), py::arg("ymin"), py::arg("ymax"))
//        .def(py::init<const Json &>())
        .def(py::init<std::string &>())

        // Clone
        .def("clone",
            [](const LRFormulaXY &self) {
                return self.clone();
            },
            py::return_value_policy::take_ownership)

        // Domain & readiness
        .def("inDomain", &LRFormulaXY::inDomain,
             py::arg("x"), py::arg("y"), py::arg("z") = 0.0)
        .def("isReady", &LRFormulaXY::isReady)

        // Parameters
        .def("getRmax", &LRFormulaXY::getRmax)

        // Evaluation
        .def("eval", &LRFormulaXY::eval,
             py::arg("x"), py::arg("y"), py::arg("z") = 0.0)

        .def("evalDrvX", &LRFormulaXY::evalDrvX,
             py::arg("x"), py::arg("y"), py::arg("z") = 0.0)
        .def("evalDrvY", &LRFormulaXY::evalDrvY,
             py::arg("x"), py::arg("y"), py::arg("z") = 0.0)

        // Data fitting
        .def("fitData", &LRFormulaXY::fitData)
        .def("addData", &LRFormulaXY::addData)
        .def("addDataPy", &LRFormulaXY::addDataPy)
        .def("doFit", &LRFormulaXY::doFit)
        .def("clearData", &LRFormulaXY::clearData)
        .def_readwrite("qtol", &LRFormulaXY::qtol)
        .def_readwrite("fvec", &LRFormulaXY::fvec)
        .def_readwrite("fjac", &LRFormulaXY::fjac)

        // Spline + JSON
//        .def("getSpline",
//             &LRFormulaXY::getSpline,
//             py::return_value_policy::reference_internal)
        .def("type", &LRFormulaXY::type)
//        .def("ToJsonObject", &LRFormulaXY::ToJsonObject)

//  VFormula-related calls
    // setters
        .def("SetParameter", &LRFormulaXY::SetParameter)
        .def("SetExpression", &LRFormulaXY::SetExpression)

        // call this to update the virtual machine if parameter names of expression were changed
        .def("InitVF", &LRFormulaXY::InitVF)

        // getters
        .def("GetParameter", &LRFormulaXY::GetParameter)
        .def("GetParVector", &LRFormulaXY::GetParVector)

        // VM-debugging
        .def("GetPrg", &LRFormulaXY::GetPrg)
        .def("GetConstMap", &LRFormulaXY::GetConstMap)
        .def("GetVarMap", &LRFormulaXY::GetVarMap)
        .def("GetOperMap", &LRFormulaXY::GetOperMap)
        .def("GetFuncMap", &LRFormulaXY::GetFuncMap) 
        .def("GetVersion", &LRFormulaXY::GetVersion)
        
        .def("GetJsonString", &LRFormulaXY::GetJsonString)
//        .def("GetHist",
//             &LRFormulaXY::GetHist,
//             py::return_value_policy::reference_internal)

        // Ratio
//        .def("GetRatio", &LRFormulaXY::GetRatio)
        ;
}
