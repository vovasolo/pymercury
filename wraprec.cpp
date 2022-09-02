#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "reconstructor.h"
#include "reconstructor_mp.h"

namespace py = pybind11;

PYBIND11_MODULE(mercury, m) {
    py::class_<RecCoG>(m, "RecCoG")
        .def(py::init<std::string &>())

        .def("ProcessEvent", &RecCoG::ProcessEvent)
    
        .def("getRecX", &RecCoG::getRecX)
        .def("getRecY", &RecCoG::getRecY)

        .def("setCogAbsCutoff", &RecCoG::setCogAbsCutoff)
        .def("setCogRelCutoff", &RecCoG::setCogRelCutoff)
        .def("getSumSignal", &RecCoG::getSumSignal)
        ;

    py::class_<RecLS>(m, "RecLS")
        .def(py::init<std::string &, bool>())

        .def("ProcessEvent", &RecLS::ProcessEvent)

        .def("getGuessX", &RecLS::getGuessX)
        .def("getGuessY", &RecLS::getGuessY)
        .def("getGuessE", &RecLS::getGuessE)
        .def("getRecStatus", &RecLS::getRecStatus)
        .def("getDof", &RecLS::getDof)        
        .def("getRecX", &RecLS::getRecX)
        .def("getRecY", &RecLS::getRecY)
        .def("getRecE", &RecLS::getRecE)        
        .def("getRecMin", &RecLS::getRecMin)
        .def("getChi2Min", &RecLS::getChi2Min)
        .def("getCovXX", &RecLS::getCovXX)
        .def("getCovYY", &RecLS::getCovYY)
        .def("getCovXY", &RecLS::getCovXY)

        .def("setCogAbsCutoff", &RecLS::setCogAbsCutoff)
        .def("setCogRelCutoff", &RecLS::setCogRelCutoff)
        .def("setRecAbsCutoff", &RecLS::setRecAbsCutoff)
        .def("setRecRelCutoff", &RecLS::setRecRelCutoff)
        .def("setRecCutoffRadius", &RecLS::setRecCutoffRadius)
        .def("setEnergyCalibration", &RecLS::setEnergyCalibration)   

        .def("setRMstepX", &RecLS::setRMstepX)
        .def("setRMstepY", &RecLS::setRMstepY)
        .def("setRMstepEnergy", &RecLS::setRMstepEnergy)
        .def("setRMmaxFuncCalls", &RecLS::setRMmaxFuncCalls)
        .def("setRMmaxIterations", &RecLS::setRMmaxIterations)
        .def("setRMtolerance", &RecLS::setRMtolerance)

        .def("setMinuitVerbosity", &RecLS::setMinuitVerbosity)

        .def("setGuessPosition", &RecLS::setGuessPosition)
        .def("setGuessPositionAuto", &RecLS::setGuessPositionAuto)
        .def("setGuessEnergy", &RecLS::setGuessEnergy)
        .def("setGuessEnergyAuto", &RecLS::setGuessEnergyAuto)

        .def("getSumSignal", &RecLS::getSumSignal)
        .def("getSumLRF", &RecLS::getSumLRF)
        .def("getSumActiveSignal", &RecLS::getSumActiveSignal)
        .def("getSumActiveLRF", &RecLS::getSumActiveLRF)

        .def("setAutoE", &RecLS::setAutoE) 
        ;

    py::class_<RecML>(m, "RecML")
        .def(py::init<std::string &>())

        .def("ProcessEvent", &RecML::ProcessEvent)

        .def("getGuessX", &RecML::getGuessX)
        .def("getGuessY", &RecML::getGuessY)
        .def("getGuessE", &RecML::getGuessE)
        .def("getRecStatus", &RecML::getRecStatus)
        .def("getDof", &RecML::getDof)        
        .def("getRecX", &RecML::getRecX)
        .def("getRecY", &RecML::getRecY)
        .def("getRecE", &RecML::getRecE)        
        .def("getRecMin", &RecML::getRecMin)
        .def("getChi2Min", &RecML::getChi2Min)
        .def("getCovXX", &RecML::getCovXX)
        .def("getCovYY", &RecML::getCovYY)
        .def("getCovXY", &RecML::getCovXY)

        .def("setCogAbsCutoff", &RecML::setCogAbsCutoff)
        .def("setCogRelCutoff", &RecML::setCogRelCutoff)
        .def("setRecAbsCutoff", &RecML::setRecAbsCutoff)
        .def("setRecRelCutoff", &RecML::setRecRelCutoff)
        .def("setRecCutoffRadius", &RecML::setRecCutoffRadius)
        .def("setEnergyCalibration", &RecML::setEnergyCalibration)   

        .def("setRMstepX", &RecML::setRMstepX)
        .def("setRMstepY", &RecML::setRMstepY)
        .def("setRMstepEnergy", &RecML::setRMstepEnergy)
        .def("setRMmaxFuncCalls", &RecML::setRMmaxFuncCalls)
        .def("setRMmaxIterations", &RecML::setRMmaxIterations)
        .def("setRMtolerance", &RecML::setRMtolerance)

        .def("setMinuitVerbosity", &RecML::setMinuitVerbosity)

        .def("setGuessPosition", &RecML::setGuessPosition)
        .def("setGuessPositionAuto", &RecML::setGuessPositionAuto)
        .def("setGuessEnergy", &RecML::setGuessEnergy)
        .def("setGuessEnergyAuto", &RecML::setGuessEnergyAuto)

        .def("getSumSignal", &RecML::getSumSignal)
        .def("getSumLRF", &RecML::getSumLRF)
        .def("getSumActiveSignal", &RecML::getSumActiveSignal)
        .def("getSumActiveLRF", &RecML::getSumActiveLRF)

        .def("setAutoE", &RecLS::setAutoE) 
        ;

        py::class_<ReconstructorMP>(m, "RecCoG_MP")
        .def(py::init<std::string &, int>())

        .def("ProcessEvents", &ReconstructorMP::ProcessEvents)
        
        .def("getRecX", &ReconstructorMP::getRecX)
        .def("getRecY", &ReconstructorMP::getRecY)

        .def("setCogAbsCutoff", &ReconstructorMP::setCogAbsCutoff)
        .def("setCogRelCutoff", &ReconstructorMP::setCogRelCutoff) 
        ;

        py::class_<RecLS_MP>(m, "RecLS_MP")
        .def(py::init<std::string &, int, bool>())

        .def("ProcessEvents", &RecLS_MP::ProcessEvents)

        .def("getRecStatus", &RecLS_MP::getRecStatus)
        .def("getDof", &RecLS_MP::getDof)        
        .def("getRecX", &RecLS_MP::getRecX)
        .def("getRecY", &RecLS_MP::getRecY)
        .def("getRecE", &RecLS_MP::getRecE)        
        .def("getRecMin", &RecLS_MP::getRecMin)
        .def("getChi2Min", &RecLS_MP::getChi2Min)
        .def("getCovXX", &RecLS_MP::getCovXX)
        .def("getCovYY", &RecLS_MP::getCovYY)
        .def("getCovXY", &RecLS_MP::getCovXY)

        .def("setCogAbsCutoff", &RecLS_MP::setCogAbsCutoff)
        .def("setCogRelCutoff", &RecLS_MP::setCogRelCutoff)
        .def("setRecAbsCutoff", &RecLS_MP::setRecAbsCutoff)
        .def("setRecRelCutoff", &RecLS_MP::setRecRelCutoff)
        .def("setRecCutoffRadius", &RecLS_MP::setRecCutoffRadius)
        .def("setEnergyCalibration", &RecLS_MP::setEnergyCalibration) 

        .def("setRMstepX", &RecLS_MP::setRMstepX)
        .def("setRMstepY", &RecLS_MP::setRMstepY)
        .def("setRMstepEnergy", &RecLS_MP::setRMstepEnergy)
        .def("setRMmaxFuncCalls", &RecLS_MP::setRMmaxFuncCalls)
        .def("setRMmaxIterations", &RecLS_MP::setRMmaxIterations)
        .def("setRMtolerance", &RecLS_MP::setRMtolerance)

        .def("setMinuitVerbosity", &RecLS_MP::setMinuitVerbosity)

        .def("setAutoE", &RecLS_MP::setAutoE) 
        ; 

        py::class_<RecML_MP>(m, "RecML_MP")
        .def(py::init<std::string &, int>())

        .def("ProcessEvents", &RecML_MP::ProcessEvents)

        .def("getRecStatus", &RecML_MP::getRecStatus)
        .def("getDof", &RecML_MP::getDof)        
        .def("getRecX", &RecML_MP::getRecX)
        .def("getRecY", &RecML_MP::getRecY)
        .def("getRecE", &RecML_MP::getRecE)        
        .def("getRecMin", &RecML_MP::getRecMin)
        .def("getChi2Min", &RecML_MP::getChi2Min)
        .def("getCovXX", &RecML_MP::getCovXX)
        .def("getCovYY", &RecML_MP::getCovYY)
        .def("getCovXY", &RecML_MP::getCovXY)

        .def("setCogAbsCutoff", &RecML_MP::setCogAbsCutoff)
        .def("setCogRelCutoff", &RecML_MP::setCogRelCutoff)
        .def("setRecAbsCutoff", &RecML_MP::setRecAbsCutoff)
        .def("setRecRelCutoff", &RecML_MP::setRecRelCutoff)
        .def("setRecCutoffRadius", &RecML_MP::setRecCutoffRadius)
        .def("setEnergyCalibration", &RecML_MP::setEnergyCalibration) 

        .def("setRMstepX", &RecML_MP::setRMstepX)
        .def("setRMstepY", &RecML_MP::setRMstepY)
        .def("setRMstepEnergy", &RecML_MP::setRMstepEnergy)
        .def("setRMmaxFuncCalls", &RecML_MP::setRMmaxFuncCalls)
        .def("setRMmaxIterations", &RecML_MP::setRMmaxIterations)
        .def("setRMtolerance", &RecML_MP::setRMtolerance)

        .def("setMinuitVerbosity", &RecML_MP::setMinuitVerbosity)

        .def("setAutoE", &RecML_MP::setAutoE) 
        ;
    
}
