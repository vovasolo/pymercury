#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "lrmodel.h"
#include "bspline123d.h"
#include "bsfit123.h"

namespace py = pybind11;

PYBIND11_MODULE(lrmodel, m) {
    py::class_<LRModel>(m, "LRModel")
        .def(py::init<int>())
        .def(py::init<std::string &>())
        .def("Version", &LRModel::Version)
        .def("SetDisabled", &LRModel::SetDisabled,
            "Disables a sensor")
        .def("SetEnabled", &LRModel::SetEnabled,
            "Enables a sensor")
        .def("IsDisabled", &LRModel::IsDisabled,
            "Returns True if the sensor is disabled")
        .def("SetDisabledList", &LRModel::SetDisabledList,
            "Disables all sensors in the list")
        .def("SetEnabledList", &LRModel::SetEnabledList,
            "Enables all sensors in the list")
        .def("GetDisabledList", &LRModel::GetDisabledList,
            "Returns list of all disabled sensors")
        .def("GetEnabledList", &LRModel::GetEnabledList,
            "Returns list of all enabled sensors")
        .def("SetGain", &LRModel::SetGain,
            "Sets gain for a sensor",
            py::arg("id"), py::arg("gain"))
        .def("SetGroup", &LRModel::SetGroup,
            "Sets group for a sensor",
            py::arg("id"), py::arg("gid"))
        .def("SetTransform", (bool (LRModel::*) (int, std::string &)) &LRModel::SetTransform,
            "Sets transform for a sensor from provided JSON string")
        .def("GetGain", &LRModel::GetGain)
        .def("GetGroup", &LRModel::GetGroup)
        .def("GetTransform", &LRModel::GetJsonTransform)
        .def("GetX", &LRModel::GetX)
        .def("GetY", &LRModel::GetY)
        .def("GetAllX", &LRModel::GetAllX)
        .def("GetAllY", &LRModel::GetAllY)
        .def("GetGroupX", &LRModel::GetGroupX)
        .def("GetGroupY", &LRModel::GetGroupY)
        .def("GetRadius", &LRModel::GetRadius)
        .def("GetPhi", &LRModel::GetPhi)
        .def("GetDistance", &LRModel::GetDistance)

        .def("GroupMembers", &LRModel::GroupMembers)
        .def("GetSensorCount", &LRModel::GetSensorCount)
        .def("GetGroupCount", &LRModel::GetGroupCount)
        .def("GetGroupMembersCount", &LRModel::GetGroupMembersCount)
        .def("SensorExists", &LRModel::SensorExists)
        .def("GroupExists", &LRModel::GroupExists)
        .def("ClearAll", &LRModel::ClearAll)
        .def("AddSensor", &LRModel::AddSensor)

// low lovel group control
        .def("ResetGroups", &LRModel::ResetGroups)
        .def("CreateGroup", &LRModel::CreateGroup)
        .def("DissolveGroup", &LRModel::DissolveGroup)
        .def("AddToGroup", (bool (LRModel::*) (int, int, std::string &)) &LRModel::AddToGroup)
        .def("RemoveFromGroup", &LRModel::RemoveFromGroup)
// high level group control
        .def("SetTolerance", &LRModel::SetTolerance)
        .def("MakeGroupsCommon", &LRModel::MakeGroupsCommon)
        .def("MakeGroupsByRadius", &LRModel::MakeGroupsByRadius)
        .def("MakeGroupsRectangle", &LRModel::MakeGroupsRectangle)
        .def("MakeGroupsSquare", &LRModel::MakeGroupsSquare)
        .def("MakeGroupsHexagon", &LRModel::MakeGroupsHexagon)
        .def("MakeGroupsNgon", &LRModel::MakeGroupsNgon)

// set & get LRFs (passed as JSON strings) 
        .def("SetLRF", &LRModel::SetJsonLRF)
        .def("GetLRF", &LRModel::GetJsonLRF)
        .def("SetGroupLRF", &LRModel::SetGroupJsonLRF)
        .def("GetGroupLRF", &LRModel::GetGroupJsonLRF)      
        .def("SetDefaultLRF", &LRModel::SetDefaultJsonLRF) 

// access to profile histograms in the fitters
        .def("GetHistBins", &LRModel::GetHistBins)
        .def("GetHistNdim", &LRModel::GetHistNdim)
        .def("GetHistBinsX", &LRModel::GetHistBinsX)        
        .def("GetHistBinsY", &LRModel::GetHistBinsY) 
        .def("GetHistBinsZ", &LRModel::GetHistBinsZ)
        .def("GetHistMeans", &LRModel::GetHistMeans)
        .def("GetHistSigmas", &LRModel::GetHistSigmas)        
        .def("GetHistWeights", &LRModel::GetHistWeights) 
        .def("GetHistXCenters", &LRModel::GetHistXCenters)
        .def("GetHistYCenters", &LRModel::GetHistYCenters)
        .def("GetHistZCenters", &LRModel::GetHistZCenters)
        .def("GetHistFit", &LRModel::GetHistFit)                       

// evaluation
        .def("InDomain", (bool (LRModel::*) (int, double, double, double)) &LRModel::InDomain)
        .def("Eval", (double (LRModel::*) (int, double, double, double)) &LRModel::Eval)
	    .def("Eval", (std::vector<double> (LRModel::*) (int, const std::vector <double> &, const std::vector <double> &, const std::vector <double> &)) &LRModel::Eval)
        .def("EvalDrvX", (double (LRModel::*) (int, double, double, double)) &LRModel::EvalDrvX)
        .def("EvalDrvY", (double (LRModel::*) (int, double, double, double)) &LRModel::EvalDrvY)
        .def("InDomainAll", &LRModel::InDomainAll)
        .def("EvalAll", &LRModel::EvalAll)
        .def("EvalDrvXAll", &LRModel::EvalDrvXAll)
        .def("EvalDrvYAll", &LRModel::EvalDrvYAll)
        .def("InDomainList", &LRModel::InDomainList)
        .def("EvalList", &LRModel::EvalList)
        .def("EvalDrvXList", &LRModel::EvalDrvXList)
        .def("EvalDrvYList", &LRModel::EvalDrvYList) 
        .def("EvalAxial", &LRModel::EvalAxial)       
// fitting
        .def("FitNotBinnedData", &LRModel::FitNotBinnedData)
//        .def("AddFitData", &LRModel::AddFitData, py::arg(), py::arg(), py::arg("ignore_group") = false)
        .def("AddFitData", &LRModel::AddFitData)
        .def("FitSensor", &LRModel::FitSensor)
        .def("FitGroup", &LRModel::FitGroup)
        .def("ClearSensorFitData", &LRModel::ClearSensorFitData)
        .def("ClearGroupFitData", &LRModel::ClearGroupFitData)
        .def("ClearAllFitData", &LRModel::ClearAllFitData)

// Sum of LRFs (used in light collection correction)
        .def("GetLRFSum", 
            (double (LRModel::*) (double, double, std::vector <bool> *)) &LRModel::GetLRFSum, 
            "Get sum of LRFs at the point (x,y), with optional boolean exclusion mask (True means excluded)",
            py::arg("x"), py::arg("y"), py::arg("mask") = new std::vector <bool>)
        .def("GetLRFSum", 
            (double (LRModel::*) (double, double, double, std::vector <bool> *)) &LRModel::GetLRFSum, 
            "Get sum of LRFs at the point (x,y,z), with optional boolean exclusion mask (True means excluded)",
            py::arg("x"), py::arg("y"), py::arg("z"), py::arg("mask") = new std::vector <bool>)
        .def("GetLRFSum", 
            (std::vector <double> (LRModel::*) (std::vector <double> &, std::vector <double> &, std::vector<std::vector <bool> > *)) &LRModel::GetLRFSum, 
            "Get sum of LRFs at the point (x,y), with optional boolean exclusion mask (True means excluded) - vector version",
            py::arg("x"), py::arg("y"), py::arg("mask") = new std::vector<std::vector <bool> >)
        .def("GetLRFSum", 
            (std::vector <double> (LRModel::*) (std::vector <double> &, std::vector <double> &, std::vector <double> &, std::vector<std::vector <bool> > *)) &LRModel::GetLRFSum, 
            "Get sum of LRFs at the point (x,y,z), with optional boolean exclusion mask (True means excluded) - vector version",
            py::arg("x"), py::arg("y"), py::arg("z"), py::arg("mask") = new std::vector<std::vector <bool> >)
        .def("GetLRFSumList", &LRModel::GetLRFSumList)

// Light collection correction factor
        .def("SetRefPoint", 
            &LRModel::SetRefPoint, 
            "Set reference position for calculation of light correction factors",
            py::arg("x"), py::arg("y"), py::arg("z") = 0)
        .def("GetCorrFactor", 
            (double (LRModel::*) (double, double, std::vector <bool> *)) &LRModel::GetCorrFactor, 
            "Get correction factor at the point (x,y), with optional boolean exclusion mask (True means excluded)",
            py::arg("x"), py::arg("y"), py::arg("mask") = new std::vector <bool>)
        .def("GetCorrFactor", 
            (double (LRModel::*) (double, double, double, std::vector <bool> *)) &LRModel::GetCorrFactor, 
            "Get correction factor at the point (x,y,z), with optional boolean exclusion mask (True means excluded)",
            py::arg("x"), py::arg("y"), py::arg("z"), py::arg("mask") = new std::vector <bool>)
        .def("GetCorrFactor", 
            (std::vector <double> (LRModel::*) (std::vector <double> &, std::vector <double> &, std::vector<std::vector <bool> > *)) &LRModel::GetCorrFactor, 
            "Get correction factors at the point (x,y), with optional boolean exclusion mask (True means excluded) - vector version",
            py::arg("x"), py::arg("y"), py::arg("mask") = new std::vector<std::vector <bool> >)
        .def("GetCorrFactor", 
            (std::vector <double> (LRModel::*) (std::vector <double> &, std::vector <double> &, std::vector <double> &, std::vector<std::vector <bool> > *)) &LRModel::GetCorrFactor, 
            "Get correction factors at the point (x,y,z), with optional boolean exclusion mask (True means excluded) - vector version",
            py::arg("x"), py::arg("y"), py::arg("z"), py::arg("mask") = new std::vector<std::vector <bool> >)
        .def("GetCorrFactorList", &LRModel::GetCorrFactorList)

// Compute area corrected for light collection
        .def("GetCorrectedArea", (double (LRModel::*) (std::vector <double> &, std::vector <bool> &, double, double, double)) &LRModel::GetCorrectedArea)
        .def("GetCorrectedArea", (std::vector <double> (LRModel::*) (std::vector<std::vector <double> > &, std::vector<std::vector <bool> > &, std::vector <double> &, std::vector <double> &, std::vector <double> &)) &LRModel::GetCorrectedArea)

// export/save        
        .def("ModelAsJson", &LRModel::GetJsonString)
        .def("SensorAsJson", &LRModel::SensorGetJsonString)
        .def("GroupAsJson", &LRModel::GroupGetJsonString)

// reporting 
        .def("GetLRFError", &LRModel::GetLRFError)
        .def("GetError", &LRModel::GetError)
        ;
        
// gain estimator
    py::class_<GainEstimator>(m, "GainEstimator")
        .def(py::init<std::string &>())
        .def("AddData", &GainEstimator::AddData)
        .def("GetRelativeGain", &GainEstimator::GetRelativeGain)
        .def("GetRelativeGainsList", &GainEstimator::GetRelativeGainsList)
        ;  

// ============================ Spline123 =============================

    py::class_<Bspline1d>(m, "Spline1D")
        .def(py::init<double, double, int>())
        .def(py::init<std::string &>())
        .def("Eval", (std::vector <double> (Bspline1d::*)(std::vector <double> &) const) &Bspline1d::Eval)
        .def("Eval", (double (Bspline1d::*)(double) const) &Bspline1d::Eval)
        .def("EvalDrv", (std::vector <double> (Bspline1d::*)(std::vector <double> &) const) &Bspline1d::EvalDrv)
        .def("EvalDrv", (double (Bspline1d::*)(double) const) &Bspline1d::EvalDrv)
        .def("GetCoef", (std::vector<double> (Bspline1d::*)() const) &Bspline1d::GetCoef)
        .def("GetCoef", (double (Bspline1d::*)(int i) const) &Bspline1d::GetCoef)
        .def("SetCoef", (bool (Bspline1d::*)(std::vector <double> &)) &Bspline1d::SetCoef)
        .def("GetJsonString", (std::string (Bspline1d::*)() const) &Bspline1d::GetJsonString)
        .def("Basis", (double (Bspline1d::*)(double x, int n) const) &Bspline1d::Basis)
        .def("Basis", (std::vector <double> (Bspline1d::*) (std::vector <double> &vx, int n) const) &Bspline1d::Basis)
        .def("BasisDrv", (double (Bspline1d::*)(double x, int n) const) &Bspline1d::BasisDrv)
        .def("BasisDrv", (std::vector <double> (Bspline1d::*) (std::vector <double> &vx, int n) const) &Bspline1d::BasisDrv)
        ;

    py::class_<Bspline2d>(m, "Spline2D")
        .def(py::init<double, double, int, double, double, int>())
        .def(py::init<std::string &>())
        .def("Eval", (std::vector <double> (Bspline2d::*)(std::vector <double> &, std::vector <double> &) const) &Bspline2d::Eval)
        .def("Eval", (double (Bspline2d::*)(double, double) const) &Bspline2d::Eval)
        .def("EvalDrvX", (std::vector <double> (Bspline2d::*)(std::vector <double> &, std::vector <double> &) const) &Bspline2d::EvalDrvX)
        .def("EvalDrvX", (double (Bspline2d::*)(double, double) const) &Bspline2d::EvalDrvX)
        .def("EvalDrvY", (std::vector <double> (Bspline2d::*)(std::vector <double> &, std::vector <double> &) const) &Bspline2d::EvalDrvY)
        .def("EvalDrvY", (double (Bspline2d::*)(double, double) const) &Bspline2d::EvalDrvY)
        .def("GetCoef", (std::vector<double> (Bspline2d::*)() const) &Bspline2d::GetCoef)
        .def("GetCoef", (double (Bspline2d::*)(int i) const) &Bspline2d::GetCoef)
        .def("SetCoef", (bool (Bspline2d::*)(std::vector <double> &)) &Bspline2d::SetCoef)
        .def("GetJsonString", (std::string (Bspline2d::*)() const) &Bspline2d::GetJsonString)
        .def("Basis", (double (Bspline2d::*)(double, double, int, int) const) &Bspline2d::Basis)
        .def("Basis", (std::vector <double> (Bspline2d::*) (std::vector <double> &, std::vector <double> &, int, int) const) &Bspline2d::Basis)
        .def("BasisDrvX", (double (Bspline2d::*)(double, double, int, int) const) &Bspline2d::BasisDrvX)
        .def("BasisDrvX", (std::vector <double> (Bspline2d::*) (std::vector <double> &, std::vector <double> &, int, int) const) &Bspline2d::BasisDrvX)
        .def("BasisDrvY", (double (Bspline2d::*)(double, double, int, int) const) &Bspline2d::BasisDrvY)
        .def("BasisDrvY", (std::vector <double> (Bspline2d::*) (std::vector <double> &, std::vector <double> &, int, int) const) &Bspline2d::BasisDrvY)
        ;

    py::class_<Bspline3d>(m, "Spline3D")
        .def(py::init<double, double, int, double, double, int, double, double, int>())
        .def(py::init<std::string &>())
        .def("Eval", (std::vector <double> (Bspline3d::*)(std::vector <double> &, std::vector <double> &, std::vector <double> &) const) &Bspline3d::Eval)
        .def("Eval", (double (Bspline3d::*)(double, double, double) const) &Bspline3d::Eval)
        .def("GetJsonString", (std::string (Bspline3d::*)() const) &Bspline3d::GetJsonString)
        ;

    py::class_<BSfit1D>(m, "Fit1D")
        .def(py::init<double, double, int>())
//    	.def(py::init<Bspline1d*>())
        .def("AddData", (void (BSfit1D::*)(std::vector <double> &, std::vector <double> &)) &BSfit1D::AddData)
        .def("SetBinning", &BSfit1D::SetBinning)
        .def("SetBinningAuto", &BSfit1D::SetBinningAuto)
    	.def("Fit", (bool (BSfit1D::*)(std::vector <double> &, std::vector <double> &)) &BSfit1D::Fit)
        .def("BinnedFit", &BSfit1D::BinnedFit)
        .def("GetCoef", &BSfit1D::GetCoef)
        .def("MakeSpline", &BSfit1D::MakeSpline)
        .def("FitAndMakeSpline", &BSfit1D::FitAndMakeSpline)
    // read the underlying profile histogram 
        .def("GetBins", &BSfit1D::GetBins)
        .def("GetMeans", &BSfit1D::GetMeans)
        .def("GetSigmas", &BSfit1D::GetSigmas)
        .def("GetWeights", &BSfit1D::GetWeights)
        .def("GetXCenters", &BSfit1D::GetXCenters)
        .def("SetMinWeight", &BSfit1D::SetMinWeight)
        .def("SetMissingFactor", &BSfit1D::SetMissingFactor)
        ;      

    py::class_<ConstrainedFit1D>(m, "ConstrainedFit1D")
    // base class functions
        .def(py::init<double, double, int>())
//    	.def(py::init<Bspline1d*>())
        .def("AddData", (void (ConstrainedFit1D::*)(std::vector <double> &, std::vector <double> &)) &ConstrainedFit1D::AddData)
        .def("SetBinning", &ConstrainedFit1D::SetBinning)
        .def("SetBinningAuto", &ConstrainedFit1D::SetBinningAuto)
    	.def("Fit", (bool (ConstrainedFit1D::*)(std::vector <double> &, std::vector <double> &)) &ConstrainedFit1D::Fit)
        .def("BinnedFit", &ConstrainedFit1D::BinnedFit)
        .def("GetCoef", &ConstrainedFit1D::GetCoef)
        .def("MakeSpline", &ConstrainedFit1D::MakeSpline)
        .def("FitAndMakeSpline", &ConstrainedFit1D::FitAndMakeSpline)
    // read the underlying profile histogram
        .def("GetBins", &ConstrainedFit1D::GetBins)
        .def("GetMeans", &ConstrainedFit1D::GetMeans)
        .def("GetSigmas", &ConstrainedFit1D::GetSigmas)
        .def("GetWeights", &ConstrainedFit1D::GetWeights)
        .def("GetXCenters", &ConstrainedFit1D::GetXCenters)
        .def("SetMinWeight", &ConstrainedFit1D::SetMinWeight)
    // constraints
        .def("SetMinimum", &ConstrainedFit1D::SetMinimum)  
        .def("SetMaximum", &ConstrainedFit1D::SetMaximum) 
        .def("ForceNonNegative", &ConstrainedFit1D::ForceNonNegative)
        .def("ForceNonIncreasing", &ConstrainedFit1D::ForceNonIncreasing)
        .def("ForceNonDecreasing", &ConstrainedFit1D::ForceNonDecreasing)
        .def("FixAt", &ConstrainedFit1D::FixAt)
        .def("FixDrvAt", &ConstrainedFit1D::FixDrvAt)
        .def("FixLeft", &ConstrainedFit1D::FixLeft)
        .def("FixDrvLeft", &ConstrainedFit1D::FixDrvLeft)
        .def("FixRight", &ConstrainedFit1D::FixRight)
        .def("FixDrvRight", &ConstrainedFit1D::FixDrvRight)       
        ;

    py::class_<BSfit2D>(m, "Fit2D")
        .def(py::init<double, double, int, double, double, int>())
//    	.def(py::init<Bspline2d*>())
        .def("AddData", (void (BSfit2D::*)(std::vector <double> &, std::vector <double> &, std::vector <double> &)) &BSfit2D::AddData)
        .def("SetBinning", &BSfit2D::SetBinning)
        .def("SetBinningAuto", &BSfit2D::SetBinningAuto)
    	.def("Fit", (bool (BSfit2D::*)(std::vector <double> &, std::vector <double> &, std::vector <double> &)) &BSfit2D::Fit)
        .def("BinnedFit", &BSfit2D::BinnedFit)
        .def("GetCoef", &BSfit2D::GetCoef)
        .def("MakeSpline", &BSfit2D::MakeSpline)
        .def("FitAndMakeSpline", &BSfit2D::FitAndMakeSpline)
    // read the underlying profile histogram
        .def("GetBins", &BSfit2D::GetBins)
        .def("GetMeans", &BSfit2D::GetMeans)
        .def("GetSigmas", &BSfit2D::GetSigmas)
        .def("GetWeights", &BSfit2D::GetWeights)
        .def("GetXCenters", &BSfit2D::GetXCenters)
        .def("GetYCenters", &BSfit2D::GetYCenters)
        .def("SetMinWeight", &BSfit2D::SetMinWeight)
        .def("SetMissingFactor", &BSfit2D::SetMissingFactor)
        ;

    py::class_<ConstrainedFit2D>(m, "ConstrainedFit2D")
    // base class functions
        .def(py::init<double, double, int, double, double, int>())
//    	.def(py::init<Bspline2d*>())
        .def("AddData", (void (ConstrainedFit2D::*)(std::vector <double> &, std::vector <double> &, std::vector <double> &)) &ConstrainedFit2D::AddData)
        .def("SetBinning", &ConstrainedFit2D::SetBinning)
        .def("SetBinningAuto", &ConstrainedFit2D::SetBinningAuto)
    	.def("Fit", (bool (ConstrainedFit2D::*)(std::vector <double> &, std::vector <double> &, std::vector <double> &)) &ConstrainedFit2D::Fit)
        .def("BinnedFit", &ConstrainedFit2D::BinnedFit)
        .def("GetCoef", &ConstrainedFit2D::GetCoef)
        .def("MakeSpline", &ConstrainedFit2D::MakeSpline)
        .def("FitAndMakeSpline", &ConstrainedFit2D::FitAndMakeSpline)
    // read the underlying profile histogram
        .def("GetBins", &ConstrainedFit2D::GetBins)
        .def("GetMeans", &ConstrainedFit2D::GetMeans)
        .def("GetSigmas", &ConstrainedFit2D::GetSigmas)
        .def("GetWeights", &ConstrainedFit2D::GetWeights)
        .def("GetXCenters", &ConstrainedFit2D::GetXCenters)
        .def("GetYCenters", &ConstrainedFit2D::GetYCenters)
    // constraints
        .def("SetMinimum", &ConstrainedFit2D::SetMinimum)  
        .def("SetMaximum", &ConstrainedFit2D::SetMaximum)    
        .def("ForceNonNegative", &ConstrainedFit2D::ForceNonNegative)
        .def("ForceNonIncreasingX", &ConstrainedFit2D::ForceNonIncreasingX)
        .def("SetSlopeY", &ConstrainedFit2D::SetSlopeY)
        .def("ForceTopDown", &ConstrainedFit2D::ForceTopDown)
        .def("ForceFlatTopX", &ConstrainedFit2D::ForceFlatTopX)    
        ; 

    py::class_<BSfit3D>(m, "Fit3D")
        .def(py::init<double, double, int, double, double, int, double, double, int>())
        .def("AddData", (void (BSfit3D::*)(std::vector <double> &, std::vector <double> &, std::vector <double> &, std::vector <double> &)) &BSfit3D::AddData)
        .def("SetBinning", &BSfit3D::SetBinning)
        .def("SetBinningAuto", &BSfit3D::SetBinningAuto)
    	.def("Fit", (bool (BSfit3D::*)(std::vector <double> &, std::vector <double> &, std::vector <double> &, std::vector <double> &)) &BSfit3D::Fit)
        .def("BinnedFit", &BSfit3D::BinnedFit)
        .def("MakeSpline", &BSfit3D::MakeSpline)
        .def("FitAndMakeSpline", &BSfit3D::FitAndMakeSpline)
    // read the underlying profile histogram
        .def("GetBins", &BSfit3D::GetBins)
        .def("GetMeans", &BSfit3D::GetMeans)
        .def("GetSigmas", &BSfit3D::GetSigmas)
        .def("GetWeights", &BSfit3D::GetWeights)
        .def("GetXCenters", &BSfit3D::GetXCenters)
        .def("GetYCenters", &BSfit3D::GetYCenters)
        .def("GetZCenters", &BSfit3D::GetZCenters)
        .def("SetMinWeight", &BSfit3D::SetMinWeight)
        .def("SetMissingFactor", &BSfit3D::SetMissingFactor)
        ;

    py::class_<ConstrainedFit3D>(m, "ConstrainedFit3D")
        .def(py::init<double, double, int, double, double, int, double, double, int>())
        .def("AddData", (void (ConstrainedFit3D::*)(std::vector <double> &, std::vector <double> &, std::vector <double> &, std::vector <double> &)) &ConstrainedFit3D::AddData)
        .def("SetBinning", &ConstrainedFit3D::SetBinning)
        .def("SetBinningAuto", &ConstrainedFit3D::SetBinningAuto)
    	.def("Fit", (bool (ConstrainedFit3D::*)(std::vector <double> &, std::vector <double> &, std::vector <double> &, std::vector <double> &)) &ConstrainedFit3D::Fit)
        .def("BinnedFit", &ConstrainedFit3D::BinnedFit)
        .def("MakeSpline", &ConstrainedFit3D::MakeSpline)
        .def("FitAndMakeSpline", &ConstrainedFit3D::FitAndMakeSpline)
    // read the underlying profile histogram
        .def("GetBins", &ConstrainedFit3D::GetBins)
        .def("GetMeans", &ConstrainedFit3D::GetMeans)
        .def("GetSigmas", &ConstrainedFit3D::GetSigmas)
        .def("GetWeights", &ConstrainedFit3D::GetWeights)
        .def("GetXCenters", &ConstrainedFit3D::GetXCenters)
        .def("GetYCenters", &ConstrainedFit3D::GetYCenters)
        .def("GetZCenters", &ConstrainedFit3D::GetZCenters)
        .def("SetMinWeight", &ConstrainedFit3D::SetMinWeight)
        .def("SetMissingFactor", &ConstrainedFit3D::SetMissingFactor)
    // constraints
        .def("SetMinimum", &ConstrainedFit3D::SetMinimum)  
        .def("SetMaximum", &ConstrainedFit3D::SetMaximum)    
        .def("ForceNonNegative", &ConstrainedFit3D::ForceNonNegative)
        ;        

// ============================ ProfileHist =============================  
    py::class_<ProfileHist1D>(m, "ProfileHist1D")
        .def(py::init<int, double, double>())
        .def("Clear", &ProfileHist1D::Clear)
        .def("Fill", (bool (ProfileHist1D::*)(std::vector <double> &, std::vector <double> &)) &ProfileHist1D::Fill)
        .def("GetBinEntries", &ProfileHist1D::GetBinEntries)
        .def("GetBinMean", &ProfileHist1D::GetBinMean)
        .def("GetBinSigma", &ProfileHist1D::GetBinSigma)
        .def("GetMeans", &ProfileHist1D::GetMeans)
        .def("GetSigmas", &ProfileHist1D::GetSigmas)
        .def("GetWeights", &ProfileHist1D::GetWeights)
        .def("GetXCenters", &ProfileHist1D::GetXCenters)
        ;

    py::class_<ProfileHist2D>(m, "ProfileHist2D")
        .def(py::init<int, double, double, int, double, double>())
        .def("Clear", &ProfileHist2D::Clear)
        .def("Fill", (bool (ProfileHist2D::*)(std::vector <double> &, std::vector <double> &, std::vector <double> &)) &ProfileHist2D::Fill)
        .def("GetBinEntries", &ProfileHist2D::GetBinEntries)
        .def("GetBinMean", &ProfileHist2D::GetBinMean)
        .def("GetBinSigma", &ProfileHist2D::GetBinSigma)
        .def("GetMeans", &ProfileHist2D::GetMeans)
        .def("GetSigmas", &ProfileHist2D::GetSigmas)
        .def("GetWeights", &ProfileHist2D::GetWeights)
        .def("GetXCenters", &ProfileHist2D::GetXCenters)
        .def("GetYCenters", &ProfileHist2D::GetYCenters)        
        ; 

    py::class_<ProfileHist3D>(m, "ProfileHist3D")
        .def(py::init<int, double, double, int, double, double, int, double, double>())
        .def("Clear", &ProfileHist3D::Clear)
        .def("Fill", (bool (ProfileHist3D::*)(std::vector <double> &, std::vector <double> &, std::vector <double> &, std::vector <double> &)) &ProfileHist3D::Fill)
        .def("GetBinEntries", &ProfileHist3D::GetBinEntries)
        .def("GetBinMean", &ProfileHist3D::GetBinMean)
        .def("GetBinSigma", &ProfileHist3D::GetBinSigma)
        .def("GetMeans", &ProfileHist3D::GetMeans)
        .def("GetSigmas", &ProfileHist3D::GetSigmas)
        .def("GetWeights", &ProfileHist3D::GetWeights)
        .def("GetXCenters", &ProfileHist3D::GetXCenters)
        .def("GetYCenters", &ProfileHist3D::GetYCenters)        
        .def("GetZCenters", &ProfileHist3D::GetZCenters) 
        ;             
}

