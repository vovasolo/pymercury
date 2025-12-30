#ifndef LRFORMULAV_H
#define LRFORMULAV_H

#include "lrf.h"
#include "profileHist.h"
#include <cmath>

#include <Eigen/Core>
#include <cstddef>
#include <string>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>
#include <vector>
#include <map>

#include "vformula.h"

// Define a functor for LM optimization
// this particualar structure is needed by LM optimizer
// plus some additional stuff required by numerical differentiation module
template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct FunctorV {
    typedef _Scalar Scalar;
    enum {
        InputsAtCompileTime = NX,
        ValuesAtCompileTime = NY
    };
    typedef Eigen::Matrix<Scalar, NX, 1> InputType;
    typedef Eigen::Matrix<Scalar, NY, 1> ValueType;
    typedef Eigen::Matrix<Scalar, NY, NX> JacobianType;

    int m_inputs, m_values;
    FunctorV(int inputs, int values) : m_inputs(inputs), m_values(values) {}
    int inputs() const { return m_inputs; }
    int values() const { return m_values; }
};

// Universal functor
struct UniFunctorV : FunctorV<double> {
    Eigen::VectorXd x;
    Eigen::VectorXd y;
    VFormula *vf;

    UniFunctorV(const Eigen::VectorXd& x_, const Eigen::VectorXd& y_, VFormula *vf_)
        : FunctorV<double>(vf_->GetConstCount()-1, x_.size()), x(x_), y(y_), vf(vf_){}

    // Compute residuals: f(p) = model(p) - y
    int operator()(const Eigen::VectorXd& p, Eigen::VectorXd& fvec) const {
        for (size_t i=1; i<vf->GetConstCount(); i++)
            vf->SetConstant(i, p[i-1]);

        fvec = vf->Eval(x.array()) - y.array();
        return 0;
    }
};

// Universal functor with weights
struct UniFunctorWV : FunctorV<double> {
    Eigen::VectorXd x;
    Eigen::VectorXd y;
    Eigen::VectorXd w;
    VFormula *vf;

    UniFunctorWV(const Eigen::VectorXd& x_, const Eigen::VectorXd& y_, const Eigen::VectorXd& w_, VFormula *vf_)
        : FunctorV<double>(vf_->GetConstCount()-1, x_.size()), x(x_), y(y_), w(w_), vf(vf_){}

    // Compute weighted residuals: f(p) = (model(p) - y) * sqrt(w)
    int operator()(const Eigen::VectorXd& p, Eigen::VectorXd& fvec) const {
        for (size_t i=1; i<vf->GetConstCount(); i++)
            vf->SetConstant(i, p[i-1]);

        fvec = (vf->Eval(x.array()) - y.array()) * sqrt(w.array());
        return 0;
    }
};

class LRFormulaV : public LRF
{
public:
 LRFormulaV(double x0, double y0, double rmax);
 LRFormulaV(const Json &json);
 LRFormulaV(std::string &json_str);    
//    LRFormula1(const LRFormula1 &obj); // copy constructor
    LRFormulaV();

    virtual LRFormulaV* clone() const;

    virtual bool inDomain(double x, double y, double z=0.) const;
    virtual bool isReady () const;
    virtual double getRmax() const { return rmax; }
    virtual double eval(double x, double y, double z=0.) const;
    double evalAxial(double r) const;
//    double evalDrv(double r) const;
    virtual double evalDrvX(double x, double y, double z=0.) const {return 0.;}
    virtual double evalDrvY(double x, double y, double z=0.) const {return 0.;}

    virtual bool fitData(const std::vector <LRFdata> &data);
    virtual bool addData(const std::vector <LRFdata> &data);
    virtual bool doFit();
    virtual void clearData() { if (h1) h1->Clear(); }

    virtual std::string type() const { return std::string("FormulaV"); }
    virtual void ToJsonObject(Json_object &json) const;

    void SetOrigin(double x, double y) {x0 = x; y0 = y; Init();}
    void SetRmin(double r) { rmin = std::max(r, 0.); Init();}
    void SetRmax(double r) { rmax = r; Init();}

//  VFormula-related calls
    // setters
    void SetParameter(std::string name, double val);
    void SetExpression(std::string expr) {expression = expr;}
    // call this to update the virtual machine if parameter names of expression were changed
    std::string InitVF();  

    // getters
    double GetParameter(std::string name);
    std::vector<double> GetParVector() const {return parvals;}

// VM-debugging
    std::vector<std::string> GetPrg() {return vf->GetPrg();}
    std::vector<std::string> GetConstMap() {return vf->GetConstMap();}
    std::vector<std::string> GetVarMap() {return vf->GetVarMap();}
    std::vector<std::string> GetOperMap() {return vf->GetOperMap();}
    std::vector<std::string> GetFuncMap() {return vf->GetFuncMap();} 
    std::string GetVersion() {return std::string("1.0");}   

// calculation of radius + provision for compression
    double R(double x, double y) const {return sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0));}
    double R2(double x, double y) const {return (x-x0)*(x-x0)+(y-y0)*(y-y0);}

// public getters
    double GetRmin() const {return rmin;}
    double GetOriginX() const {return x0;}
    double GetOriginY() const {return y0;}
    ProfileHist *GetHist() {return h1;}
    

// relative gain calculation
    double GetRatio(LRF* other) const; 

protected:
    void Init();

protected:
    double x0 = 0., y0 = 0.;  // center
    double rmin = 0.;    // domain
    double rmax = 0.;	// domain
    double rmin2;   // domain
    double rmax2;	// domain
    bool init_done = false;

// prof. histogram used for binned fitting
    ProfileHist1D *h1 = nullptr; 
    int nbins = 20;

// VFormula
    VFormula *vf = nullptr;
    std::string expression = std::string("");   // expression to parse
// vectors because parameters passed in a vector to functors
    std::vector<std::string> parnames;             // parameter names
    std::vector<double> parvals;                   // parameter values
    std::map<std::string, size_t> parmap;          // map for bookkeeping

// safeguards
    double rzerod = 1e-6; // derivatives are assumed to be zero below this radius 
};

#endif // LRFORMULAV_H
