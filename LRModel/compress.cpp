#include "compress.h"
#include "json11.hpp"

Compress1d* Compress1d::Factory(const Json &json)
{
    if (json["method"].is_string() && json["method"].string_value() == "dualslope") {
        Compress1d *comp = new DualSlopeCompress(json);
        if (comp->fValid)
            return comp;
        else {
            delete comp;
            return nullptr;
        }
    } else if (json["method"].is_string() && json["method"].string_value() == "quadratic") {
        Compress1d *comp = new QuadraticCompress();
        return comp;
    } else {
        return nullptr;
    }
}

void DualSlopeCompress::Init()
{
    if (r0 < 0. || k <= 1.)
        return;

    a = (k+1)/(k-1);
    lam2 = std::max(lam*lam, 1e-6); // hack to avoid division by zero when lam=0
    b = sqrt(r0*r0+lam2)+a*r0;

    fValid = true;
}

DualSlopeCompress::DualSlopeCompress(double k, double r0, double lam) :
    k(k), r0(r0), lam(lam)
{
    Init();
}

DualSlopeCompress::DualSlopeCompress(const Json &json)
{
    if (!json["r0"].is_number() || !json["k"].is_number() || !json["lam"].is_number())
        return;

    r0 = json["r0"].number_value();
    k = json["k"].number_value();
    lam = json["lam"].number_value();

    Init();
}

double DualSlopeCompress::Rho(double r) const
{
    double dr = r - r0;
    return std::max(0., b + dr*a - sqrt(dr*dr + lam2));
}

double DualSlopeCompress::RhoDrv(double r) const
{
    double dr = r - r0;
// Andrey has spotted that this can lead to division by zero if lam2=0
// Mitigated by adding a small positive lower limit to lam2 in Init() 
    return dr/sqrt(dr*dr + lam2) + a;
}

void  DualSlopeCompress::ToJsonObject(Json_object &json) const
{
    json["method"] = "dualslope";
    json["r0"] = r0;
    json["lam"] = lam;
    json["k"] = k;
}

void  QuadraticCompress::ToJsonObject(Json_object &json) const
{
    json["method"] = "quadratic";
}
