#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include "lrf.h"
#include "lrmodel.h"
#include "lrfaxial.h"
//#include "lrfaxial3d.h"
//#include "lrfcomp.h"
// #include "reconstructor.h"
// #include "reconstructor_mp.h"
#include <cmath>

int main()
{  
    std::string jsonlrf("{\"type\": \"Axial\", \"rmin\": 0, \"rmax\": 7.605551275463986, \"nint\": 15, \"x0\": 0, \"y0\": 0, \"compression\": {\"method\": \"dualslope\", \"k\": 5, \"lam\": 0.5, \"r0\": 1}, \"constraints\": [\"non-negative\", \"non-increasing\", \"flattop\"]}");
    LRF* lrf = LRF::mkFromJson(jsonlrf);
    std::string json2 = lrf->GetJsonString();
    std::cout << json2 << std::endl << std::endl;

    LRF* lrf2 = LRF::mkFromJson(json2);
    std::cout << lrf->GetJsonString() << std::endl << std::endl;

    std::cout << "valid: " << lrf2->isValid() << ", ready: " << lrf2->isReady() << std::endl;

    LRF* lrfc = lrf2->clone();
    std::cout << "valid: " << lrfc->isValid() << ", ready: " << lrfc->isReady() << std::endl;
    std::cout << lrf2->eval(0,1,0) << std::endl;
    std::cout << lrfc->eval(0,1,0) << std::endl;
/*
    LRF *lrf1 = ((LRFcomp*)lrf)->GetLayer(0)->clone();
    LRF *lrf2 = ((LRFcomp*)lrf)->GetLayer(1)->clone();
    LRF *lrf3 = LRF::mkFromJson(lrf2->GetJsonString());

    std::cout << ((LRFaxial3d*)lrf2)->isReady() << std::endl;
    std::cout << ((LRFaxial3d*)lrf3)->clone()->isReady() << std::endl;
    std::cout << lrf2->eval(0, 0, 900) << std::endl;
    std::cout << lrf3->clone()->eval(0, 0, 900) << std::endl;

    std::cout << lrf->eval(0, 0, 900) << std::endl;
*/    
    return 0;
}

