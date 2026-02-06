#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include "lrf.h"
#include "lrmodel.h"
#include "lrfaxial.h"
#include "lrfxy.h"
//#include "lrfaxial3d.h"
//#include "lrfcomp.h"
// #include "reconstructor.h"
// #include "reconstructor_mp.h"
#include <cmath>

int main()
{  
    std::string json1("{\"constraints\": [\"non-negative\"], \"type\": \"XY\", \"xmax\": 4, \"xmin\": -4, \"nintx\": 15, \"ymax\": 4, \"ymin\": -4, \"ninty\": 15}");
    LRF* lrf = LRF::mkFromJson(json1);
    if (!lrf) {
        std::cout << "Failed to create LRF\n";
        return -1;
    }
    std::cout << "valid: " << lrf->isValid() << ", ready: " << lrf->isReady() << std::endl;

    std::string json2 = lrf->GetJsonString();
    std::cout << json2 << std::endl << std::endl;      

    LRF* lrf2 = LRF::mkFromJson(json2);
    std::cout << lrf2->GetJsonString() << std::endl << std::endl;
    std::cout << "valid: " << lrf2->isValid() << ", ready: " << lrf2->isReady() << std::endl;

    LRF* lrfc = lrf2->clone();
    std::cout << "valid: " << lrfc->isValid() << ", ready: " << lrfc->isReady() << std::endl;
    std::cout << "Normal: " << lrf2->eval(0,1,0) << std::endl;
    std::cout << "Cloned: " << lrfc->eval(0,1,0) << std::endl;
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

