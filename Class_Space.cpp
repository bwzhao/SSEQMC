////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author:          Bowen Zhao <bwzhao@bu.edu>
// Organization:    Boston University
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "Class_Space.h"


void SSE::Class_Space::Set_Class_Space(const SSE::Class_Lattice &_which_Lattice) {
    this->Array_Spin.resize(_which_Lattice.Get_NumSite());

    // Initialize the spin array
//    static std::random_device ran_dev;
//    static std::uniform_real_distribution<SSE::type_DataFloat> Ran_Int(0, 1);
//
//    for (auto &spin : this->Array_Spin) {
//        if (Ran_Int(ran_dev) >= 0.5){
//            spin = 1;
//        }
//        else{
//            spin = -1;
//        }
//    }

    std::vector<type_DataInt> SitaA_Array;
    std::vector<type_DataInt> SitaB_Array;
    for (type_DataInt index_x = 0; index_x != _which_Lattice.Get_ParaHamil().Num_X; ++index_x) {
        for (type_DataInt index_y = 0; index_y != _which_Lattice.Get_ParaHamil().Num_Y; ++index_y) {
            type_DataInt nl_site = index_x + index_y * _which_Lattice.Get_ParaHamil().Num_X;
            if ((index_x + index_y) % 2 == 0) {
                SitaA_Array.push_back(nl_site);
            }
            else {
                SitaB_Array.push_back(nl_site);
            }
        }
    }

    static std::random_device ran_dev;
    static std::uniform_real_distribution<SSE::type_DataFloat> Ran_Int(0, 1);

    std::shuffle(SitaA_Array.begin(), SitaA_Array.end(), ran_dev);
    std::shuffle(SitaB_Array.begin(), SitaB_Array.end(), ran_dev);

    for(type_DataInt index_site = 0; index_site != SitaA_Array.size(); ++index_site){
        auto siteA = SitaA_Array[index_site];
        auto siteB = SitaB_Array[index_site];

        //Site
        if (Ran_Int(ran_dev) >= 0.5){
            Array_Spin[siteA] = -1;
            Array_Spin[siteB] = 1;
        }
        else{
            Array_Spin[siteA] = 1;
            Array_Spin[siteB] = -1;
        }
    }
}
