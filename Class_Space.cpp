////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author:          Bowen Zhao <bwzhao@bu.edu>
// Organization:    Boston University
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "Class_Space.h"


void SSE::Class_Space::Set_Class_Space(const SSE::Class_Lattice &_which_Lattice) {
    this->Array_Spin.resize(_which_Lattice.Get_NumSite());

    // Initialize the spin array
    static std::random_device ran_dev;
    static std::uniform_real_distribution<SSE::type_DataFloat> Ran_Int(0, 1);

    for (auto &spin : this->Array_Spin) {
        if (Ran_Int(ran_dev) >= 0.5){
            spin = 1;
        }
        else{
            spin = -1;
        }
    }
}
