////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author:          Bowen Zhao <bwzhao@bu.edu>
// Organization:    Boston University
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "Config.h"
#include <random>
#include <ctime>

namespace SSE{
    // One engine shared by every update, so the updates do not draw from identically seeded copies.
    inline std::mt19937_64 &Get_Engine(){
        static std::mt19937_64 engine;
        return engine;
    }

    // Mix hardware entropy with the time and the CPU index, so jobs started in the same (or adjacent)
    // second still get different seeds. Must be called before any random number is drawn.
    inline void Seed_Engine(type_DataInt _which_cpu){
        std::random_device ran_dev;
        std::seed_seq seed{ran_dev(), ran_dev(), ran_dev(), ran_dev(),
                           static_cast<unsigned>(std::time(nullptr)),
                           static_cast<unsigned>(_which_cpu)};
        Get_Engine().seed(seed);
    }
}
