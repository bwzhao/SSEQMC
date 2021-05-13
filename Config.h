////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author:          Bowen Zhao <bwzhao@bu.edu>
// Organization:    Boston University
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once
#include <cstdint>
#include <vector>
#include <string>
#include <ctime>

namespace SSE{
    // Data type
    using type_DataFloat = double;
    using type_DataInt = std::int_fast32_t;

    using type_NumSite = std::int_fast32_t;

    // Data structure type
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Hamiltonian Parameter type
    struct type_ParaHamil{
        type_DataFloat ty;
        type_DataFloat tx;
        type_DataFloat Jy;
        type_DataFloat Jx;
        type_DataFloat Q2;
        type_DataFloat Q2s;
        type_DataFloat JfA;
        type_DataFloat JfB;
        type_DataFloat Q3;
        type_DataFloat Q3s;

        type_NumSite Num_X;
        type_NumSite Num_Y;
        type_DataFloat beta;
    };

    constexpr type_DataInt MODELTYPE_tJ_PBC = 0;



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Model specific type
    using type_SiteinOper = std::vector<SSE::type_NumSite>;

    // Type of vertex (opeartor with legs)
    constexpr type_DataInt NUM_TYPE_OPER = 3;
    constexpr type_DataInt OPERTYPE_ID = -1;
    constexpr type_DataInt OPERTYPE_t = 0;
    constexpr type_DataInt OPERTYPE_J = 1;
    constexpr type_DataInt OPERTYPE_Q2 = 2;


    // DetailedType of vertex (opeartor with legs)
    constexpr type_DataInt NUM_DETAILEDTYPE = 6;
    constexpr type_DataInt DETAILEDTYPE_ID = -1;
    constexpr type_DataInt DETAILEDTYPE_ty = 0;
    constexpr type_DataInt DETAILEDTYPE_tx = 1;
    constexpr type_DataInt DETAILEDTYPE_Jy = 2;
    constexpr type_DataInt DETAILEDTYPE_Jx = 3;
    constexpr type_DataInt DETAILEDTYPE_Q2y = 4;
    constexpr type_DataInt DETAILEDTYPE_Q2x = 5;

    // Properties:


    // Number of Sites
    constexpr type_DataInt MAX_NUMSITE = 4;
    constexpr type_DataInt NUMSITE_t = 2;
    constexpr type_DataInt NUMSITE_J = 2;
    constexpr type_DataInt NUMSITE_Q2 = 4;
    const std::vector<type_DataInt> ARRAY_NUMSITE = {NUMSITE_t, NUMSITE_J, NUMSITE_Q2};


    constexpr type_DataInt NUM_LEG_PER_SITE = 2;

    // Check pair
//    const std::vector<std::vector<type_DataInt>> CHECKPAIR_t = {{0, 1}};
//    const std::vector<std::vector<type_DataInt>> CHECKPAIR_J = {{0, 1}};
//    const std::vector<std::vector<type_DataInt>> CHECKPAIR_Q2 = {{0, 1}, {2, 3}};
//    const std::vector<std::vector<std::vector<type_DataInt>>> Array_CHECKPAIR = {{CHECKPAIR_t,
//                                                                                     CHECKPAIR_J,
//                                                                                     CHECKPAIR_Q2}};

    // Update site:
//    const std::vector<std::vector<type_DataInt>> UPDATESITE_t = {};
//    const std::vector<std::vector<type_DataInt>> UPDATESITE_J = {{}, {0}};
//    const std::vector<std::vector<type_DataInt>> UPDATESITE_Q2 = {{}, {0}, {1}, {0, 1}};
//    const  std::vector<std::vector<std::vector<type_DataInt>>> Array_UPDATESITE = {{UPDATESITE_t,
//                                                                                           UPDATESITE_J,
//                                                                                           UPDATESITE_Q2}};

    const std::vector<type_DataInt> T_t = {{1, 0, 3, 2}};

    const std::vector<type_DataInt> T_J = {{2, 3, 0, 1}};

    const std::vector<type_DataInt> T_Q2 = {{2, 3, 0, 1, 2 + 4, 3 + 4, 0 + 4, 1 + 4}};


    // Measurement
    const std::vector<std::string> Vec_Name{"E",
                                            "n",
                                            "n^2",
                                            "m_z^2",
                                            "m_z^4",
                                            "m_c^2",
                                            "m_c^4",
                                            "Phi_4_c",
                                            "m_s^2",
                                            "m_s^4",
                                            "rho_s",
                                            "chi_u"
    };

    const int INDEX_E = 0;
    const int INDEX_N = 1;
    const int INDEX_N2 = 2;

    const int start_M = 3;
    const int INDEX_MZ2 = start_M + 0;
    const int INDEX_MZ4 = start_M + 1;
    const int INDEX_MC2 = start_M + 2;
    const int INDEX_MC4 = start_M + 3;
    const int INDEX_PHI4C = start_M + 4;
    const int INDEX_MS2 = start_M + 5;
    const int INDEX_MS4 = start_M + 6;

    const int start_rhos = 10;
    const int INDEX_RHOS = start_rhos + 0;
    const int INDEX_CHIU = start_rhos + 1;

    const std::vector<std::string> Vec_OperNum{"n_t",
                                               "n_J",
                                               "n_Q2",
                                        };

    const int INDEX_Nt = 0;
    const int INDEX_NJ = 1;
    const int INDEX_NQ2 = 2;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Other parameters
    constexpr type_DataInt MIN_INITIAL_CUTOFF = 20;
    // Index of left and right space

    constexpr type_DataFloat PI = 3.14159265358979323846;
}
