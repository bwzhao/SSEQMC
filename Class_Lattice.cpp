////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author:          Bowen Zhao <bwzhao@bu.edu>
// Organization:    Boston University
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "Class_Lattice.h"
SSE::Class_Lattice::Class_Lattice(type_ParaHamil _ParaHamil,
                                  type_DataInt _which_cpu,
                                  type_DataInt _bc,
                                  type_DataInt _Num_Segment):
        ParaHamil(_ParaHamil),
        Which_BC(_bc),
        Num_Site(_ParaHamil.Num_X * _ParaHamil.Num_Y),
        Which_CPU(_which_cpu),
        Ratio_Sum(0.)
{
    //Initalize the Sublattice
    for (type_NumSite index_site = 0; index_site != Num_Site; ++index_site) {
        auto index_x = index_site % ParaHamil.Num_X;
        auto index_y = index_site / ParaHamil.Num_X;
        auto index_sitephase = index_x + index_y;
        if (index_sitephase % 2 == 0) {
            Map_Site_Sublattice.push_back(1);
        } else {
            Map_Site_Sublattice.push_back(-1);
        }
    }

    //Initalize the plaquettePhase
    const auto Label_Qs = SSE::Class_SiteLabel(PI, PI);
    const auto Label_Qx = SSE::Class_SiteLabel(PI, 0.);
    const auto Label_Qy = SSE::Class_SiteLabel(0., PI);
    for (type_NumSite index_y = 0; index_y != _ParaHamil.Num_Y; ++index_y) {
        for (type_NumSite index_x = 0; index_x != _ParaHamil.Num_X; ++index_x) {
            SSE::Class_SiteLabel temp_SiteLabel(index_x, index_y);
            if (index_x + index_y == 0) {
                Map_Site_Qs.push_back(1);
            }
            else{
                Map_Site_Qs.push_back(-1);
            }
            if (index_x == 0) {
                Map_Site_Qx.push_back(1);
            }
            else{
                Map_Site_Qx.push_back(-1);
            }

            if (index_y == 0) {
                Map_Site_Qy.push_back(1);
            }
            else{
                Map_Site_Qy.push_back(-1);
            }
        }
    }

    // Map site-diffsite
    for (type_NumSite index_site = 0; index_site != Num_Site; ++index_site) {
        auto index_x = index_site % ParaHamil.Num_X;
        auto index_y = index_site / ParaHamil.Num_X;
        std::vector<SSE::type_DataInt> temp_vec;
        for (type_NumSite index_diff = 0; index_diff != Num_Site; ++index_diff) {
            auto index_diff_x = index_diff % ParaHamil.Num_X;
            auto index_diff_y = index_diff / ParaHamil.Num_X;

            auto new_index_x = (index_x + index_diff_x + ParaHamil.Num_X) % ParaHamil.Num_X;
            auto new_index_y = (index_y + index_diff_y + ParaHamil.Num_Y) % ParaHamil.Num_Y;

            temp_vec.emplace_back(new_index_x + new_index_y * ParaHamil.Num_X);
        }
        Map_Site_DiffSite.emplace_back(temp_vec);
    }

    //Initialize Oper_Site
    for (SSE::type_NumSite index_y = 0; index_y != _ParaHamil.Num_Y; ++index_y) {
        for (SSE::type_NumSite index_x = 0; index_x != _ParaHamil.Num_X; ++index_x) {
            auto nl_unit_x0y0 = index_x + index_y * _ParaHamil.Num_X;

            auto nl_unit_x1y0 = ((index_x + 1) % _ParaHamil.Num_X) + index_y * _ParaHamil.Num_X;
            auto nl_unit_x0y1 = (index_x + ((index_y + 1) % _ParaHamil.Num_Y) * _ParaHamil.Num_X);

            auto nl_unit_x2y0 = ((index_x + 2) % _ParaHamil.Num_X) + index_y * _ParaHamil.Num_X;
            auto nl_unit_x0y2 = (index_x + ((index_y + 2) % _ParaHamil.Num_Y) * _ParaHamil.Num_X);

            auto nl_unit_x1y1 = ((index_x + 1) % _ParaHamil.Num_X) + ((index_y + 1) % _ParaHamil.Num_Y) * _ParaHamil.Num_X;
            auto nl_unit_x1y2 = ((index_x + 1) % _ParaHamil.Num_X) + ((index_y + 2) % _ParaHamil.Num_Y) * _ParaHamil.Num_X;
            auto nl_unit_x2y1 = ((index_x + 2) % _ParaHamil.Num_X) + ((index_y + 1) % _ParaHamil.Num_Y) * _ParaHamil.Num_X;
            auto nl_unit_x2y2 = ((index_x + 2) % _ParaHamil.Num_X) + ((index_y + 2) % _ParaHamil.Num_Y) * _ParaHamil.Num_X;
            auto nl_unit_x3y2 = ((index_x + 3) % _ParaHamil.Num_X) + ((index_y + 2) % _ParaHamil.Num_Y) * _ParaHamil.Num_X;
            auto nl_unit_x2y3 = ((index_x + 2) % _ParaHamil.Num_X) + ((index_y + 3) % _ParaHamil.Num_Y) * _ParaHamil.Num_X;

            auto nl_unit_xm1y1 = ((index_x - 1 + _ParaHamil.Num_X) % _ParaHamil.Num_X) + ((index_y + 1) % _ParaHamil.Num_Y) * _ParaHamil.Num_X;
            auto nl_unit_xm1y2 = ((index_x - 1 + _ParaHamil.Num_X) % _ParaHamil.Num_X) + ((index_y + 2) % _ParaHamil.Num_Y) * _ParaHamil.Num_X;
            auto nl_unit_xm2y2 = ((index_x - 2 + _ParaHamil.Num_X) % _ParaHamil.Num_X) + ((index_y + 2) % _ParaHamil.Num_Y) * _ParaHamil.Num_X;

            auto nl_unit_x1ym1 = ((index_x + 1) % _ParaHamil.Num_X) + ((index_y - 1 + _ParaHamil.Num_Y) % _ParaHamil.Num_Y) * _ParaHamil.Num_X;
            auto nl_unit_x2ym1 = ((index_x + 2) % _ParaHamil.Num_X) + ((index_y - 1 + _ParaHamil.Num_Y) % _ParaHamil.Num_Y) * _ParaHamil.Num_X;
            auto nl_unit_x2ym2 = ((index_x + 2) % _ParaHamil.Num_X) + ((index_y - 2 + _ParaHamil.Num_Y) % _ParaHamil.Num_Y) * _ParaHamil.Num_X;

            Map_Site_x1y0.emplace_back(nl_unit_x1y0);
            Map_Site_x0y1.emplace_back(nl_unit_x0y1);

            if (Which_BC == MODELTYPE_tJ_PBC) {
                // t-terms
                // 0: t term: site(x, y)-site(x, y+1)
                Map_Oper_Site.emplace_back(SSE::type_SiteinOper{{nl_unit_x0y0, nl_unit_x0y1}});

                // 1: t term: site(x, y)-site(x+1, y)
                Map_Oper_Site.emplace_back(SSE::type_SiteinOper{{nl_unit_x0y0, nl_unit_x1y0}});

                // J-terms
                // 2: J term: site(x, y)-site(x, y+1)
                Map_Oper_Site.emplace_back(SSE::type_SiteinOper{{nl_unit_x0y0, nl_unit_x0y1}});
                // 3: J term: site(x, y)-site(x+1, y)
                Map_Oper_Site.emplace_back(SSE::type_SiteinOper{{nl_unit_x0y0, nl_unit_x1y0}});

                // Q2-terms
                // 4: Q term: site(x, y)-site(x, y+1)-site(x+1, y)-site(x+1,y+1)
                Map_Oper_Site.emplace_back(SSE::type_SiteinOper{{nl_unit_x0y0,
                                                                 nl_unit_x0y1,
                                                                 nl_unit_x1y0,
                                                                 nl_unit_x1y1}});
                // 5: Q term: site(x, y)-site(x+1, y)-site(x+1, y)-site(x+1,y+1)
                Map_Oper_Site.emplace_back(SSE::type_SiteinOper{{nl_unit_x0y0,
                                                                 nl_unit_x1y0,
                                                                 nl_unit_x0y1,
                                                                 nl_unit_x1y1}});
            }
        }
    }

    //Initialize the OperWeight/ Opertype / Mat
    if (Which_BC == MODELTYPE_tJ_PBC) {
        for (type_NumSite index_unit = 0; index_unit != Num_Site; ++index_unit) {
            Map_Oper_OperWeight.emplace_back(ParaHamil.ty);
            Map_Oper_OperWeight.emplace_back(ParaHamil.tx);
            Map_Oper_OperWeight.emplace_back((ParaHamil.Jy) / 2.);
            Map_Oper_OperWeight.emplace_back((ParaHamil.Jx) / 2.);
            Map_Oper_OperWeight.emplace_back((ParaHamil.Q2) / 4.);
            Map_Oper_OperWeight.emplace_back((ParaHamil.Q2) / 4.);

            if (ParaHamil.Num_Y == 1){
                if (ParaHamil.ty != 0) {
                    throw std::runtime_error("ty not consistant with 1D system");
                }
                if (ParaHamil.Jy != 0) {
                    throw std::runtime_error("Jy not consistant with 1D system");
                }
                if (ParaHamil.Q2 != 0) {
                    throw std::runtime_error("Q3 not consistant with 1D system");
                }
            }

            Map_Oper_Opertype.emplace_back(OPERTYPE_t);
            Map_Oper_Opertype.emplace_back(OPERTYPE_t);
            Map_Oper_Detailedtype.emplace_back(DETAILEDTYPE_ty);
            Map_Oper_Detailedtype.emplace_back(DETAILEDTYPE_tx);

            Map_Oper_Opertype.emplace_back(OPERTYPE_J);
            Map_Oper_Opertype.emplace_back(OPERTYPE_J);
            Map_Oper_Detailedtype.emplace_back(DETAILEDTYPE_Jy);
            Map_Oper_Detailedtype.emplace_back(DETAILEDTYPE_Jx);

            Map_Oper_Opertype.emplace_back(OPERTYPE_Q2);
            Map_Oper_Opertype.emplace_back(OPERTYPE_Q2);
            Map_Oper_Detailedtype.emplace_back(DETAILEDTYPE_Q2y);
            Map_Oper_Detailedtype.emplace_back(DETAILEDTYPE_Q2x);
        }
    }
    Num_Bond = Map_Oper_Opertype.size();

    //Initialize Ratio_Sum
    Ratio_Sum = std::accumulate(Map_Oper_OperWeight.begin(), Map_Oper_OperWeight.end(), 0.) * (_ParaHamil.beta / _Num_Segment);

    auto norm = std::accumulate(Map_Oper_OperWeight.begin(), Map_Oper_OperWeight.end(), 0.);
    for (auto & ele: Map_Oper_OperWeight){
        ele /= norm;
    }
}


