////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author:          Bowen Zhao <bwzhao@bu.edu>
// Organization:    Boston University
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once
#include <vector>
#include "Config.h"
#include "Class_Lattice.h"
#include <random>
#include <array>
#include <numeric>

namespace SSE{
    class Class_Space {
    private:
        // Store the spin configuration
        std::vector<SSE::type_DataInt> Array_Spin;
    public:
        void Set_Class_Space(const SSE::Class_Lattice & _which_Lattice);

        SSE::type_DataInt Get_Spin(type_DataInt _which_site) const;
        void Set_Spin(type_DataInt _which_site, type_DataInt _which_value);
        void Flip_Spin(type_DataInt _which_site);

        bool Is_Same(const SSE::Class_Space &_other_space) const;

        //Measurement
        SSE::type_DataFloat Measure_mz(const SSE::Class_Lattice & _which_Lattice) const;
        SSE::type_DataFloat Measure_mu(const SSE::Class_Lattice & _which_Lattice) const;
        std::array<SSE::type_DataFloat, 4> Measure_md(const SSE::Class_Lattice & _which_Lattice) const;

    };
}


inline SSE::type_DataInt SSE::Class_Space::Get_Spin(type_DataInt _which_site) const{
    return Array_Spin[_which_site];
}

inline void SSE::Class_Space::Flip_Spin(SSE::type_DataInt _which_site) {
    Array_Spin[_which_site] *= -1;
}

inline bool SSE::Class_Space::Is_Same(const SSE::Class_Space &_other_space) const {
    return Array_Spin == _other_space.Array_Spin;
}

inline void SSE::Class_Space::Set_Spin(SSE::type_DataInt _which_site, SSE::type_DataInt _which_value) {
    Array_Spin[_which_site] = _which_value;
}

inline SSE::type_DataFloat SSE::Class_Space::Measure_mz(const SSE::Class_Lattice &_which_Lattice) const {

    SSE::type_DataFloat sum_mz = std::inner_product(Array_Spin.cbegin(), Array_Spin.cend(), _which_Lattice.Get_ArraySublattice().cbegin(), 0);
    sum_mz /= _which_Lattice.Get_NumSite() * 2;

    return sum_mz;
}

inline SSE::type_DataFloat SSE::Class_Space::Measure_mu(const SSE::Class_Lattice &_which_Lattice) const {

    SSE::type_DataFloat sum_mu = std::accumulate(Array_Spin.cbegin(), Array_Spin.cend(),0);
    sum_mu /= _which_Lattice.Get_NumSite() * 0.5;

    return sum_mu;
}

inline std::array<SSE::type_DataFloat, 4> SSE::Class_Space::Measure_md(const SSE::Class_Lattice &_which_Lattice) const {
    SSE::type_DataInt sum_mdc_x = 0;
    SSE::type_DataInt sum_mdc_y = 0;
    SSE::type_DataInt sum_mds_x = 0;
    SSE::type_DataInt sum_mds_y = 0;

    for (SSE::type_NumSite index_site = 0; index_site != _which_Lattice.Get_NumSite(); ++index_site) {
        const auto &index_site_Tx = _which_Lattice.Get_Map_Site_x1y0()[index_site];
        const auto &index_site_Ty = _which_Lattice.Get_Map_Site_x0y1()[index_site];

        //bond0 y
        SSE::type_DataInt temp_Dy = Array_Spin[index_site] * Array_Spin[index_site_Ty];
        SSE::type_DataInt temp_Dx = Array_Spin[index_site] * Array_Spin[index_site_Tx];

        sum_mdc_x += temp_Dx * _which_Lattice.Get_Qx()[index_site];
        sum_mdc_y += temp_Dy * _which_Lattice.Get_Qy()[index_site];
        sum_mds_x += temp_Dx * _which_Lattice.Get_Qs()[index_site];
        sum_mds_y += temp_Dy * _which_Lattice.Get_Qs()[index_site];
    }

    return std::array<SSE::type_DataFloat, 4>{{static_cast<type_DataFloat>(sum_mdc_x) / 4./ _which_Lattice.Get_NumSite() ,
                                                      static_cast<type_DataFloat>(sum_mdc_y)  / 4. / _which_Lattice.Get_NumSite(),
                                                      static_cast<type_DataFloat>(sum_mds_x) / 4. / _which_Lattice.Get_NumSite(),
                                                      static_cast<type_DataFloat>(sum_mds_y) / 4. / _which_Lattice.Get_NumSite()}};
}
