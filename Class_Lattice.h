////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author:          Bowen Zhao <bwzhao@bu.edu>
// Organization:    Boston University
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <vector>
#include "Config.h"
#include <numeric>
#include <random>
#include "Class_SiteLabel.h"
#include <ctime>
#include <algorithm>
#include <stdexcept>

namespace SSE{
    class Class_Lattice {
    private:
        //size of the lattice
        const type_ParaHamil ParaHamil;

        //Boundary Condition
        // Definition based on different models
        const type_DataInt Which_BC;

        const type_NumSite Num_Site;
        type_NumSite Num_Bond;
        const type_DataInt Which_CPU;
        type_DataFloat Ratio_Sum;

        // index_oper->site under this operator mapping
        std::vector<type_SiteinOper> Map_Oper_Site;

        // index_oper->weight mapping
        std::vector<SSE::type_DataFloat> Map_Oper_OperWeight;

        // index_oper->opertype mapping
        std::vector<SSE::type_DataInt> Map_Oper_Opertype;

        // index_oper->Detailedopertype mapping
        std::vector<SSE::type_DataInt> Map_Oper_Detailedtype;

        // index_oper sub-lattice mapping
        std::vector<SSE::type_DataInt> Map_Site_Sublattice;

        // (index_oper, index_operdiff) -> index_oper mpaaing
        std::vector<std::vector<SSE::type_DataInt>> Map_Site_DiffSite;

        // Some usefully array
        std::vector<SSE::type_DataInt> Map_Site_x1y0;
        std::vector<SSE::type_DataInt> Map_Site_x0y1;
        std::vector<SSE::type_DataInt> Map_Site_Qs;
        std::vector<SSE::type_DataInt> Map_Site_Qx;
        std::vector<SSE::type_DataInt> Map_Site_Qy;
        std::vector<SSE::type_DataInt> Array_Diasites;

    public:
        explicit Class_Lattice(type_ParaHamil _ParaHamil,
                               type_DataInt _which_cpu,
                               type_DataInt _bc,
                               type_DataInt _Num_Segment);

        type_DataInt Get_NumSite() const;
        type_DataFloat Get_Ratiosum() const;
        type_DataInt Get_CPU() const;
        type_DataInt Get_OperSite(type_NumSite _which_Index, type_NumSite _which_Site) const;
        const type_ParaHamil &Get_ParaHamil() const;
        const std::vector<SSE::type_DataInt> &Get_ArraySublattice() const;
        const std::vector<SSE::type_DataInt> &Get_Map_Site_x1y0() const;
        const std::vector<SSE::type_DataInt> &Get_Map_Site_x0y1() const;
        const std::vector<SSE::type_DataInt> &Get_Qs() const;
        const std::vector<SSE::type_DataInt> &Get_Qx() const;
        const std::vector<SSE::type_DataInt> &Get_Qy() const;
        SSE::type_NumSite Get_Sitediff(type_NumSite _which_site, type_NumSite _diff) const;
        const std::vector<SSE::type_DataInt> &Get_array_Sitediff(type_NumSite _diff) const;
        const std::vector<SSE::type_DataInt> &Get_array_Diasites() const;

        // About update
        type_DataInt Get_RandomBond() const;
        type_DataInt Get_BondType(type_NumSite _which_Index) const;
        type_DataInt Get_NumDetailedType(type_DataInt _which_DetailedType) const;
        type_DataInt Get_DetailedType(type_NumSite _which_index) const;

    };
}

inline SSE::type_DataInt SSE::Class_Lattice::Get_NumSite() const{
    return Num_Site;
}

inline SSE::type_DataFloat SSE::Class_Lattice::Get_Ratiosum() const {
    return Ratio_Sum;
}

inline SSE::type_DataInt SSE::Class_Lattice::Get_CPU() const {
    return Which_CPU;
}

inline SSE::type_DataInt SSE::Class_Lattice::Get_RandomBond() const{
    static std::mt19937_64 engine(static_cast<unsigned>(std::time(nullptr) + Which_CPU));
    static std::discrete_distribution<SSE::type_NumSite> Ran_Bond(Map_Oper_OperWeight.begin(), Map_Oper_OperWeight.end());

    return Ran_Bond(engine);
}

inline SSE::type_DataInt SSE::Class_Lattice::Get_BondType(type_NumSite _which_Index) const{
    return Map_Oper_Opertype[_which_Index];
}

inline SSE::type_DataInt SSE::Class_Lattice::Get_OperSite(type_NumSite _which_Index, type_NumSite _which_Site) const {
    return this->Map_Oper_Site[_which_Index][_which_Site];
}

inline SSE::type_DataInt SSE::Class_Lattice::Get_NumDetailedType(SSE::type_DataInt _which_DetailedType) const{
    return std::count(Map_Oper_Detailedtype.cbegin(), Map_Oper_Detailedtype.cend(), _which_DetailedType);;
}

inline const SSE::type_ParaHamil & SSE::Class_Lattice::Get_ParaHamil() const {
    return ParaHamil;
}

inline const std::vector<SSE::type_DataInt> &SSE::Class_Lattice::Get_ArraySublattice() const {
    return Map_Site_Sublattice;
}

inline const std::vector<SSE::type_DataInt> &SSE::Class_Lattice::Get_Map_Site_x1y0() const {
    return Map_Site_x0y1;
}

inline const std::vector<SSE::type_DataInt> &SSE::Class_Lattice::Get_Map_Site_x0y1() const {
    return Map_Site_x1y0;
}

inline const std::vector<SSE::type_DataInt> &SSE::Class_Lattice::Get_Qs() const {
    return Map_Site_Qs;
}

inline const std::vector<SSE::type_DataInt> &SSE::Class_Lattice::Get_Qx() const {
    return Map_Site_Qx;
}

inline const std::vector<SSE::type_DataInt> &SSE::Class_Lattice::Get_Qy() const {
    return Map_Site_Qy;
}

inline SSE::type_DataInt SSE::Class_Lattice::Get_DetailedType(SSE::type_DataInt _which_index) const {
    return Map_Oper_Detailedtype[_which_index];
}

inline SSE::type_NumSite SSE::Class_Lattice::Get_Sitediff(type_NumSite _which_site, type_NumSite _diff) const{
    return Map_Site_DiffSite[_which_site][_diff];
}

inline const std::vector<SSE::type_DataInt> &SSE::Class_Lattice::Get_array_Sitediff(SSE::type_NumSite _diff) const {
    return Map_Site_DiffSite[_diff];
}

inline const std::vector<SSE::type_DataInt> &SSE::Class_Lattice::Get_array_Diasites() const {
    return Array_Diasites;
}
