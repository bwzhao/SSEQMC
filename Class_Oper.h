////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author:          Bowen Zhao <bwzhao@bu.edu>
// Organization:    Boston University
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once
#include "Config.h"
#include "Class_Lattice.h"
#include "Class_Space.h"

namespace SSE{
    class Class_Oper {
    private:
        // Which index of the operator
        type_DataInt Which_Index;

        // Which Type of the operator
        type_DataInt Which_Type;

        // Which status of this operator
        // Detailed definition, please see Lattice.h
        type_DataInt Which_Status;

    public:
        explicit Class_Oper();
        ~Class_Oper() = default;

        // If it is an identity operator
        bool If_Identity() const;
        bool If_Suitable_andUpdateTypeOper(SSE::Class_Space &_which_Space, const SSE::Class_Lattice &_which_Lattice);
        bool If_Dia_IfnotUpdateSpace(SSE::Class_Space &_which_Space, const SSE::Class_Lattice &_which_Lattice) const;
        bool If_Dia_IfnotUpdateSpace_MeasureNxy(SSE::Class_Space &_which_Space, type_DataInt &_Nx, type_DataInt &_Ny, const SSE::Class_Lattice &_which_Lattice) const;

        void Set_Oper_temp(SSE::type_NumSite _which_Index, const SSE::Class_Lattice &_which_Lattice);
        void Reset_fromDiatoID();

        type_DataInt Get_NumSites() const;
        type_DataInt Get_Index() const;
        type_DataInt Get_Type() const;

        //Update the operator in Loop update process
        type_DataInt Update_Oper(type_DataInt _which_Leg, bool _ifFlip);

    };
}

inline bool SSE::Class_Oper::If_Identity() const{
    return this->Which_Type == SSE::OPERTYPE_ID;
}



inline void SSE::Class_Oper::Set_Oper_temp(SSE::type_NumSite _which_Index, const SSE::Class_Lattice &_which_Lattice) {
    this->Which_Index = _which_Index;
    this->Which_Type = _which_Lattice.Get_BondType(_which_Index);
    this->Which_Status = 0;
    // Don't need to set other terms for now.
}

inline void SSE::Class_Oper::Reset_fromDiatoID() {
    Which_Index = -1;
    Which_Type = OPERTYPE_ID;
//    Which_Status = 0;
}

inline SSE::type_DataInt SSE::Class_Oper::Get_NumSites() const {
    return ARRAY_NUMSITE[Which_Type];
}

inline SSE::type_DataInt SSE::Class_Oper::Get_Index() const {
    return Which_Index;
}

inline SSE::type_DataInt SSE::Class_Oper::Get_Type() const {
    return Which_Type;
}






