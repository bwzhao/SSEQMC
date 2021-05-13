////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author:          Bowen Zhao <bwzhao@bu.edu>
// Organization:    Boston University
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once
#include "Config.h"
#include "Class_Lattice.h"
#include "Class_Space.h"
#include <fstream>
#include "Class_fMat.h"
#include <ctime>

namespace SSE{
    class Class_Oper {
    private:
        // Which index of the operator
        type_DataInt Which_Index;

        // Which Type of the operator
        type_DataInt Which_Type;

        // Which status of this operator
        // Detailed definition, please see Lattice.h
        // For t oper: 0-> up-up or down-down; 1-> up-down or down-up
        // For J oper: 0-> dia; 1-> off-dia
        // For Q2 oper: (Which_Status / Which_part % 2): 0-> dia; 1-> off-dia
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
        type_DataInt Get_NumLegs() const;
        type_DataInt Get_Index() const;
        type_DataInt Get_Type() const;

        //Update the operator in Loop update process
        type_DataInt Update_Oper(type_DataInt _which_Leg, bool _ifFlip);
        void Update_fMat(SSE::Class_fMat &_which_fMat, const SSE::Class_Lattice &_which_Lattice) const;
        void Update_fVec(std::valarray<type_DataFloat> &_which_fVec, const SSE::Class_Lattice &_which_Lattice) const;

        // Write the list:
        void Write_Oper(std::ofstream &_which_file);
        void Read_Oper(std::ifstream &_which_file, const SSE::Class_Lattice &_which_Lattice);

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

inline SSE::type_DataInt SSE::Class_Oper::Get_NumLegs() const {
    return ARRAY_NUMSITE[Which_Type] * NUM_LEG_PER_SITE;
}

inline SSE::type_DataInt SSE::Class_Oper::Get_Index() const {
    return Which_Index;
}

inline SSE::type_DataInt SSE::Class_Oper::Get_Type() const {
    return Which_Type;
}

inline void SSE::Class_Oper::Write_Oper(std::ofstream &_which_file) {
    _which_file.write(reinterpret_cast<char *>(&this->Which_Index), sizeof(decltype(this->Which_Index)));
    _which_file.write(reinterpret_cast<char *>(&Which_Status), sizeof(decltype(this->Which_Status)));
}

inline void SSE::Class_Oper::Read_Oper(std::ifstream &_which_file, const SSE::Class_Lattice &_which_Lattice) {
    _which_file.read(reinterpret_cast<char *>(&this->Which_Index), sizeof(decltype(this->Which_Index)));
    _which_file.read(reinterpret_cast<char *>(&this->Which_Status), sizeof(decltype(this->Which_Status)));
    if (this->Which_Index != -1) {
        this->Set_Oper_temp(this->Which_Index, _which_Lattice);
    }
    else{
        this->Which_Index = SSE::OPERTYPE_ID;
        this->Which_Type = DETAILEDTYPE_ID;
    }
}








