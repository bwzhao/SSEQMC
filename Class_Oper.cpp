////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author:          Bowen Zhao <bwzhao@bu.edu>
// Organization:    Boston University
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "Class_Oper.h"

SSE::Class_Oper::Class_Oper():
        Which_Index(-1),
        Which_Type(OPERTYPE_ID),
        Which_Status(0)
{}

bool SSE::Class_Oper::If_Suitable_andUpdateTypeOper(SSE::Class_Space &_which_Space, const SSE::Class_Lattice &_which_Lattice){
    if (this->Which_Type == SSE::OPERTYPE_t) {
        auto spin_0 = _which_Space.Get_Spin(_which_Lattice.Get_OperSite(this->Which_Index, 0));
        auto spin_1 = _which_Space.Get_Spin(_which_Lattice.Get_OperSite(this->Which_Index, 1));
        if (spin_0 == spin_1) {
            this->Which_Status = 0;
        }
        else{
            this->Which_Status = 1;
        }

        return true;
    }
    else if (this->Which_Type == SSE::OPERTYPE_J) {
        auto spin_0 = _which_Space.Get_Spin(_which_Lattice.Get_OperSite(this->Which_Index, 0));
        auto spin_1 = _which_Space.Get_Spin(_which_Lattice.Get_OperSite(this->Which_Index, 1));

        if (spin_0 == spin_1) {
            return false;
        }

        // Change the status to dia-dia
        this->Which_Status = 0;
        return true;
    }
    else {
        auto spin_0 = _which_Space.Get_Spin(_which_Lattice.Get_OperSite(this->Which_Index, 0));
        auto spin_1 = _which_Space.Get_Spin(_which_Lattice.Get_OperSite(this->Which_Index, 1));

        if (spin_0 == spin_1) {
            return false;
        }

        auto spin_2 = _which_Space.Get_Spin(_which_Lattice.Get_OperSite(this->Which_Index, 2));
        auto spin_3 = _which_Space.Get_Spin(_which_Lattice.Get_OperSite(this->Which_Index, 3));

        if (spin_2 == spin_3) {
            return false;
        }

        // Change the status to dia-dia
        this->Which_Status = 0;
        return true;
    }
}

bool SSE::Class_Oper::If_Dia_IfnotUpdateSpace(SSE::Class_Space &_which_Space, const SSE::Class_Lattice &_which_Lattice) const {
    if (this->Which_Type == SSE::OPERTYPE_t) {
        return true;
    }
    else if (this->Which_Type == SSE::OPERTYPE_J) {
        if (Which_Status == 0) {
            return true;
        }
        else {
            const auto index_Site_0 = _which_Lattice.Get_OperSite(this->Which_Index, 0);
            const auto index_Site_1 = _which_Lattice.Get_OperSite(this->Which_Index, 1);
            _which_Space.Flip_Spin(index_Site_0);
            _which_Space.Flip_Spin(index_Site_1);

            return false;
        }
    }
    else{
        if (Which_Status == 0) {
            return true;
        }
        else {
            if (Which_Status % 2 != 0) {
                const auto index_Site_0 = _which_Lattice.Get_OperSite(this->Which_Index, 0);
                const auto index_Site_1 = _which_Lattice.Get_OperSite(this->Which_Index, 1);
                _which_Space.Flip_Spin(index_Site_0);
                _which_Space.Flip_Spin(index_Site_1);
            }
            if (Which_Status / 2 != 0) {
                const auto index_Site_0 = _which_Lattice.Get_OperSite(this->Which_Index, 2);
                const auto index_Site_1 = _which_Lattice.Get_OperSite(this->Which_Index, 3);
                _which_Space.Flip_Spin(index_Site_0);
                _which_Space.Flip_Spin(index_Site_1);
            }
            return false;
        }
    }
}

bool SSE::Class_Oper::If_Dia_IfnotUpdateSpace_MeasureNxy(SSE::Class_Space &_which_Space, type_DataInt &_Nx,
                                                         type_DataInt & _Ny,
                                                         const SSE::Class_Lattice &_which_Lattice) const {
    if (this->Which_Type == SSE::OPERTYPE_t) {
        return true;
    }
    else if (this->Which_Type == SSE::OPERTYPE_J) {
        if (Which_Status == 0) {
            return true;
        }
        else {
            const auto index_Site_0 = _which_Lattice.Get_OperSite(this->Which_Index, 0);
            const auto index_Site_1 = _which_Lattice.Get_OperSite(this->Which_Index, 1);
            _which_Space.Flip_Spin(index_Site_0);
            _which_Space.Flip_Spin(index_Site_1);
            if (_which_Lattice.Get_DetailedType(Which_Index) == DETAILEDTYPE_Jx) {
                _Nx +=  _which_Space.Get_Spin(index_Site_1);
            }
            else{
                _Ny += _which_Space.Get_Spin(index_Site_1);
            }

            return false;
        }
    }
    else{
        if (Which_Status == 0) {
            return true;
        }
        else {
            if (Which_Status % 2 != 0) {
                const auto index_Site_0 = _which_Lattice.Get_OperSite(this->Which_Index, 0);
                const auto index_Site_1 = _which_Lattice.Get_OperSite(this->Which_Index, 1);
                _which_Space.Flip_Spin(index_Site_0);
                _which_Space.Flip_Spin(index_Site_1);
                if (_which_Lattice.Get_DetailedType(Which_Index) == DETAILEDTYPE_Q2x) {
                    _Nx += _which_Space.Get_Spin(index_Site_1);
                }
                else{
                    _Ny += _which_Space.Get_Spin(index_Site_1);
                }
            }
            if (Which_Status / 2 != 0) {
                const auto index_Site_0 = _which_Lattice.Get_OperSite(this->Which_Index, 2);
                const auto index_Site_1 = _which_Lattice.Get_OperSite(this->Which_Index, 3);
                _which_Space.Flip_Spin(index_Site_0);
                _which_Space.Flip_Spin(index_Site_1);

                if (_which_Lattice.Get_DetailedType(Which_Index) == DETAILEDTYPE_Q2x) {
                    _Nx += _which_Space.Get_Spin(index_Site_1);
                }
                else{
                    _Ny += _which_Space.Get_Spin(index_Site_1);
                }
            }
            return false;
        }
    }
}

SSE::type_DataInt SSE::Class_Oper::Update_Oper(SSE::type_DataInt _which_Leg, bool _ifFlip) {


    if (this->Which_Type == SSE::OPERTYPE_t) {
        auto other_leg =  SSE::T_t[_which_Leg];
        if (_ifFlip) {
            Which_Status = (Which_Status + 1) % 2;
        }
        return other_leg;
    }
        // Depend on different operators
    else if (this->Which_Type == SSE::OPERTYPE_J){

        auto other_leg = SSE::T_J[_which_Leg];
        if (_ifFlip) {
            Which_Status = (Which_Status + 1) % 2;
        }
        return other_leg;
    }
    else{
        auto other_leg = SSE::T_Q2[_which_Leg];
        if (_ifFlip) {
            if (_which_Leg <= 3){
                Which_Status = (Which_Status + 1) % 2 + Which_Status / 2 * 2;
            }
            else{
                Which_Status = (Which_Status + 2) % 4;
            }

        }
        return other_leg;
    }
}

void SSE::Class_Oper::Update_fMat(SSE::Class_fMat &_which_fMat, const SSE::Class_Lattice &_which_Lattice) const {
    if (this->Which_Type == SSE::OPERTYPE_J) {
        const auto index_Site_0 = _which_Lattice.Get_OperSite(this->Which_Index, 0);
        const auto index_Site_1 = _which_Lattice.Get_OperSite(this->Which_Index, 1);
        //If it is a diagonal operator
        if (Which_Status == 0) {
            _which_fMat.Set_RowZero(index_Site_0);
            _which_fMat.Set_RowZero(index_Site_1);
        }
        else{
            auto temp_weight = 2. * _which_Lattice.Get_ParaHamil().tx / _which_Lattice.Get_ParaHamil().Jx;
            _which_fMat.Multi_Row(index_Site_0, temp_weight);
            _which_fMat.Multi_Row(index_Site_1, temp_weight);
            _which_fMat.Swap_Row(index_Site_0, index_Site_1);
        }
    }
    else if (this->Which_Type == SSE::OPERTYPE_t) {
        // If it is the up/up matrix
        if (Which_Status == 0) {
            const auto index_Site_0 = _which_Lattice.Get_OperSite(this->Which_Index, 0);
            const auto index_Site_1 = _which_Lattice.Get_OperSite(this->Which_Index, 1);

            _which_fMat.SwapMinus_Row(index_Site_0, index_Site_1);
        }
    }

}

void SSE::Class_Oper::Update_fVec(std::valarray<type_DataFloat> &_which_fVec,
                                  const SSE::Class_Lattice &_which_Lattice) const {
    if (this->Which_Type == SSE::OPERTYPE_J) {
        const auto index_Site_0 = _which_Lattice.Get_OperSite(this->Which_Index, 0);
        const auto index_Site_1 = _which_Lattice.Get_OperSite(this->Which_Index, 1);
        //If it is a diagonal operator
        if (Which_Status == 0) {
            _which_fVec[index_Site_0] = 0;
            _which_fVec[index_Site_1] = 0;
        }
        else{
            auto temp_weight = _which_Lattice.Get_ParaHamil().tx / (_which_Lattice.Get_ParaHamil().Jx / 2 +
                    2 * _which_Lattice.Get_ParaHamil().Q2 / 4);
            auto temp_val_0 = _which_fVec[index_Site_0] * temp_weight;
            auto temp_val_1 = _which_fVec[index_Site_1] * temp_weight;

            _which_fVec[index_Site_0] = temp_val_1;
            _which_fVec[index_Site_1] = temp_val_0;
        }
    }
    else if (this->Which_Type == SSE::OPERTYPE_t) {
        // If it is the up/up matrix
        if (Which_Status == 0) {
            const auto index_Site_0 = _which_Lattice.Get_OperSite(this->Which_Index, 0);
            const auto index_Site_1 = _which_Lattice.Get_OperSite(this->Which_Index, 1);

            auto temp_val_0 = _which_fVec[index_Site_0];
            auto temp_val_1 = _which_fVec[index_Site_1];

            _which_fVec[index_Site_0] = temp_val_0 - temp_val_1;
            _which_fVec[index_Site_1] = temp_val_1 - temp_val_0;
        }
        // up/down case, still the same
    }
    // Q2 oper
    else{
        if ((Which_Status == 0) or (Which_Status == 3)) {
            const auto index_Site_0 = _which_Lattice.Get_OperSite(this->Which_Index, 0);
            const auto index_Site_1 = _which_Lattice.Get_OperSite(this->Which_Index, 1);
            const auto index_Site_2 = _which_Lattice.Get_OperSite(this->Which_Index, 2);
            const auto index_Site_3 = _which_Lattice.Get_OperSite(this->Which_Index, 3);

            _which_fVec[index_Site_0] = 0.;
            _which_fVec[index_Site_1] = 0.;
            _which_fVec[index_Site_2] = 0.;
            _which_fVec[index_Site_3] = 0.;
        }
        else {
            if (Which_Status % 2 != 0) {
                const auto index_Site_0 = _which_Lattice.Get_OperSite(this->Which_Index, 0);
                const auto index_Site_1 = _which_Lattice.Get_OperSite(this->Which_Index, 1);

                auto temp_weight = _which_Lattice.Get_ParaHamil().tx / (_which_Lattice.Get_ParaHamil().Jx / 2 +
                                                                        2 * _which_Lattice.Get_ParaHamil().Q2 / 4);
                auto temp_val_0 = _which_fVec[index_Site_0] * temp_weight;
                auto temp_val_1 = _which_fVec[index_Site_1] * temp_weight;

                _which_fVec[index_Site_0] = temp_val_1;
                _which_fVec[index_Site_1] = temp_val_0;
            }
            else {
                const auto index_Site_0 = _which_Lattice.Get_OperSite(this->Which_Index, 2);
                const auto index_Site_1 = _which_Lattice.Get_OperSite(this->Which_Index, 3);
                auto temp_weight = _which_Lattice.Get_ParaHamil().tx / (_which_Lattice.Get_ParaHamil().Jx / 2 +
                                                                        2 * _which_Lattice.Get_ParaHamil().Q2 / 4);
                auto temp_val_0 = _which_fVec[index_Site_0] * temp_weight;
                auto temp_val_1 = _which_fVec[index_Site_1] * temp_weight;

                _which_fVec[index_Site_0] = temp_val_1;
            }
        }
    }
}