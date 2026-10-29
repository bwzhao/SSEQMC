////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author:          Bowen Zhao <bwzhao@bu.edu>
// Organization:    Boston University
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "Class_SSE.h"
#include <cstdint>
#include <cstdio>
#include <iterator>


SSE::Class_SSE::Class_SSE(SSE::type_ParaHamil _para_Hamil, SSE::type_DataInt _which_cpu,
                          SSE::type_DataInt _bc,
                          const std::string &_file_config, const std::string &_file_data, const std::string &_file_corr,
                          SSE::type_DataInt _num_Segment, SSE::type_DataInt _num_TimeSlices):
        ThisLattice(_para_Hamil, _which_cpu, _bc, _num_Segment),
        Num_Segment(_num_Segment),
        Array_Oper(_num_Segment, std::vector<SSE::Class_Oper>(std::max(4 * static_cast<type_DataInt>(ThisLattice.Get_Ratiosum()), MIN_INITIAL_CUTOFF))),
        Array_NumOper(_num_Segment, 0),
        Num_Measure_TimeSlices(_num_TimeSlices),
        Str_FileConfig(_file_config),
        Str_FileData(_file_data),
        Str_FileCorr(_file_corr),
        Str_NameClass("SSE_"),
        Array_SpaceLinkList(ThisLattice.Get_NumSite(), std::vector<class_LegLink>(2)),
        Array_OperLinkList(_num_Segment),
        Array_OperMapping(_num_Segment)
{
    // Seed the shared random engine before anything draws from it
    SSE::Seed_Engine(_which_cpu);

    // Initialize Space
    ThisSpace.Set_Class_Space(ThisLattice);

    // Initialize nameclass
    char temp_chars_nameclass[200];
    snprintf(temp_chars_nameclass, 200, "%.6f_%.6f_%.6f_%.6f_%.6f_%.6f_%.6f_%.6f_%.6f_%.6f_%d_%d_%lf_%d_%d_%d",
             _para_Hamil.ty,
             _para_Hamil.tx,
             _para_Hamil.Jy,
             _para_Hamil.Jx,
             _para_Hamil.Q2,
             _para_Hamil.Q2s,
             _para_Hamil.JfA,
             _para_Hamil.JfB,
             _para_Hamil.Q3,
             _para_Hamil.Q3s,
             _para_Hamil.Num_X,
             _para_Hamil.Num_Y,
             _para_Hamil.beta,
             _which_cpu,
             _bc,
             _num_Segment);
    Str_NameFile = std::string(temp_chars_nameclass);

    std::cout << "Successfully Create a Calculate Class:" << std::endl;
    std::cout << "#parameter:" << std::endl
              << "#Haml_tA = " << _para_Hamil.ty << std::endl
              << "#Haml_tB = " << _para_Hamil.tx << std::endl
              << "#Haml_JA = " << _para_Hamil.Jy << std::endl
              << "#Haml_JB = " << _para_Hamil.Jx << std::endl
              << "#Haml_Q2 = " << _para_Hamil.Q2 << std::endl
              << "#Haml_Q2s = " << _para_Hamil.Q2s << std::endl
              << "#Haml_JfA = " << _para_Hamil.JfA << std::endl
              << "#Haml_JfB = " << _para_Hamil.JfB << std::endl
              << "#Haml_Q3 = " << _para_Hamil.Q3 << std::endl
              << "#Haml_Q3s = " << _para_Hamil.Q3s << std::endl
              << "#Lx = " << _para_Hamil.Num_X << std::endl
              << "#Ly = " << _para_Hamil.Num_Y << std::endl
              << "#beta = " << _para_Hamil.beta << std::endl
              << "#BC = " << _bc << std::endl
              << "#Num_Segment = " << _num_Segment << std::endl
              << "#Dir_Config = " << _file_config << std::endl
              << "#Dir_Data = " << _file_data << std::endl
              << "#Dir_Corr = " << _file_corr << std::endl;

    // Initializing measurement
    SSE::Class_DataMeasurement<SSE::type_DataFloat> temp_Class_DataMeasurement;
    for (const auto &which_key: SSE::Vec_Name) {
        MeasureData.emplace_back(temp_Class_DataMeasurement);
        for (const auto &which_operkey: SSE::Vec_OperNum) {
            MeasureData.emplace_back(temp_Class_DataMeasurement);
        }
    }
    // Plus the operator numbers themselves, stored after all the quantities in Measure_woDynamics
    for (const auto &which_operkey: SSE::Vec_OperNum) {
        MeasureData.emplace_back(temp_Class_DataMeasurement);
    }

    // Initialize f measurement
    SSE::Class_CorrMeasurement<SSE::type_DataFloat> temp_Array_Corr;
    for (SSE::type_NumSite index_tau = 0; index_tau != this->Num_Measure_TimeSlices; ++index_tau) {
        Corr_ftau.emplace_back(temp_Array_Corr);
    }
}


void SSE::Class_SSE::DiagonalUpdate() {
    // Random number generator
    auto & engine = SSE::Get_Engine();
    static std::uniform_real_distribution<SSE::type_DataFloat> Ran_Double(0., 1.);
#ifndef NDEBUG
    auto temp_space = ThisSpace;
#endif

    // Propagating space
    // Loop through segments
    for (type_DataInt index_Segment = 0; index_Segment != Num_Segment; ++index_Segment) {
        // Loop through operators in each segment
        auto & which_Segment = Array_Oper[index_Segment];
        auto & current_NumOper_Segment = Array_NumOper[index_Segment];
        const type_DataInt & max_NumOper_Segment = which_Segment.size();

        for (auto & which_oper : which_Segment) {
            // If it is an identity operator: ID -> OP
            if (which_oper.If_Identity()) {
                auto ratio_ID_to_OP = ThisLattice.Get_Ratiosum() / (max_NumOper_Segment - current_NumOper_Segment);
                // If it can update
                if(Ran_Double(engine) < ratio_ID_to_OP){
                    auto which_bond = ThisLattice.Get_RandomBond();
                    which_oper.Set_Oper_temp(which_bond, ThisLattice);

                    if (which_oper.If_Suitable_andUpdateTypeOper(ThisSpace, ThisLattice)) {
                        ++current_NumOper_Segment;
                    }
                    else{
                        which_oper.Reset_fromDiatoID();
                    }
                }
                // If it cannot update, NOT NECESSARY!!!
//                else{
//                    which_oper.Reset_fromDiatoID();
//                }
            }
                // If it is a real operator: OP -> ID
            else {
                // If it is a diagonal operator
                if (which_oper.If_Dia_IfnotUpdateSpace(ThisSpace, ThisLattice)) {
                    auto ratio_OP_to_ID = (max_NumOper_Segment - current_NumOper_Segment + 1) / ThisLattice.Get_Ratiosum();
                    // If this oper can update back to identity
                    if(Ran_Double(engine) < ratio_OP_to_ID){
                        which_oper.Reset_fromDiatoID();
                        --current_NumOper_Segment;
                    }
//                    else{}
                }
//                 If it is an off-diagonal operator
//                else{}
            }
        }
    }

    //Just for debug
#ifndef NDEBUG
    if (!temp_space.Is_Same(ThisSpace)) {
        throw std::runtime_error("Error! In Diagonal Update, The propgated space and right space are not the same");
    }
#endif

}

void SSE::Class_SSE::Generate_Linklist() {
    std::vector<SSE::class_LegLink> first_leg_array(ThisLattice.Get_NumSite(), SSE::class_LegLink());
    std::vector<SSE::class_LegLink> last_leg_array(ThisLattice.Get_NumSite(), SSE::class_LegLink());

    //for every time slide
    // Loop through segments
    for (SSE::type_DataInt index_Segment = 0; index_Segment != Num_Segment; ++index_Segment) {
        // Loop through operators in each segment
        auto & which_Segment = Array_Oper[index_Segment];
        auto & which_OperMapping_Segment = Array_OperMapping[index_Segment];
        auto & which_OperLink_list_Segment = Array_OperLinkList[index_Segment];

        // Remember it is the index in the pure operator list
        SSE::type_DataInt index_oper_inSegment = 0;

//        for (SSE::type_DataInt index_oper_inSegment = 0; index_oper_inSegment != which_Segment.size(); ++index_oper_inSegment) {
        for (auto & which_oper : which_Segment){
            // Exclude identity operators
            if (which_oper.If_Identity()) {
                continue;
            }
            which_OperMapping_Segment[index_oper_inSegment] = & which_oper;

            auto & temp_Operlink = which_OperLink_list_Segment[index_oper_inSegment];

            for (SSE::type_DataInt index_site_inBond = 0; index_site_inBond != which_oper.Get_NumSites(); ++index_site_inBond) {
                auto index_this_site = ThisLattice.Get_OperSite(which_oper.Get_Index(), index_site_inBond);

                auto this_site_left_leg = NUM_LEG_PER_SITE * index_site_inBond + 1;
                auto & temp_LeftLink = temp_Operlink[this_site_left_leg];
                auto this_site_right_leg = NUM_LEG_PER_SITE * index_site_inBond + 0;
                auto & temp_RightLink = temp_Operlink[this_site_right_leg];

                const auto & temp_last_leg = last_leg_array[index_this_site];

                if (temp_last_leg.If_Unlinked()) {        //it is the first leg
                    first_leg_array[index_this_site].Set_LegLink(index_Segment, index_oper_inSegment, this_site_left_leg);
                }
                else {                            //it is not the first leg
//                    temp_LeftLink.Set_LegLink(temp_last_leg);
                    temp_LeftLink = temp_last_leg;

                    auto & LegLink_linked = Array_OperLinkList[temp_last_leg.Which_Segment][temp_last_leg.Which_Time][temp_last_leg.Which_smLeg];
                    LegLink_linked.Set_LegLink(index_Segment, index_oper_inSegment, this_site_left_leg);
                }
                last_leg_array[index_this_site].Set_LegLink(index_Segment, index_oper_inSegment, this_site_right_leg);
            }

            // Counting plus one
            ++index_oper_inSegment;
        }
    }

    //link the last leg and the first leg
    for (SSE::type_NumSite index_site = 0; index_site != ThisLattice.Get_NumSite(); ++index_site) {
        auto & this_site_first_leg = first_leg_array[index_site];
        auto & this_site_last_leg = last_leg_array[index_site];

        auto & temp_LeftLink = Array_SpaceLinkList[index_site][1];
        auto & temp_RightLink = Array_SpaceLinkList[index_site][0];

        //this site the space is connected with operator
        if (!this_site_first_leg.If_Unlinked()) {
            //link the leg with the spin
            Array_OperLinkList[this_site_first_leg.Which_Segment][this_site_first_leg.Which_Time]
            [this_site_first_leg.Which_smLeg].Set_LegLink(-2 - 0, index_site, 0);

            temp_RightLink.Set_LegLink(this_site_first_leg);

            Array_OperLinkList[this_site_last_leg.Which_Segment][this_site_last_leg.Which_Time]
                    [this_site_last_leg.Which_smLeg].Set_LegLink(-2 - 1, index_site, 1);
            temp_LeftLink.Set_LegLink(this_site_last_leg);
        }
            //not connected with operator
        else{
            temp_LeftLink.Reset_LegLink();
            temp_RightLink.Reset_LegLink();
        }
    }
}

void SSE::Class_SSE::Flip_Update() {
    auto & engine = SSE::Get_Engine();
    static std::uniform_real_distribution<SSE::type_DataFloat> Ran_Double(0., 1.);

    for (int nl_site = 0; nl_site != ThisLattice.Get_NumSite(); ++nl_site) {
        if (this->If_CanFlip(nl_site)) {
            const auto FlipRate = Ran_Double(engine);
            if (FlipRate <= 0.5) {
                ThisSpace.Flip_Spin(nl_site);
            }
        }
    }
}

void SSE::Class_SSE::Adjust_Cutoff() {
    auto sum_NumOper = std::accumulate(Array_NumOper.begin(), Array_NumOper.end(), 0.);
    SSE::type_DataInt max_NumOper_inSegment = *(std::max_element(Array_NumOper.begin(), Array_NumOper.end()));
    SSE::type_DataInt NumTime_inSegment = Array_Oper.back().size();

    for (SSE::type_DataInt index_Segment = 0; index_Segment != Num_Segment; ++index_Segment) {
        Array_OperLinkList[index_Segment].resize(Array_NumOper[index_Segment]);
        Array_OperMapping[index_Segment].resize(Array_NumOper[index_Segment]);
    }

    if (((max_NumOper_inSegment * 1.2) > NumTime_inSegment) || (sum_NumOper * 1.5 / Num_Segment > NumTime_inSegment) ) {
        auto new_NumTime = std::max(static_cast<SSE::type_DataInt>(sum_NumOper * 1.5 / Num_Segment),
                                    static_cast<SSE::type_DataInt>(max_NumOper_inSegment * 1.2)) + 1;
//        auto Oper_Array_Size = static_cast<int>(Oper_Loc_Array.size());
        for (SSE::type_DataInt index_Segment = 0; index_Segment != Num_Segment; ++index_Segment) {
            Array_Oper[index_Segment].resize(new_NumTime);
        }
//        for (auto & which_Segment : Array_Oper) {
//            which_Segment.resize(new_NumTime);
//        }
    }
}

void SSE::Class_SSE::Loop_Update() {
    //Choose a random color
    auto & engine = SSE::Get_Engine();
    static std::uniform_real_distribution<SSE::type_DataFloat> Ran_Double(0., 1.);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Loop through all the legs in each segment
    // For each segment
    for (type_DataInt index_start_Segment = 0; index_start_Segment != Num_Segment; ++index_start_Segment) {
        auto &which_start_Segment_Oper = Array_OperMapping[index_start_Segment];
        auto &which_start_Segment_LegLinkList = Array_OperLinkList[index_start_Segment];
        // For each operators
        for (type_DataInt index_start_oper = 0; index_start_oper != which_start_Segment_Oper.size(); ++index_start_oper) {
            auto &which_start_oper = *(which_start_Segment_Oper[index_start_oper]);

            auto & which_StartOperLink = which_start_Segment_LegLinkList[index_start_oper];
            // For each legs
            for (SSE::type_DataInt index_start_leg = 0; index_start_leg != which_start_oper.Get_NumLegs(); ++index_start_leg) {
                // Only consider the unupdated legs
                if (which_StartOperLink[index_start_leg].If_Unupdated()) {
                    SSE::class_LegLink COPY_StartLeg = {index_start_Segment, index_start_oper, index_start_leg};

                    bool IfFlip = true;
                    if (Ran_Double(engine) <= 0.5){
                        IfFlip = false;
                    };
                    // Starting from one unupdated-leg and update all the legs within a loop.

                    // The loop will stop at the initial leg
                    auto Insmleg_Start = index_start_leg;
                    auto Outsmleg_Start = which_start_oper.Update_Oper(Insmleg_Start, IfFlip);

                    SSE::class_LegLink current_Leg = which_StartOperLink[Outsmleg_Start];

                    // Reset the linkList
                    which_StartOperLink[Insmleg_Start].Reset_LegLink();
                    which_StartOperLink[Outsmleg_Start].Reset_LegLink();

                    // This one must be linked to an operator
                    current_Leg = Find_Nextinleg_UpdateSpace(current_Leg, IfFlip);
                    while (current_Leg != COPY_StartLeg) {
                        auto & current_oper = *(Array_OperMapping[current_Leg.Which_Segment][current_Leg.Which_Time]);

                        auto Insmleg_Current = current_Leg.Which_smLeg;
                        auto Outsmleg_Current = current_oper.Update_Oper(Insmleg_Current, IfFlip);

                        auto & current_OperLink = Array_OperLinkList[current_Leg.Which_Segment][current_Leg.Which_Time];

                        current_Leg = current_OperLink[Outsmleg_Current];

                        // Reset the linkList
                        current_OperLink[Insmleg_Current].Reset_LegLink();
                        current_OperLink[Outsmleg_Current].Reset_LegLink();

                        current_Leg = Find_Nextinleg_UpdateSpace(current_Leg, IfFlip);
                    }
                }
            }
        }
    }
}

SSE::class_LegLink SSE::Class_SSE::Find_Nextinleg_UpdateSpace(const SSE::class_LegLink &_which_LegLink,
                                                              bool _ifFlip) {
    SSE::class_LegLink temp_return;
    //If it is linked to an actual operator
    if (_which_LegLink.If_OperLink()) {
        return _which_LegLink;
    }
        //If it is linked to a space, and only for periodic boundary condition
    else{
        const auto & which_site = _which_LegLink.Which_Time;

        if (_ifFlip) {
            ThisSpace.Flip_Spin(which_site);
        }

        return Array_SpaceLinkList[which_site][(_which_LegLink.Which_smLeg + 1) % 2];

        // Not necessary for reset the leglink
//        which_Leftspace.Get_LegLink(which_site, INDEX_LEFTSPACE).Reset_LegLink();
//        which_Rightspace.Get_LegLink(which_site, INDEX_RIGHTSPACE).Reset_LegLink();
    }
}

void SSE::Class_SSE::CheckLinkList() const {
    for (type_DataInt index_start_Segment = 0; index_start_Segment != Num_Segment; ++index_start_Segment) {
        auto &which_start_Segment_LegLinkList = Array_OperLinkList[index_start_Segment];
        if (Array_NumOper[index_start_Segment] != Array_OperMapping[index_start_Segment].size()){
            throw std::runtime_error("Oper Number Error");
        }
        for (type_DataInt index_start_oper = 0;
             index_start_oper != which_start_Segment_LegLinkList.size(); ++index_start_oper) {
            auto &which_StartOperLink = which_start_Segment_LegLinkList[index_start_oper];
            // For each legs
            for (SSE::type_DataInt index_start_leg = 0;
                 index_start_leg != which_StartOperLink.size(); ++index_start_leg) {
                auto & Which_LegLinkList = which_StartOperLink[index_start_leg];
                if (Which_LegLinkList.Which_Segment >= 0){
                    auto & Other_LegLinkList = Array_OperLinkList[Which_LegLinkList.Which_Segment][Which_LegLinkList.Which_Time][Which_LegLinkList.Which_smLeg];
                    if ((Other_LegLinkList.Which_Time != index_start_oper) | (Other_LegLinkList.Which_smLeg != index_start_leg)) {
                        throw std::runtime_error("Oper Link Error");
                    }
                }
                else if (Which_LegLinkList.Which_Segment == -1) {
                    throw std::runtime_error("Not fully linked");
                }
                else{ // If Segment <= -2
                    const auto & WhichOper = *(Array_OperMapping[index_start_Segment][index_start_oper]);
                    if ((ThisLattice.Get_OperSite(WhichOper.Get_Index(), index_start_leg / 2)) !=
                        Which_LegLinkList.Which_Time) {
                        throw std::runtime_error("Wrong Sites one Oper");
                    }
                    const auto Ohter_LegLinkList = Array_SpaceLinkList[Which_LegLinkList.Which_Time][Which_LegLinkList.Which_smLeg];
                    if ((Ohter_LegLinkList.Which_Time != index_start_oper) | (Ohter_LegLinkList.Which_smLeg != index_start_leg)){
                        throw std::runtime_error("Oper-Space Link Error");
                    }
                }
            }
        }
    }
}

void SSE::Class_SSE::Measure_woDynamics() {
    // store the temp data
    std::vector<SSE::Class_DataMeasurement<SSE::type_DataFloat>> temp_MeasureData;

    for (const auto &which_key: SSE::Vec_Name) {
        temp_MeasureData.emplace_back(SSE::Class_DataMeasurement<SSE::type_DataFloat>());
    }
    std::vector<SSE::type_DataInt> temp_NumOper;
    for (const auto &which_key: SSE::Vec_OperNum) {
        temp_NumOper.emplace_back(0);
    }

    // Measure_wSzDynamics Energy, which is the number of operators
    auto sum_NumOper = std::accumulate(Array_NumOper.begin(), Array_NumOper.end(), 0.);
    auto Num_tyOper = ThisLattice.Get_NumDetailedType(DETAILEDTYPE_ty);
    auto Num_txOper = ThisLattice.Get_NumDetailedType(DETAILEDTYPE_tx);

    auto temp_Energy = -(static_cast<SSE::type_DataFloat>(sum_NumOper) / ThisLattice.Get_ParaHamil().beta) /
                       static_cast<SSE::type_DataFloat>(ThisLattice.Get_NumSite())
                       + Num_tyOper * ThisLattice.Get_ParaHamil().ty / static_cast<SSE::type_DataFloat>(ThisLattice.Get_NumSite())
                       + Num_txOper * ThisLattice.Get_ParaHamil().tx / static_cast<SSE::type_DataFloat>(ThisLattice.Get_NumSite());


    // Count the operator number
    SSE::type_DataInt count_time = 0;
    auto Num_Oper = std::vector<SSE::type_DataInt>(Vec_OperNum.size(), 0);
    SSE::type_DataInt Nx = 0;
    SSE::type_DataInt Ny = 0;

    // Loop through segments
    for (type_DataInt index_Segment = 0; index_Segment != Num_Segment; ++index_Segment) {
        // Loop through operators in each segment
        auto &which_operSegment = Array_Oper[index_Segment];

        // Measure the order parameters
        for (auto &which_oper : which_operSegment) {
            if (!which_oper.If_Identity()) {
                Num_Oper[which_oper.Get_Type()] += 1;

                if (count_time % ThisLattice.Get_NumSite() == 0) {
                    const auto temp_mz = ThisSpace.Measure_mz(ThisLattice);
                    temp_MeasureData[SSE::INDEX_MZ2].AppendValue(temp_mz * temp_mz);
                    temp_MeasureData[SSE::INDEX_MZ4].AppendValue(temp_mz * temp_mz * temp_mz * temp_mz);

                    const auto &temp_md = ThisSpace.Measure_md(ThisLattice);
                    auto temp_mc2 = temp_md[0] * temp_md[0] + temp_md[1] * temp_md[1];
                    auto temp_ms2 = temp_md[2] * temp_md[2] + temp_md[3] * temp_md[3];

                    if (temp_mc2 != 0) {
                        auto temp_costheta_c = temp_md[0] / sqrt(temp_mc2);
                        auto temp_Phi4_c =
                                8 * temp_costheta_c * temp_costheta_c * temp_costheta_c * temp_costheta_c
                                - 8 * temp_costheta_c * temp_costheta_c + 1;
                        temp_MeasureData[SSE::INDEX_PHI4C].AppendValue(temp_Phi4_c);
                    }
//
                    temp_MeasureData[SSE::INDEX_MC2].AppendValue(temp_mc2);
                    temp_MeasureData[SSE::INDEX_MC4].AppendValue(temp_mc2 * temp_mc2);

                    temp_MeasureData[SSE::INDEX_MS2].AppendValue(temp_ms2);
                    temp_MeasureData[SSE::INDEX_MS4].AppendValue(temp_ms2 * temp_ms2);
                }
                // Just use this function to update the space
                which_oper.If_Dia_IfnotUpdateSpace_MeasureNxy(ThisSpace, Nx, Ny, ThisLattice);
                ++count_time;
            }
        }
    }

    temp_MeasureData[SSE::INDEX_E].AppendValue(temp_Energy);
    temp_MeasureData[SSE::INDEX_N].AppendValue(sum_NumOper);
    temp_MeasureData[SSE::INDEX_N2].AppendValue(sum_NumOper * sum_NumOper);

    // Measure rho_s and chi_u
    SSE::type_DataFloat temp_Wx2 = Nx*Nx;
    SSE::type_DataFloat temp_Wy2 = Ny*Ny;

    const auto val_rho_s = (temp_Wx2 + temp_Wy2) / (ThisLattice.Get_ParaHamil().beta * 2);

    temp_MeasureData[SSE::INDEX_RHOS].AppendValue(val_rho_s);

    const auto val_mu = ThisSpace.Measure_mu(ThisLattice);
    const auto val_chi_u = val_mu * val_mu * ThisLattice.Get_ParaHamil().beta * ThisLattice.Get_NumSite();
    temp_MeasureData[SSE::INDEX_CHIU].AppendValue(val_chi_u);

    temp_NumOper[INDEX_Nt] = Num_Oper[OPERTYPE_t];
    temp_NumOper[INDEX_NJ] = Num_Oper[OPERTYPE_J];
    temp_NumOper[INDEX_NQ2] = Num_Oper[OPERTYPE_Q2];


    for (int index_Name = 0; index_Name != temp_MeasureData.size(); ++index_Name) {
        auto & which_MeaClass = temp_MeasureData[index_Name];
        if (which_MeaClass.Is_Empty()) {
            which_MeaClass.ClearValue();
            continue;
        }
        const auto temp_val = which_MeaClass.Get_AveValue();

        which_MeaClass.ClearValue();
        MeasureData[index_Name * (SSE::Vec_OperNum.size() + 1)].AppendValue(temp_val);

        for (int index_Oper = 0; index_Oper != SSE::Vec_OperNum.size(); ++index_Oper){
//        for (const auto &which_operkey: SSE::Vec_OperNum) {
            const auto temp_val_O = temp_val * temp_NumOper[index_Oper];
            MeasureData[index_Name * (SSE::Vec_OperNum.size() + 1) + 1 + index_Oper].AppendValue(temp_val_O);
        }
    }

    for (int index_Oper = 0; index_Oper != SSE::Vec_OperNum.size(); ++index_Oper){
        MeasureData[temp_MeasureData.size() * (SSE::Vec_OperNum.size() + 1) + index_Oper].AppendValue(temp_NumOper[index_Oper]);
    }
}

void SSE::Class_SSE::WriteBins() {
    // Outfile
    std::ofstream outfile;
    outfile.open(Str_FileData + Str_NameClass + Str_NameFile + ".txt", std::ofstream::out | std::ofstream::app);
//    outfile.setf(std::ios::fixed);
//    outfile.precision(16);

    // For all the quantities
    for (auto & which_MeaClass: MeasureData){
        auto temp_val = which_MeaClass.Get_AveValue();

        outfile << temp_val << "\t";
        which_MeaClass.ClearValue();
    }
    outfile << std::endl;
    outfile.close();
}

// Config file layout, all integers fixed-width in native byte order (little-endian on x86 and ARM):
// magic, int64 sweep, int64 Num_Site, int64 Num_Segment, int64 segment length,
// int32 spins, then (int32 index, int32 status) for every slot of every segment.
namespace {
    const char CONFIG_MAGIC[8] = {'S', 'S', 'E', 'Q', 'M', 'C', '0', '1'};

    void Write_Int64(std::ofstream &_which_file, std::int64_t _value) {
        _which_file.write(reinterpret_cast<const char *>(&_value), sizeof(_value));
    }

    bool Read_Int64(std::ifstream &_which_file, std::int64_t &_value) {
        _which_file.read(reinterpret_cast<char *>(&_value), sizeof(_value));
        return static_cast<bool>(_which_file);
    }
}

void SSE::Class_SSE::Write_Config(SSE::type_DataInt which_sweep) {
    // Write to a temporary file and rename it over the old one, so a job killed while writing keeps the last good config
    const auto str_filename = Str_FileConfig + Str_NameClass + Str_NameFile + "_config.bin";
    const auto str_filename_temp = str_filename + ".tmp";

    std::ofstream file_data(str_filename_temp, std::ios::out | std::ios::binary | std::ios::trunc);
    file_data.write(CONFIG_MAGIC, sizeof(CONFIG_MAGIC));
    Write_Int64(file_data, which_sweep);
    Write_Int64(file_data, ThisLattice.Get_NumSite());
    Write_Int64(file_data, Num_Segment);
    Write_Int64(file_data, Array_Oper.back().size());

    // Write Spin
    for (const auto &which_spin: ThisSpace.Get_Array_Spin()) {
        const std::int32_t temp_spin = static_cast<std::int32_t>(which_spin);
        file_data.write(reinterpret_cast<const char *>(&temp_spin), sizeof(temp_spin));
    }

    // Write the operator configuration
    for (const auto &which_Segment: Array_Oper) {
        for (const auto &which_oper: which_Segment) {
            which_oper.Write_Oper(file_data);
        }
    }
    file_data.close();

    if (!file_data) {
        std::cerr << "Warning! Could not write the config file " << str_filename_temp << std::endl;
        std::remove(str_filename_temp.c_str());
        return;
    }
    if (std::rename(str_filename_temp.c_str(), str_filename.c_str()) != 0) {
        std::cerr << "Warning! Could not rename " << str_filename_temp << " to " << str_filename << std::endl;
    }
}

SSE::type_DataInt SSE::Class_SSE::Read_Config(SSE::type_DataInt _total_sweep) {
    const auto str_filename = Str_FileConfig + Str_NameClass + Str_NameFile + "_config.bin";
    std::ifstream file_data(str_filename, std::ios::in | std::ios::binary);

    if (!file_data.good()) {
        std::cout << "There is no warmup file!" << std::endl;
        if (std::ifstream(Str_FileConfig + Str_NameClass + Str_NameFile + "_info.txt").good()) {
            std::cout << "Found a config in the old _info.txt/_data.bin format, which is no longer read." << std::endl;
        }
        return 0;
    }

    const auto reject_config = [&str_filename](const std::string &_reason) -> SSE::type_DataInt {
        std::cerr << "Warning! Ignoring the config file " << str_filename << ": " << _reason
                  << ". Starting from scratch." << std::endl;
        return 0;
    };

    // Header
    char temp_magic[sizeof(CONFIG_MAGIC)];
    file_data.read(temp_magic, sizeof(temp_magic));
    if (!file_data || !std::equal(std::begin(temp_magic), std::end(temp_magic), std::begin(CONFIG_MAGIC))) {
        return reject_config("not a config file of this version");
    }

    std::int64_t which_sweep = 0;
    std::int64_t num_site = 0;
    std::int64_t num_segment = 0;
    std::int64_t size_segment = 0;
    if (!Read_Int64(file_data, which_sweep) || !Read_Int64(file_data, num_site)
        || !Read_Int64(file_data, num_segment) || !Read_Int64(file_data, size_segment)) {
        return reject_config("truncated header");
    }
    if ((which_sweep < 0) || (num_site != ThisLattice.Get_NumSite()) || (num_segment != Num_Segment)
        || (size_segment <= 0)) {
        return reject_config("header does not match the lattice size or number of segments");
    }

    // Read into temporaries, and only replace the current configuration once the whole file checks out
    std::vector<SSE::type_DataInt> temp_Array_Spin(num_site);
    for (auto &which_spin: temp_Array_Spin) {
        std::int32_t temp_spin = 0;
        file_data.read(reinterpret_cast<char *>(&temp_spin), sizeof(temp_spin));
        if (!file_data || ((temp_spin != 1) && (temp_spin != -1))) {
            return reject_config("truncated or invalid spins");
        }
        which_spin = temp_spin;
    }

    std::vector<std::vector<SSE::Class_Oper>> temp_Array_Oper(Num_Segment, std::vector<SSE::Class_Oper>(size_segment));
    std::vector<SSE::type_DataInt> temp_Array_NumOper(Num_Segment, 0);
    for (type_DataInt index_Segment = 0; index_Segment != Num_Segment; ++index_Segment) {
        for (auto &which_oper: temp_Array_Oper[index_Segment]) {
            if (!which_oper.Read_Oper(file_data, ThisLattice)) {
                return reject_config("truncated or invalid operators");
            }
            if (!which_oper.If_Identity()) {
                ++temp_Array_NumOper[index_Segment];
            }
        }
    }
    if (file_data.peek() != std::ifstream::traits_type::eof()) {
        return reject_config("unexpected data at the end of the file");
    }

    // The operators must act on allowed spin states, and propagating through the whole string must come back
    // to the starting state
    SSE::Class_Space temp_Space = ThisSpace;
    temp_Space.Get_Array_Spin() = temp_Array_Spin;
    for (const auto &which_Segment: temp_Array_Oper) {
        for (const auto &which_oper: which_Segment) {
            if (which_oper.If_Identity()) {
                continue;
            }
            auto temp_oper = which_oper;
            if (!temp_oper.If_Suitable_andUpdateTypeOper(temp_Space, ThisLattice)) {
                return reject_config("operator acts on a forbidden spin state");
            }
            which_oper.If_Dia_IfnotUpdateSpace(temp_Space, ThisLattice);
        }
    }
    if (temp_Space.Get_Array_Spin() != temp_Array_Spin) {
        return reject_config("operator string is not periodic in imaginary time");
    }

    ThisSpace.Get_Array_Spin() = temp_Array_Spin;
    Array_Oper = std::move(temp_Array_Oper);
    Array_NumOper = std::move(temp_Array_NumOper);

    return static_cast<SSE::type_DataInt>(which_sweep);
}

void SSE::Class_SSE::Measure_Corrf() {
    // Initial spin configuration vector

    static SSE::Class_fMat temp_fMat(ThisLattice.Get_NumSite());
    // Note it taks the value of 0 or 1
    temp_fMat.Reset( ThisSpace.Get_Array_Spin());

    static std::vector<std::valarray<type_DataFloat>> array_vec_Gq(ThisLattice.Get_NumSite(), std::valarray<type_DataFloat>(Num_Measure_TimeSlices));

    SSE::type_DataInt val_taudiff = 0;
    bool if_Measure = true;
    for (type_DataInt index_Segment = 0; index_Segment != Num_Segment; ++index_Segment){
        // Too few segments left for a full set of time slices: stop measuring, but keep propagating the space
        // through the rest of the operator string, so it ends up back at tau = 0
        if ((Num_Measure_TimeSlices - val_taudiff) > (Num_Segment - index_Segment)) {
            if_Measure = false;
        }

        // Update ftau
        for (auto &which_oper : Array_Oper[index_Segment]) {
            if (!which_oper.If_Identity()) {
                // Just use this function to update the space
                which_oper.If_Dia_IfnotUpdateSpace(ThisSpace, ThisLattice);

                // Update fMat only if it is not the last slice
                if (if_Measure) {
                    which_oper.Update_fMat(temp_fMat, ThisLattice);
                }
            }
        }
        if (!if_Measure) {
            continue;
        }

        // value of tau plus one
        ++val_taudiff;

        // Measure
//        std::vector<type_DataFloat> temp_Gtau(ThisLattice.Num_Site, 0.);
        temp_fMat.Measure( Corr_ftau[val_taudiff - 1], ThisLattice);

        // If it measured the full time slices
        if (val_taudiff == Num_Measure_TimeSlices) {
            // Restart to count the time slices
            temp_fMat.Reset(ThisSpace.Get_Array_Spin());

            val_taudiff = 0;
        }
    }
}

void SSE::Class_SSE::Measure_Corrf_Simple() {
// Initial spin configuration vector
    // Note it taks the value of 0 or 1
    auto temp_Gx = std::valarray<type_DataFloat>(0., ThisLattice.Get_NumSite());
    temp_Gx[0] = 0.5;

//    static std::vector<std::valarray<type_DataFloat>> array_vec_Gq(ThisLattice.Get_NumSite(), std::valarray<type_DataFloat>(Num_Measure_TimeSlices));

    SSE::type_DataInt val_taudiff = 0;
    bool if_Measure = true;
    for (type_DataInt index_Segment = 0; index_Segment != Num_Segment; ++index_Segment){
        // Same as Measure_Corrf: keep propagating the space back to tau = 0 after measuring stops
        if ((Num_Measure_TimeSlices - val_taudiff) > (Num_Segment - index_Segment)) {
            if_Measure = false;
        }

        // Update ftau
        for (auto &which_oper : Array_Oper[index_Segment]) {
            if (!which_oper.If_Identity()) {
                // Just use this function to update the space
                which_oper.If_Dia_IfnotUpdateSpace(ThisSpace, ThisLattice);

                // Update fMat only if it is not the last slice
                if (if_Measure) {
                    which_oper.Update_fVec(temp_Gx, ThisLattice);
                }
            }
        }
        if (!if_Measure) {
            continue;
        }

        // value of tau plus one
        ++val_taudiff;

        Corr_ftau[val_taudiff - 1].AppendValue(temp_Gx);

        // Measure
        if (val_taudiff == Num_Measure_TimeSlices) {
            temp_Gx = 0.;
            temp_Gx[0] = 0.5;

            val_taudiff = 0;
        }
    }
}


void SSE::Class_SSE::Write_Corrf() {
    std::ofstream outfile;
    outfile.open(Str_FileCorr + Str_NameClass + Str_NameFile + ".txt", std::ofstream::out | std::ofstream::app);
//    outfile.setf(std::ios::fixed);
//    outfile.precision(10);

    // For all q-values
    for (type_DataInt index_tau = 0; index_tau != Num_Measure_TimeSlices; ++index_tau) {
        auto temp_Ave = Corr_ftau[index_tau].Get_AveValue();
        Corr_ftau[index_tau].ClearValue();
        for (type_NumSite index_site = 0; index_site != ThisLattice.Get_NumSite(); ++index_site){
            SSE::type_DataFloat val_tau = static_cast<SSE::type_DataFloat>(index_tau + 1) * ThisLattice.Get_ParaHamil().beta / this->Array_NumOper.size();
            outfile << index_site << '\t' << val_tau << '\t' << temp_Ave[index_site] << std::endl;
        }
    }
    outfile.close();
}

