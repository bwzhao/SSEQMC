////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author:          Bowen Zhao <bwzhao@bu.edu>
// Organization:    Boston University
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once
#include <iostream>
#include "Config.h"
#include "Class_Lattice.h"
#include "Class_Space.h"
#include "Class_Oper.h"
#include "Class_DataMeasurement.h"
#include "Class_CorrMeasurement.h"
#include "Class_LegLink.h"
#include <fstream>

namespace SSE{
    class Class_SSE {
    private:
        // Lattice
        const Class_Lattice ThisLattice;

        //Oper Array
        SSE::type_DataInt Num_Segment;
        std::vector<std::vector<SSE::Class_Oper>> Array_Oper;
        std::vector<std::vector<Class_Oper *>> Array_OperMapping;
        std::vector<std::vector<std::vector<class_LegLink>>> Array_OperLinkList;

        std::vector<SSE::type_DataInt> Array_NumOper;

        //Space, including spins
        SSE::Class_Space ThisSpace;
        std::vector<std::vector<class_LegLink>> Array_SpaceLinkList;

        //Stored MeasureData:
        std::vector<SSE::Class_DataMeasurement<SSE::type_DataFloat>> MeasureData;

        // Index: (index_site, index_slice)
        std::vector<SSE::Class_CorrMeasurement<SSE::type_DataFloat>> Corr_ftau;

        type_DataInt Num_Measure_TimeSlices;

        //Some stored names
        std::string Str_NameFile;
        const std::string Str_FileConfig;
        const std::string Str_FileData;
        const std::string Str_FileCorr;
        const std::string Str_NameClass;

    public:
        // Constructor & Destructor
        explicit Class_SSE(type_ParaHamil _para_Hamil,
                           type_DataInt _which_cpu,
                           type_DataInt _bc,
                           const std::string & _file_config,
                           const std::string & _file_data,
                           const std::string & _file_corr,
                           type_DataInt _num_Segment,
                           type_DataInt _num_TimeSlices
        );
        ~Class_SSE() = default;
        void Print_Time();

        // Update procedure
        void DiagonalUpdate();
        void Adjust_Cutoff();
        void Generate_Linklist();
        void Flip_Update();
        void Loop_Update();

        void Measure_woDynamics();
        void WriteBins();

        // Functions related to lilnklist
        bool If_CanFlip(type_NumSite _which_site);
        SSE::class_LegLink Find_Nextinleg_UpdateSpace(const class_LegLink &_which_LegLink, bool _ifFlip);
        void CheckLinkList() const;
    };
}


inline void SSE::Class_SSE::Print_Time() {
    std::time_t rawtime;
    struct tm *timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    std::cout << std::asctime(timeinfo) << std::endl;
}


inline bool SSE::Class_SSE::If_CanFlip(type_NumSite _which_site) {
    return Array_SpaceLinkList[_which_site][0].Which_Segment == -1;
}
