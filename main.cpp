////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author:          Bowen Zhao <bwzhao@bu.edu>
// Organization:    Boston University
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <vector>
#include "Config.h"
#include "Class_SSE.h"

int main(int argc, char *argv[]) {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Pre-calculation
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //input parameters
    std::vector<std::string> para_main;
    for (decltype(argc) nl_arg = 0; nl_arg != argc; ++nl_arg) {
        std::string temp_string(argv[nl_arg]);
        para_main.push_back(temp_string);
    }
    // Class_Lattice Properties
    // ParaHamiltonian:
    const int index_StartPara = 11;
    SSE::type_ParaHamil ParaHamil = {
            std::stod(para_main[index_StartPara + 0]), //ty
            std::stod(para_main[index_StartPara + 1]), //tx
            std::stod(para_main[index_StartPara + 2]), //Jy
            std::stod(para_main[index_StartPara + 3]), //Jx
            std::stod(para_main[index_StartPara + 4]), //Q2
            std::stod(para_main[index_StartPara + 5]), //Q2s
            std::stod(para_main[index_StartPara + 6]), //JfA
            std::stod(para_main[index_StartPara + 7]), //JfB
            std::stod(para_main[index_StartPara + 8]), //Q3
            std::stod(para_main[index_StartPara + 9]), //Q3s

            std::stoi(para_main[index_StartPara + 10]), //Num_X
            std::stoi(para_main[index_StartPara + 11]), //Num_Y
            std::stod(para_main[index_StartPara + 12]) //beta
    };
    // Which_BC: boundary conditions
    SSE::type_DataInt Which_BC = std::stoi(para_main[1]);

    // Num_Warmup: minimum number sweeps of warmup
    // Num_SweepinBin: number sweeps in each measuring "bin"
    // Num_Bin: number of bin to measure this time
    SSE::type_DataInt Num_Warmup = std::stoi(para_main[2]);
    SSE::type_DataInt Num_SweepinBin = std::stoi(para_main[3]);
    SSE::type_DataInt Num_Bin = std::stoi(para_main[4]);

    // num_cpu:
    // value = 0: just for warm up,
    // value > 0: for i^th measurement
    SSE::type_DataInt Which_Cpu = std::stoi(para_main[5]);
    SSE::type_DataInt Num_Segment = std::stoi(para_main[6]);
    SSE::type_DataInt Num_MeasureSegment = std::stoi(para_main[7]);

    std::string file_config = para_main[8];
    std::string file_data = para_main[9];
    std::string file_corr = para_main[10];

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Start to Calculate
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Generate the Calculate class
    SSE::Class_SSE LittleLion(ParaHamil,
                              Which_Cpu,
                              Which_BC,
                              file_config,
                              file_data,
                              file_corr,
                              Num_Segment,
                              Num_MeasureSegment
    );
//    // Read the warm-up file, if any. Or current_sweep = 0;
    SSE::type_DataInt current_sweep = 0;
    current_sweep = LittleLion.Read_Config(Num_Warmup);
//
//    // Start to do warm-up
//    // The system must have enough warm-up steps to start doing measurement.
    if (current_sweep < Num_Warmup) {
        std::cout << "Still need to Warmup" << std::endl;
        LittleLion.Print_Time();
        for (auto index_sweep = current_sweep; index_sweep != Num_Warmup + 1; ++index_sweep) {
            LittleLion.DiagonalUpdate();
            LittleLion.Adjust_Cutoff();
            LittleLion.Generate_Linklist();
//            LittleLion.CheckLinkList();
            LittleLion.Flip_Update();
            LittleLion.Loop_Update();
            // Store the configurations.
            if (Num_Warmup / 20 != 0) {
                if (index_sweep / (Num_Warmup / 20) * (Num_Warmup / 20) == index_sweep) {
                    LittleLion.Write_Config(index_sweep);
                }
            }
        }
    }
//
//    // Start to meausre
    std::cout << "Finished Warmup" << std::endl;
    LittleLion.Print_Time();

    for (int index_bin = 0; index_bin != Num_Bin; ++index_bin) {
        for (int nl_sweep = 0; nl_sweep != Num_SweepinBin; nl_sweep++) {
            LittleLion.DiagonalUpdate();
            LittleLion.Adjust_Cutoff();
            LittleLion.Generate_Linklist();
            LittleLion.Flip_Update();
            LittleLion.Loop_Update();
//            LittleLion.Measure_woDynamics();
            LittleLion.Measure_Corrf();
        }
//        LittleLion.WriteBins();
        LittleLion.Write_Corrf();

        // Store the configurations
        if (Num_Bin / 10 != 0) {
            if (index_bin / (Num_Bin / 10) * (Num_Bin / 10) == index_bin) {
                std::cout << "#bin #" << index_bin << "/" << Num_Bin << " complete! " << std::endl;
                LittleLion.Print_Time();
//                LittleLion.Write_Config(
//                        (index_bin + 1) * Num_SweepinBin + current_sweep + std::max(0, int(Num_Warmup - current_sweep)));
            }
        }
    }
    LittleLion.Print_Time();
}
