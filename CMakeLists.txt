cmake_minimum_required(VERSION 3.17)
project(20210327_NewProgram)

set(CMAKE_CXX_STANDARD 11)
#find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(20210327_NewProgram
        main.cpp
        Config.h
        Class_SSE.cpp
        Class_SSE.h
        Class_Lattice.cpp
        Class_Lattice.h
        Class_CorrMeasurement.h
        Class_DataMeasurement.h
        Class_Oper.cpp
        Class_Oper.h
        Class_Space.cpp
        Class_Space.h
        Class_LegLink.h
        Class_SiteLabel.h
        Class_fMat.h)

#target_link_libraries (20210327_NewProgram Eigen3::Eigen)