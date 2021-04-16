////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author:          Bowen Zhao <bwzhao@bu.edu>
// Organization:    Boston University
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#pragma once

#include <map>
#include <string>
#include <vector>
#include "Config.h"
#include <numeric>

namespace SSE {
    template <typename T>
    class Class_DataMeasurement {
    private:
        T Data;
        SSE::type_DataInt Length;

    public:
        Class_DataMeasurement():Length(0){};
        ~Class_DataMeasurement() = default;

        // Add a measurement value to a specific quantity(key)
        void AppendValue(const T _which_value){
            if (Length == 0) {
                Data = _which_value;
            }
            else{
                Data += _which_value;
            }
            ++Length;
        }

        T Get_AveValue(){
            return (Data / Length);
        }

        void ClearValue(){
            Data = 0;
            Length = 0;
        }

        bool Is_Empty(){
            return (Length == 0);
        }

        SSE::type_DataInt Size() const{
            return Length;
        }

    };
}

//namespace SSE {
//    template <typename T>
//    class Class_DataMeasurement {
//    private:
//        std::vector<T> Data;
//
//    public:
//        Class_DataMeasurement() = default;
//        ~Class_DataMeasurement() = default;
//
//        // Add a measurement value to a specific quantity(key)
//        void AppendValue(const T _which_value){
//            Data.emplace_back(_which_value);
//        }
//        T Get_AveValue(){
//            return std::accumulate(Data.begin(), Data.end(), 0.) / Data.size();
//        }
//
//        void ClearValue(){
//            Data.clear();
//        }
//
//    };
//}



