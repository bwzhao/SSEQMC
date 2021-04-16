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
#include <cmath>
#include <complex>
#include <algorithm>
#include <valarray>

namespace SSE {
    template <typename T>
    class Class_CorrMeasurement {
    private:
        std::valarray<T> Data;
        SSE::type_DataInt Length;

    public:
        Class_CorrMeasurement():Length(0){};
//        Class_CorrMeasurement(SSE::type_NumSegment _length):
//                Length(_length) {}
        ~Class_CorrMeasurement() = default;

        // Add a measurement value to a specific quantity(key)
//        void AppendValue(const std::vector<T> &_which_value){
//            if (Length == 0) {
//                Data = _which_value;
//            }
//            else{
//                for (SSE::type_NumSegment index_r = 0; index_r != Data.size(); ++index_r) {
//                    Data[index_r] += _which_value[index_r];
//                }
//            }
//            ++Length;
//        }

        void AppendValue(const std::valarray<T> &_which_value){
            if (Length == 0) {
                Data = _which_value;
            }
            else{
                Data += _which_value;
//                for (SSE::type_NumSegment index_r = 0; index_r != Data.size(); ++index_r) {
//                    Data[index_r] += _which_value[index_r];
//                }
            }
            ++Length;
        }

        const std::valarray<T> & Get_AveValue(){
//            std::vector<T> temp_return_vec(Data);
//            for (SSE::type_NumSegment index_r = 0; index_r != Data.size(); ++index_r) {
//                temp_return_vec[index_r] /= Length;
//            }
//            std::cout << Length << std::endl;
            Data /= Length;
            return Data;
        }

        void ClearValue(){
            Data = 0;
            Length = 0;
        }
    };
}

