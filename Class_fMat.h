//
// Created by Bowen Zhao (bwzhao@bu.edu) on 7/17/20.
//
#pragma once

#include <vector>
#include "Config.h"
#include <valarray>

namespace SSE {
    template <typename T>
    class Class_fMat {
    private:
        std::vector<std::valarray<T>> Mat;
        std::vector<std::valarray<T>> Mat_Id;
        type_NumSite MatLinSize;
    public:
        // Constructor & destructor
        // Default
        Class_fMat() = default;
        // Initialize it as a identity matrix
        explicit Class_fMat(SSE::type_NumSite _size):
                Mat(_size, std::valarray<T>(_size)),
                MatLinSize(_size)
        {
            for (type_NumSite index_1 = 0; index_1 != MatLinSize; ++index_1) {
                for (type_NumSite index_2 = 0; index_2 != MatLinSize; ++index_2) {
                    if (index_1 == index_2) {
                        Mat[index_1][index_2] = 1;
                    }
                    else{
                        Mat[index_1][index_2] = 0;
                    }
                }
            }
            Mat_Id = Mat;
        }
        ~Class_fMat() = default;

        // Get element
        inline T& Get(type_NumSite _index_row, type_NumSite _index_col){
            return this->Mat[_index_row][_index_col];
        }

        inline void Set(type_NumSite _index_row, type_NumSite _index_col, T _val){
            this->Mat[_index_row][_index_col] = _val;
        }

        inline void Set_Row(type_NumSite _index_row, T _val){
            this->Mat[_index_row] = _val;
        }

        inline void Multi_Row(type_NumSite _index_row, T _val){
            this->Mat[_index_row] *= _val;
        }

        inline void Swap_Row(type_NumSite _index_row1, type_NumSite _index_row2){
            Mat[_index_row1].swap(Mat[_index_row2]);
        }

        inline void SwapMinus_Row(type_NumSite _index_row1, type_NumSite _index_row2){
            const auto temp_row_1 = Mat[_index_row1];
            Mat[_index_row1] -= Mat[_index_row2];
            Mat[_index_row2] -= temp_row_1;
        }

        inline type_NumSite Get_Size(){
            return MatLinSize;
        }
//        const T& Get(type_NumSite _index_row, type_NumSite _index_col) const{
//            return Mat[_index_row + _index_col * this->MatLinSize];
//        }

        // Reset element
        inline void Reset(){
            Mat = Mat_Id;
        }
    };
};