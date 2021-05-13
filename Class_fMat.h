//
// Created by Bowen Zhao (bwzhao@bu.edu) on 7/17/20.
//
#pragma once

#include <vector>
#include "Config.h"
#include <valarray>
#include <Eigen/Eigen>
#include "Class_CorrMeasurement.h"

namespace SSE {
    class Class_fMat {
    private:
        Eigen::Matrix<type_DataFloat , Eigen::Dynamic, Eigen::Dynamic> Mat;
        std::vector<type_DataInt> Array_Index;
        type_NumSite MatLinSize;
    public:
        // Constructor & destructor
        // Default
        Class_fMat() = default;
        // Initialize it as a identity matrix
        explicit Class_fMat(SSE::type_NumSite _size):
                MatLinSize(_size),
                Array_Index(_size / 2),
                Mat(Eigen::Matrix<type_DataFloat, Eigen::Dynamic, Eigen::Dynamic>::Zero(_size, _size / 2))
        {}
        ~Class_fMat() = default;

        // Get element
        inline type_DataFloat & Get(type_NumSite _index_row, type_NumSite _index_col){
            return this->Mat(_index_row, _index_col);
        }

//        inline void Set_Row(type_NumSite _index_row, T _val){
//            this->Mat.row(_index_row) = _val;
//        }
        inline void Set_RowZero(type_NumSite _index_row){
            Mat.row(_index_row).setZero();
        }

        inline void Multi_Row(type_NumSite _index_row, type_DataFloat _val){
            this->Mat.row(_index_row) *= _val;
        }

        inline void Swap_Row(type_NumSite _index_row1, type_NumSite _index_row2){
            Mat.row(_index_row1).swap(Mat.row(_index_row2));
        }

        inline void SwapMinus_Row(type_NumSite _index_row1, type_NumSite _index_row2){
            const auto temp_row_1 = Mat.row(_index_row1);
            Mat.row(_index_row1) -= Mat.row(_index_row2);
            Mat.row(_index_row2) -= temp_row_1;
        }

        inline type_NumSite Get_Size(){
            return MatLinSize;
        }
//        const T& Get(type_NumSite _index_row, type_NumSite _index_col) const{
//            return Mat[_index_row + _index_col * this->MatLinSize];
//        }

        // Reset element
        inline void Reset(const std::vector<type_DataInt> & _arraySpin){
            type_DataInt index_col = 0;
            for (type_DataInt index_spin = 0; index_spin != _arraySpin.size(); ++index_spin) {
                Array_Index[index_col] = index_spin;
                Mat(index_spin, index_col) = 0.5;
            }
        }

        void Measure(Class_CorrMeasurement <type_DataFloat> &_which_CorrMeaa, const SSE::Class_Lattice &_which_Lattice);
    };

    inline void Class_fMat::Measure(Class_CorrMeasurement <type_DataFloat> &_which_CorrMeaa,
                                    const SSE::Class_Lattice &_which_Lattice) {
        auto temp_arrayReturn = std::valarray<type_DataFloat>(0., MatLinSize);

        for (type_DataInt index_col = 0; index_col != Array_Index.size(); ++index_col) {
            const auto & temp_col = Mat.col(index_col);
            const auto &temp_arrayDiff = _which_Lattice.Get_array_Sitediff(Array_Index[index_col]);
            for (type_DataInt index_diff = 0; index_diff != MatLinSize; ++index_diff) {
                temp_arrayReturn[index_diff] += temp_col[temp_arrayDiff[index_diff]];
            }
        }

        _which_CorrMeaa.AppendValue(temp_arrayReturn / MatLinSize);
    }

};

//class Class_fMat {
//private:
//    std::vector<std::valarray<T>> Mat;
//    std::vector<std::valarray<T>> Mat_Id;
//    type_NumSite MatLinSize;
//public:
//    // Constructor & destructor
//    // Default
//    Class_fMat() = default;
//    // Initialize it as a identity matrix
//    explicit Class_fMat(SSE::type_NumSite _size):
//            Mat(_size, std::valarray<T>(_size)),
//            MatLinSize(_size)
//    {
//        for (type_NumSite index_1 = 0; index_1 != MatLinSize; ++index_1) {
//            for (type_NumSite index_2 = 0; index_2 != MatLinSize; ++index_2) {
//                if (index_1 == index_2) {
//                    Mat[index_1][index_2] = 1;
//                }
//                else{
//                    Mat[index_1][index_2] = 0;
//                }
//            }
//        }
//        Mat_Id = Mat;
//    }
//    ~Class_fMat() = default;
//
//    // Get element
//    inline T& Get(type_NumSite _index_row, type_NumSite _index_col){
//        return this->Mat[_index_row][_index_col];
//    }
//
//    inline void Set(type_NumSite _index_row, type_NumSite _index_col, T _val){
//        this->Mat[_index_row][_index_col] = _val;
//    }
//
//    inline void Set_Row(type_NumSite _index_row, T _val){
//        this->Mat[_index_row] = _val;
//    }
//
//    inline void Multi_Row(type_NumSite _index_row, T _val){
//        this->Mat[_index_row] *= _val;
//    }
//
//    inline void Swap_Row(type_NumSite _index_row1, type_NumSite _index_row2){
//        Mat[_index_row1].swap(Mat[_index_row2]);
//    }
//
//    inline void SwapMinus_Row(type_NumSite _index_row1, type_NumSite _index_row2){
//        const auto temp_row_1 = Mat[_index_row1];
//        Mat[_index_row1] -= Mat[_index_row2];
//        Mat[_index_row2] -= temp_row_1;
//    }
//
//    inline type_NumSite Get_Size(){
//        return MatLinSize;
//    }
////        const T& Get(type_NumSite _index_row, type_NumSite _index_col) const{
////            return Mat[_index_row + _index_col * this->MatLinSize];
////        }
//
//    // Reset element
//    inline void Reset(){
//        Mat = Mat_Id;
//    }
//};