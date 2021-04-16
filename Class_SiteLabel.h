////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author:          Bowen Zhao <bwzhao@bu.edu>
// Organization:    Boston University
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#pragma once

#include "Config.h"
#include <vector>

namespace SSE{
    class Class_SiteLabel {
    private:
        std::vector<SSE::type_DataFloat> Array_Label;
        type_DataInt Num_Dimension;
    public:
        Class_SiteLabel() = default;
        // 1D case
        Class_SiteLabel(SSE::type_DataFloat _label_0){
            Array_Label.emplace_back(_label_0);
            Num_Dimension = 1;
        }
        // 2D case
        Class_SiteLabel(SSE::type_DataFloat _label_0, SSE::type_DataFloat _label_1){
            Array_Label.emplace_back(_label_0);
            Array_Label.emplace_back(_label_1);

            Num_Dimension = 2;
        }
        // general case
        Class_SiteLabel(const std::vector<SSE::type_DataFloat> & _array_Label){
            Array_Label = _array_Label;
            Num_Dimension = _array_Label.size();
        }

        SSE::type_DataFloat operator*(const Class_SiteLabel &_other_SiteLabel) const {
            SSE::type_DataFloat temp_return = 0.;
            for (decltype(this->Num_Dimension) index_dimension = 0;
                 index_dimension != Num_Dimension; ++index_dimension) {
                temp_return += this->Array_Label[index_dimension] * _other_SiteLabel.Array_Label[index_dimension];
            }
            return temp_return;
        }

        decltype(SSE::Class_SiteLabel::Num_Dimension) Get_Dim() const{
            return this->Num_Dimension;
        }

        SSE::type_DataFloat Get_Component(const int & _index_component) const {
            return this->Array_Label[_index_component];
        }
    };
}