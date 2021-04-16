////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Author:          Bowen Zhao <bwzhao@bu.edu>
// Organization:    Boston University
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#pragma once

#include "Config.h"

namespace SSE{
    struct class_LegLink {
        SSE::type_DataInt Which_Segment;
        SSE::type_DataInt Which_smLeg ;
        SSE::type_DataInt Which_Time;

        explicit class_LegLink():
                Which_Segment(-1),
                Which_smLeg(-1),
                Which_Time(-1){}

        class_LegLink(SSE::type_DataInt _which_Segment, SSE::type_DataInt _which_Time, SSE::type_DataInt _which_Leg):
                Which_Segment(_which_Segment),
                Which_smLeg(_which_Leg),
                Which_Time(_which_Time){}

        ~class_LegLink() = default;

        bool If_Unlinked() const {
            return (Which_Segment == -1);
        }

        bool If_Unupdated() const {
            return (Which_Segment != -1);
        }

        bool If_OperLink() const{
            return (Which_Segment >= 0);
        }

        void Set_LegLink(SSE::type_DataInt _which_Segment, SSE::type_DataInt _which_Time, SSE::type_DataInt _which_Leg){
            this->Which_Segment = _which_Segment;
            this->Which_Time = _which_Time;
            this->Which_smLeg = _which_Leg;
        }

        void Reset_LegLink(){
            this->Which_Segment = -1;
            this->Which_Time = -4;
            this->Which_smLeg = -5;
        }

        void Set_LegLink(const SSE::class_LegLink & _which_LegLink){
            *this = _which_LegLink;
        }

        bool operator!=(const class_LegLink &_other_LegLink){
            return (Which_Segment != _other_LegLink.Which_Segment)
                   ||(Which_Time != _other_LegLink.Which_Time)
                   ||(Which_smLeg != _other_LegLink.Which_smLeg);
        }
    };

}
