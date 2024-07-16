//
// Created by cepheid on 7/16/24.
//

#ifndef CURLS_H
#define CURLS_H

#include "enums.h"

//====== Curl Operators =======
//=============================
 template<Derivative D, FieldType FT, typename... IDXS>
 struct curl {
   static constexpr auto apply(const auto&, IDXS...) { return 0.0; }
 };

 template<FieldType FT, typename... IDXS>
 struct curl<Derivative::DX, FT, IDXS...> {
   static auto apply(const auto& f, IDXS... idxs) {
     if constexpr (FT == FieldType::E) {
       return f.backward_diff_x(idxs...);
     } else {
       return f.forward_diff_x(idxs...);
     }
   }
 };

 template<FieldType FT, typename... IDXS>
 struct curl<Derivative::DY, FT, IDXS...> {
   static auto apply(const auto& f, IDXS... idxs) {
     if constexpr (FT == FieldType::E) {
       return f.backward_diff_y(idxs...);
     } else {
       return f.forward_diff_y(idxs...);
     }
   }
 };

 template<FieldType FT, typename... IDXS>
 struct curl<Derivative::DZ, FT, IDXS...> {
   static auto apply(const auto& f, IDXS... idxs) {
     if constexpr (FT == FieldType::E) {
       return f.backward_diff_z(idxs...);
     } else {
       return f.forward_diff_z(idxs...);
     }
   }
 };


#endif //CURLS_H
