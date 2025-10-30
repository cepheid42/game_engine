#ifndef EM_DEFINITIONS_HPP
#define EM_DEFINITIONS_HPP

#include "program_params.hpp"
#include "typelist.hpp"
#include "bc_data.hpp"
#include "bc_functors.hpp"
#include "update_functors.hpp"

namespace tf::electromagnetics {

template<EMFace F, EMSide S>
using BCDataTypes = TypeList<
   PeriodicFaceBC<F, S>,   // 0
   PMLFaceBC<F, S>,        // 1
   ReflectingFaceBC<F, S>  // 2
>;

using bcdata_t = BCData<
   TypeListAt<BCSelect[0], BCDataTypes<EMFace::X, EMSide::Lo>>, // x0
   TypeListAt<BCSelect[1], BCDataTypes<EMFace::X, EMSide::Hi>>, // x1
   TypeListAt<BCSelect[2], BCDataTypes<EMFace::Y, EMSide::Lo>>, // y0
   TypeListAt<BCSelect[3], BCDataTypes<EMFace::Y, EMSide::Hi>>, // y1
   TypeListAt<BCSelect[4], BCDataTypes<EMFace::Z, EMSide::Lo>>, // z0
   TypeListAt<BCSelect[5], BCDataTypes<EMFace::Z, EMSide::Hi>>  // z1
>;

template<EMFace F, EMSide S>
using BCFuncTypes = TypeList<
   PeriodicBoundary<F, S>,  // 0
   PMLBoundary<F, S>,       // 1
   ReflectingBoundary<F, S>       // 2
>;

using EMDiffTypes = TypeList<
   noop,                   // 0
   forward_dx,             // 1
   backward_dx,            // 2
   forward_dy,             // 3
   backward_dy,            // 4
   forward_dz,             // 5
   backward_dz             // 6
>;

using boundary_t = TypeList<
   TypeListAt<BCSelect[0], BCFuncTypes<EMFace::X, EMSide::Lo>>, // x0
   TypeListAt<BCSelect[1], BCFuncTypes<EMFace::X, EMSide::Hi>>, // x1
   TypeListAt<BCSelect[2], BCFuncTypes<EMFace::Y, EMSide::Lo>>, // y0
   TypeListAt<BCSelect[3], BCFuncTypes<EMFace::Y, EMSide::Hi>>, // y1
   TypeListAt<BCSelect[4], BCFuncTypes<EMFace::Z, EMSide::Lo>>, // z0
   TypeListAt<BCSelect[5], BCFuncTypes<EMFace::Z, EMSide::Hi>>  // z1
>;

// turns differences on and off depending on simulation geometry
// note that if x_collapsed, y_collapsed, and z_collapsed are all set (0D), the field solver won't do anything!
using field_t = TypeList<
   FieldIntegrator<ExplicitUpdateFunctor<TypeListAt<y_collapsed ? 0 : 4, EMDiffTypes>, TypeListAt<z_collapsed ? 0 : 6, EMDiffTypes>>>,
   FieldIntegrator<ExplicitUpdateFunctor<TypeListAt<z_collapsed ? 0 : 6, EMDiffTypes>, TypeListAt<x_collapsed ? 0 : 2, EMDiffTypes>>>,
   FieldIntegrator<ExplicitUpdateFunctor<TypeListAt<x_collapsed ? 0 : 2, EMDiffTypes>, TypeListAt<y_collapsed ? 0 : 4, EMDiffTypes>>>,
   FieldIntegrator<ExplicitUpdateFunctor<TypeListAt<z_collapsed ? 0 : 5, EMDiffTypes>, TypeListAt<y_collapsed ? 0 : 3, EMDiffTypes>>>,
   FieldIntegrator<ExplicitUpdateFunctor<TypeListAt<x_collapsed ? 0 : 1, EMDiffTypes>, TypeListAt<z_collapsed ? 0 : 5, EMDiffTypes>>>,
   FieldIntegrator<ExplicitUpdateFunctor<TypeListAt<y_collapsed ? 0 : 3, EMDiffTypes>, TypeListAt<x_collapsed ? 0 : 1, EMDiffTypes>>>
>;

} // end namespace tf::electromagnetics
#endif //EM_DEFINITIONS_HPP