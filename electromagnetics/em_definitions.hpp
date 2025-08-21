#ifndef EM_DEFINITIONS_HPP
#define EM_DEFINITIONS_HPP

#include "em_params.hpp"
#include "typelist.hpp"
#include "bc_data.hpp"
#include "bc_functors.hpp"
#include "update_functors.hpp"

namespace tf::electromagnetics {

template<EMFace F, EMSide S>
using BCDataTypes = TypeList<
   ReflectingFaceBC<F, S>,  // 0
   PeriodicFaceBC<F, S>,    // 1
   PMLFaceBC<F, S>          // 2
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
   ReflectingBoundary,      // 0
   PeriodicBoundary<F, S>,  // 1
   PMLBoundary<F, S>        // 2
>;

using boundary_t = TypeList<
   TypeListAt<BCSelect[0], BCFuncTypes<EMFace::X, EMSide::Lo>>, // x0
   TypeListAt<BCSelect[1], BCFuncTypes<EMFace::X, EMSide::Hi>>, // x1
   TypeListAt<BCSelect[2], BCFuncTypes<EMFace::Y, EMSide::Lo>>, // y0
   TypeListAt<BCSelect[3], BCFuncTypes<EMFace::Y, EMSide::Hi>>, // y1
   TypeListAt<BCSelect[4], BCFuncTypes<EMFace::Z, EMSide::Lo>>, // z0
   TypeListAt<BCSelect[5], BCFuncTypes<EMFace::Z, EMSide::Hi>>  // z1
>;

// // 2D TEy
// using field_t = TypeList<
//    FieldIntegrator<ExplicitUpdateFunctor<noop, backward_dz>>,
//    FieldIntegrator<ExplicitUpdateFunctor<backward_dz, backward_dx>>,
//    FieldIntegrator<ExplicitUpdateFunctor<backward_dx, noop>>,
//    FieldIntegrator<ExplicitUpdateFunctor<forward_dz, noop>>,
//    FieldIntegrator<ExplicitUpdateFunctor<forward_dx, forward_dz>>,
//    FieldIntegrator<ExplicitUpdateFunctor<noop , forward_dx>>
// >;

// 3D
using field_t = TypeList<
   FieldIntegrator<ExplicitUpdateFunctor<backward_dy, backward_dz>>,
   FieldIntegrator<ExplicitUpdateFunctor<backward_dz, backward_dx>>,
   FieldIntegrator<ExplicitUpdateFunctor<backward_dx, backward_dy>>,
   FieldIntegrator<ExplicitUpdateFunctor<forward_dz, forward_dy>>,
   FieldIntegrator<ExplicitUpdateFunctor<forward_dx, forward_dz>>,
   FieldIntegrator<ExplicitUpdateFunctor<forward_dy, forward_dx>>
>;

} // end namespace tf::electromagnetics
#endif //EM_DEFINITIONS_HPP