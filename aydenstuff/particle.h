//
// Created by akis on 8/27/24.
//

#include <concepts>
#include <cmath>
#include <vector>
#include <algorithm>

#include "/home/akis/projects/triforce/libs/common/include/constants.h"
#include "control.h"
#include "tags.h"
#include "utilities.h"
#include "vector.h"

#ifndef TRIFORCE_PARTICLE_H
#define TRIFORCE_PARTICLE_H

// ===== Particle Data Concepts =====
// ==================================
namespace tf::concepts
{
  template <typename T>
  concept particle = requires {
      requires std::floating_point<typename T::value_t>;
      requires is_dimensional<typename T::location_t>;
      requires is_dimensional<typename T::momentum_t>;
      
      requires std::same_as<std::decay<decltype(T::location)>, typename T::location_t>;
      requires std::same_as<std::decay<decltype(T::old_location)>, typename T::location_t>;
      requires std::same_as<std::decay<decltype(T::momentum)>, typename T::momentum_t>;
      requires std::floating_point<std::decay<decltype(T::gamma)>>;
      requires std::floating_point<std::decay<decltype(T::weight)>>;
      requires std::same_as<std::decay<decltype(T::fcid)>, typename T::index_t>;
      requires std::same_as<std::decay<decltype(T::hcid)>, typename T::index_t>;
      requires std::convertible_to<std::decay<decltype(T::active)>, bool>;
  };
  
  template <typename T>
  concept particle_vector = requires {
      true;
  };
  
  template <typename T>
  concept particle_group = requires {
      true;
  };
  
  // ----- Move to PIC library -----
  // ---------------------------------------------------
  template <typename T>
  concept trajectory_integrator = requires {
      T::advance;
  };
  
  template <typename T>
  concept velocity_integrator = requires {
      T::advance;
  };
  
  template <typename T>
  concept particle_pusher = requires {
      T::advance;
  };
  // ---------------------------------------------------
  
  //
} // end namespace tf::concepts

// ===== Particle Data Types =====
// ===============================
namespace tf::types
{
  // ----- Detail -----
  namespace detail
  {
    // >> ParticleBase << (Satisfies concepts::particle)
    template<typename T, template<typename> typename L, template<typename> typename M, template<typename> typename I>
    struct ParticleBase {
      using value_t = T;
      using location_t = L<T>;
      using momentum_t = M<T>;
      using index_t = I<int>;
      
      ParticleBase() = default;
      
      ParticleBase(
          location_t location_,
          location_t old_location_,
          momentum_t momentum_,
          value_t gamma_,
          value_t weight_,
          index_t fcid_,
          index_t hcid_
      ) : location(location_),
          old_location(old_location_),
          momentum(momentum_),
          gamma(gamma_),
          weight(weight_),
          fcid(fcid_),
          hcid(hcid_),
          active(true) {}
      
      location_t location;
      location_t old_location;
      momentum_t momentum;
      value_t gamma;
      value_t weight;
      
      index_t fcid;
      index_t hcid;
      bool active;
      //
    };// end struct ParticleBase
    
    // >> ParticleVectorBase << (Satisfies concepts::particle_vector)
    template <concepts::particle P>
    class ParticleVectorBase {
    public:
      using particle_t = P;
      
      ParticleVectorBase(size_t i) : data(i) {};
      
      P& operator[](int i) { return data[i]; }
      const P& operator[](int i) const { return data[i]; }
      
      void append(particle_t&& particle) { data.insert(data.end(), particle); }
      void append(ParticleVectorBase<particle_t>& particles) { data.insert(data.end(), particles.data.begin(), particles.data.end()); }
      
      void remove(particle_t* particle_) { (*particle_) = std::move(data.pop_back()); }
      void remove(std::vector<particle_t*> particles_) { for (auto& particle_ : particles_) { remove(particle_); } }
      
      void remove_inactive() {
        auto is_active = [](particle_t particle) { return particle.active == true; };
        
        auto pp = std::partition(data.begin(), data.end(), is_active);
        data.erase(pp, data.end());
      }
    
    private:
      std::vector<P> data;
      
      template<typename V, typename M>
      requires std::same_as<typename V::particle_t::location_t::dimension_t, typename M::dimension_t>
      friend class ParticleMapBase;
    };
    
    // >> ParticleMapBase << (Satisfies concepts::particle_map)
    template <typename V, typename Mesh>
    requires std::same_as<typename V::particle_t::location_t::dimension_t, typename Mesh::dimension_t>
    class ParticleMapBase {
      using particle_t = typename V::particle_t;
      using location_t = typename particle_t::location_t;
      using index_t = typename particle_t::index_t;
    
    public:
      ParticleMapBase(V* particles_, Mesh* mesh_)
          : particles(particles_),
            mesh(mesh_),
            full_num_cells(mesh_->full_mesh.nc),
            half_num_cells(mesh_->half_mesh.nc),
            full_data(full_num_cells),
            half_data(half_num_cells)
      {};
      
      const std::vector<particle_t*>& particles_in_full_cell(const int scid) const { return full_data[scid]; }
      const std::vector<particle_t*>& particles_in_half_cell(const int scid) const { return half_data[scid]; }
      
      void update()
      {
        // Reset cell lists
        for (int scid = 0; scid < full_num_cells; scid++) { full_data[scid].clear(); }
        for (int scid = 0; scid < half_num_cells; scid++) { half_data[scid].clear(); }
        
        // Construct new lists
        for (const auto& particle : particles->data) {
          auto constexpr dim_minus_one = index_t::dimension_t::value - 1;
          
          auto& full_cell_vector = full_data[particle.fcid[dim_minus_one]];
          full_cell_vector.insert(full_cell_vector.end(), &particle);
          
          auto& half_cell_vector = half_data[dim_minus_one];
          half_cell_vector.insert(half_cell_vector.end(), &particle);
        }
        //
      };// end update()
      
      void coalesce() { /* fill at some point */ };
    
    private:
      int full_num_cells, half_num_cells;
      V* const particles;
      Mesh* const mesh;
      
      std::vector<std::vector<const particle_t*>> full_data, half_data;
      //
    };// end class ParticleMap
    
    // >> ParticleGroupBase << (Satisfies concepts::particle_group)
    template <concepts::particle P, template<typename> typename V, typename mesh_t>
    //requires concepts::particle_vector<V<P>> and std::same_as<typename P::location_t::dimension_t, typename mesh_t::dimension_t>
    struct ParticleGroupBase { // Satisfies concepts::particle_group
      using particle_t = P;
      using particle_vector_t = V<P>;
      using particle_map_t = ParticleMapBase<particle_vector_t, mesh_t>;
      using value_t = typename particle_t::value_t;
      
      explicit ParticleGroupBase(size_t size_, value_t mass_, value_t charge_, mesh_t* mesh_)
          : particles(size_), mass(mass_), charge(charge_),
            map(&(this->particles), mesh_)
      {};
      
      void add_particle(particle_t&& particle_) { this->particles.append(particle_); }
      void add_particle(particle_vector_t& particles_) { this->particles.append(particles_); }
      
      void remove_particle(particle_t* particle_) { this->particles.remove(particle_); }
      void remove_particle(std::vector<particle_t*>& particles_) { this->particles.remove(particles_); }
      
      void remove_inactive() { this->particles.remove_inactive(); }
      
      [[nodiscard]] bool update_this_step(const value_t simulation_dt, const size_t simulation_step)
      {
        switch (update_interval.type) {
          case control::IntervalType::StepCount:
            return control::updateThisStep(update_interval, simulation_step);
          
          case control::IntervalType::SimTime:
          case control::IntervalType::WallTime:
            return control::updateThisStep(update_interval, simulation_dt);
          
          default:
            return false;
        }
      }
      
      // Data Members
      value_t mass;
      value_t charge;
      
      particle_vector_t particles;
      particle_map_t map;
      control::Interval<value_t> update_interval;
      //
    };// end struct ParticleGroup
    //
  } // end namespace scratch::pic::detail
  
  // ----- Utilities -----
  namespace utilities
  {
  //
  } // end namespace scratch::pic::utilities

// ----- Particle Type Generators -----
template <std::floating_point fp>
using Particle2D3V = detail::ParticleBase<fp, vec2, vec3, vec3>;

template <std::floating_point fp>
using Particle3D = detail::ParticleBase<fp, vec3, vec3, vec4>;

// ----- Particle Vector -----
template <concepts::particle P>
using ParticleVector = detail::ParticleVectorBase<P>; // currently a wrapper for std::vector, can be extended

// ----- Particle Group -----
template <concepts::particle P, typename M>
requires std::same_as<typename P::location_t::dimension_t, typename M::dimension_t>
using ParticleGroup = detail::ParticleGroupBase<P, ParticleVector, M>;

// ----- Integrators & Pusher -----
// >> NullIntegrator << (Satisfies concepts::integrator)
template <std::floating_point fp, concepts::particle_group G>
struct NullIntegrator {
  static inline constexpr void advance([[maybe_unused]] G& group, [[maybe_unused]] fp dt) {
    std::cout << "Called NullIntegrator::advance()." << std::endl;
  }
  //
};// end struct NullIntegrator

// >> TrajectoryIntegrator << (Satisfies concepts::integrator)
template <std::floating_point fp, concepts::particle_group G>
class TrajectoryIntegrator {
public:
  static constexpr size_t dimension_v = G::particle_t::location_t::dimension_t::value;
  
  static inline constexpr void advance([[maybe_unused]] G& group, [[maybe_unused]] fp dt) {
    if constexpr (dimension_v == 2) {
      advance2d(group, dt);
    } else if constexpr (dimension_v == 3) {
      advance3d(group, dt);
    }
  }

private:
  static inline void advance2d([[maybe_unused]] G& group, [[maybe_unused]] fp dt) {
    std::cout << "Called TrajectoryIntegrator::advance2d()." << std::endl;
  }
  static inline void advance3d([[maybe_unused]] G& group, [[maybe_unused]] fp dt) {
    std::cout << "Called TrajectoryIntegrator::advance3d()." << std::endl;
  }
  //
};// end class TrajectoryIntegrator

// >> VelocityIntegrator << (Satisfies concepts::integrator)
template <std::floating_point fp, concepts::particle_group G>
class VelocityIntegrator {
public:
  static constexpr size_t dimension_v = G::particle_t::momentum_t::dimension_t::value;
  
  static inline constexpr void advance([[maybe_unused]] G& group, [[maybe_unused]] fp dt) {
    if constexpr (dimension_v == 2) {
      advance2d(group, dt);
    } else if constexpr (dimension_v == 3) {
      advance3d(group, dt);
    }
  }

private:
  static inline void advance2d([[maybe_unused]] G& group, [[maybe_unused]] fp dt) {
    std::cout << "Called VelocityIntegrator::advance2d()." << std::endl;
  }
  
  static inline void advance3d([[maybe_unused]] G& group, [[maybe_unused]] fp dt) {
    std::cout << "Called VelocityIntegrator::advance3d()." << std::endl;
  }
  //
};// end class VelocityIntegrator

// >> ParticlePusherBase << (Satisfies concepts::particle_pusher)
template <std::floating_point fp, concepts::particle_group G, concepts::trajectory_integrator TI, concepts::velocity_integrator VI>
struct ParticlePusherBase {
  static constexpr void advance(G& group, fp dt) {
    VI::advance(group, dt);
    TI::advance(group, dt);
  };
  //
};// end struct ParticlePusherBase

// ----- Particle Pusher Type Generators -----
template <std::floating_point fp, concepts::particle_group G>
using ParticlePusher = ParticlePusherBase<fp, G, TrajectoryIntegrator<fp, G>, VelocityIntegrator<fp, G>>;

  //
} // end namespace tf::types

// ===== Particle Data Utilities =====
// ===================================
namespace tf::utilities
{
  template <tf::concepts::particle P>
  inline auto computeLorentzFactorP(const typename P::momentum_t& momentum, const typename P::value_t mass) {
    return std::sqrt( 1.0 + momentum.length_squared() / math::SQR(mass * constants::c) );
  } // end computeLorentzFactor()
  
  template <tf::concepts::particle P>
  inline auto computeLorentzFactorV(const typename P::momentum_t& velocity) {
    return 1.0 / std::sqrt( 1.0 - velocity.length_squared() / constants::c_sqr);
  } // end computeLorentzFactor()
  
  template <tf::concepts::particle P>
  inline void update_lorentz_factor(P particle, typename P::value_t mass) {
    particle.gamma = computeLorentzFactorP<P>(particle.momentum, mass);
  } // end update_lorentz_factor()
}

#endif //TRIFORCE_PARTICLE_H
