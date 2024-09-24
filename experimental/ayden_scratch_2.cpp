#include <concepts>
#include <iostream>
#include <vector>

namespace scratch
{
  namespace tags {
    template <size_t N>
    struct Dimension {
      static constexpr size_t value = N;
    };
  } // end namespace tags
  
  namespace concepts {
    template <typename T>
    concept is_dimensional = requires { typename T::dimension_t; };
  } // end namespace concepts
  
  namespace vector
  {
    namespace detail
    {
      template<typename T, size_t N>
      struct VectorBase {
        using value_t = T;
        using dimension_t = tags::Dimension<N>;
        
        T data[N];
      };
    } // end namespace detail
    
    template<typename T>
    using vec2 = detail::VectorBase<T, 2>;
    
    template<typename T>
    using vec3 = detail::VectorBase<T, 3>;
    
    template<typename T>
    using vec4 = detail::VectorBase<T, 4>;
  } // end namespace vector
  
  namespace pic
  {
    // ===== Concepts =====
    // ====================
    namespace concepts
    {
      // ----- Particle Data -----
      template <typename T>
      concept particle = requires {
        { std::floating_point<typename T::value_t> };
        { scratch::concepts::is_dimensional<typename T::location_t> };
        { scratch::concepts::is_dimensional<typename T::momentum_t> };
        
        { std::same_as<std::decay<decltype(T::location)>, typename T::location_t> };
        { std::same_as<std::decay<decltype(T::old_location)>, typename T::location_t> };
        { std::same_as<std::decay<decltype(T::momentum)>, typename T::momentum_t> };
        { std::floating_point<std::decay<decltype(T::gamma)>> };
        { std::floating_point<std::decay<decltype(T::weight)>> };
        { std::same_as<std::decay<decltype(T::fcid)>, typename T::index_t> };
        { std::same_as<std::decay<decltype(T::hcid)>, typename T::index_t> };
        { std::convertible_to<std::decay<decltype(T::active)>, bool> };
      };
      
      template <typename T>
      concept particle_vector = requires {
        true;
      };
      
      template <typename T>
      concept particle_group = requires {
        true;
      };
      
      // ----- Particle Pusher -----
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
      
      //
    } // end namespace scratch::pic::concepts
    
    // ===== Base Types =====
    // ======================
    namespace detail
    {
      template<typename T,
          template<typename> typename L,
          template<typename> typename M,
          template<typename> typename I>
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
        
        // Data members
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
      
      template <concepts::particle P>
      using ParticleVectorBase = std::vector<P>;
      
      template <concepts::particle P, template<typename> typename V>
      requires concepts::particle_vector<V<P>>
      struct ParticleGroupBase {
        using particle_t = P;
        using particle_vector_t = V<P>;
        using value_t = typename particle_t::value_t;
        
        ParticleGroupBase(size_t size_, value_t mass_, value_t charge_)
        : particles(size_), mass(mass_), charge(charge_)
        {}
        
        // Data Members
        particle_vector_t particles;
        typename particle_t::value_t mass;
        typename particle_t::value_t charge;
        
        //
      };// end struct ParticleGroup
      //
    } // end namespace scratch::pic::detail
    
    // ===== Particle Types ======
    // ===========================
    template <std::floating_point fp>
    using Particle2D3V = detail::ParticleBase<fp, vector::vec2, vector::vec3, vector::vec3>;
    
    template <std::floating_point fp>
    using Particle3D = detail::ParticleBase<fp, vector::vec3, vector::vec3, vector::vec4>;
    
    // ===== Particle Vector =====
    // ===========================
    template <concepts::particle P>
    using ParticleVector = detail::ParticleVectorBase<P>; // currently a wrapper for std::vector, can be extended
    
    // ===== Particle Group =====
    // ==========================
    template <concepts::particle P>
    using ParticleGroup = detail::ParticleGroupBase<P, ParticleVector>;
    
    // ===== Particle Pusher =====
    // ===========================
    template <std::floating_point fp, concepts::particle_group G>
    struct NullIntegrator {
      static constexpr void advance(G&, fp) {
        std::cout << "Called NullIntegrator::advance()." << std::endl;
      }
    };
    
    template <std::floating_point fp, concepts::particle_group G>
    class TrajectoryIntegrator {
    public:
      static constexpr size_t dimension_v = G::particle_t::location_t::dimension_t::value;
      
      static constexpr void advance([[maybe_unused]] G& group, [[maybe_unused]] fp dt) {
        if constexpr (dimension_v == 2) {
          advance2d(group, dt);
        } else if constexpr (dimension_v == 3) {
          advance3d(group, dt);
        }
      }
      
    private:
      static void advance2d([[maybe_unused]] G& group, [[maybe_unused]] fp dt) {
        std::cout << "Called TrajectoryIntegrator::advance2d() dt = ." << dt << std::endl;
      }
      static void advance3d([[maybe_unused]] G& group, [[maybe_unused]] fp dt) {
        std::cout << "Called TrajectoryIntegrator::advance3d(). dt = " << dt << std::endl;
      }
      
    };
    
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
        std::cout << "Called VelocityIntegrator::advance2d(). dt = " << dt << std::endl;
      }
      
      static inline void advance3d([[maybe_unused]] G& group, [[maybe_unused]] fp dt) {
        std::cout << "Called VelocityIntegrator::advance3d(). dt = " << dt << std::endl;
      }
    };
    
    template <std::floating_point fp, concepts::particle_group G, concepts::trajectory_integrator TI, concepts::velocity_integrator VI>
    struct ParticlePusherBase {
      static constexpr void advance(G& group, fp dt) {
        std::cout << dt << std::endl;
        VI::advance(group, dt);
        TI::advance(group, dt);
      };
    };

    template <std::floating_point fp, concepts::particle_group G>
    using ParticlePusher = ParticlePusherBase<fp, G, TrajectoryIntegrator<fp, G>, VelocityIntegrator<fp, G>>;
    
    //
  } // end namespace scratch::pic
  //
} // end namespace scratch

int main()
{
  using fptype = double;
  
  using Particle = scratch::pic::Particle3D<fptype>;
//  using Particle = scratch::pic::Particle2D3V<fptype>;

  using ParticleGroup = scratch::pic::ParticleGroup<Particle>;
  using ParticlePusher = scratch::pic::ParticlePusher<fptype, ParticleGroup>;
  
  ParticleGroup particleGroup(10, 1.0e-3, -1.0e-19);
  
  ParticlePusher::advance(particleGroup, 1.0e-12);
  
  return 0;
}