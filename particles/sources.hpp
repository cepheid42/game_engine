#ifndef SOURCES_HPP
#define SOURCES_HPP

#include "particles.hpp"
#include "vec3.hpp"

namespace tf::particles
{

struct ParticleBeam {
   std::vector<Particle> particles;
   vec3<double> velocity;
   vec3<double> temperature;
   vec3<double> kvec;
   vec3<double> ppc;
   double density;
   double mass;
   std::array<std::size_t, 6> extents;

   ParticleBeam(const auto& vel_, const auto& temp_, const auto& kvec_, const auto& ppc_, const auto den_, const auto mass_, const auto& extents_)
   : velocity(vel_), temperature(temp_), kvec(kvec_), ppc(ppc_), density(den_), mass(mass_), extents(extents_)
   {
      initParticles();
   }

   void initParticles() {
      const auto weight = density * dx * dy * dz / (ppc.x * ppc.y * ppc.z);
      const auto gamma = calculateGammaV(velocity);
      for (std::size_t ci = extents[0]; ci < extents[1]; ++ci) {
         for (std::size_t cj = extents[2]; cj < extents[3]; ++cj) {
            for (std::size_t ck = extents[4]; ck < extents[5]; ++ck) {
               const auto cell_x = x_range[0] + (static_cast<double>(ci) * dx);
               const auto cell_y = y_range[0] + (static_cast<double>(cj) * dy);
               const auto cell_z = z_range[0] + (static_cast<double>(ck) * dz);

               for (std::size_t i = 0; i < ppc.x; ++i) {
                  for (std::size_t j = 0; j < ppc.y; ++j) {
                     for (std::size_t k = 0; k < ppc.z; ++k) {
                        const vec3 loc {
                           cell_x + dx * (static_cast<double>(i) + 0.5) / ppc.x,
                           cell_y + dy * (static_cast<double>(j) + 0.5) / ppc.y,
                           cell_z + dz * (static_cast<double>(k) + 0.5) / ppc.z
                        };

                        particles.emplace_back(
                           loc,
                           loc,
                           velocity,
                           weight,
                           gamma,
                           false
                        );
                     } // end for(k)
                  } // end for(j)
               } // end for(i)
            } // end for(ck)
         } // end for(cj)
      } // end for(ci)
   } // end initParticles()

   void apply(auto& gpart) {
      std::ranges::copy(particles, gpart.end());
   }
}; // end struct ParticleBeam


} // end namespace tf::particles



#endif //SOURCES_HPP
