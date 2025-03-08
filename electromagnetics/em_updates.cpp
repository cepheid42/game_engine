#include <cstdint>
#include <iostream>

#include "em_updates.hpp"

namespace tf::electromagnetics {
  void FieldIntegrator::apply() {
    for (std::size_t i = ; i < ; ++i) {
      for (std::size_t j = ; j < ; ++j) {
        for (std::size_t k = ; k < ; ++k) {
          std::cout << "Not Implemented" << std::endl;
        }
      }
    }
  }

  void ExplicitUpdate::apply() {}
}
