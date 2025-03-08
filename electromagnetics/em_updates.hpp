#ifndef EM_UPDATES_HPP
#define EM_UPDATES_HPP

namespace tf::electromagnetics {
  struct FieldIntegrator {
    static void apply();
  };

  struct ExplicitUpdate {
    static void apply();
  };

}



#endif //EM_UPDATES_HPP
