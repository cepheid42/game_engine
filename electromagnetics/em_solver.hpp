#ifndef EM_SOLVER_HPP
#define EM_SOLVER_HPP


namespace tf::electromagnetics {
  class EMSolver {
    void updateE();
    void updateH();
    void advance();
  };
}

#endif //EM_SOLVER_HPP
