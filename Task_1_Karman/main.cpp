#include <iostream>
#include <vector>
#include <functional>
#include <cmath>

#include "solver.hpp"

vector func( double ksi, const vector &y, double s ) {
  std::vector<double> res(y.size());
  res[0] = y[1];
  res[1] = y[4] * y[1] + pow(y[0], 2.) - pow(y[2], 2.) + pow(s, 2.);
  res[2] = y[3];
  res[3] = 2 * y[0] * y[2] + y[4] * y[3];
  res[4] = -2 * y[0];
  return res;
}

int main() {
//    double alpha = 0.510214, beta = -0.61591;
  double alpha = 0.39, beta = -0.41;

  double s = 0.5;
  double ksi_0 = 0.;
  double ksi_n = 3;

  int n = (ksi_n - ksi_0) / 1e-3;
  std::vector<double> initial_conditions = {0., alpha, 1., beta, 0};
  std::vector<double> result = func(ksi_0, initial_conditions, s);
//    for (auto c: result) {
//        std::cout << c << std::endl;
//    }

  solver Solver(initial_conditions, s, func, ksi_0, ksi_n, n);
  Solver.solve_system();
//    Solver.find_infinity();
//  Solver.find_streamlines();
  return 0;

}
