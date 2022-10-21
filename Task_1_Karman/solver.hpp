#pragma once

#include <vector>
#include <functional>
#include <cmath>
#include <fstream>

enum class solver_types { Shooting_method, Runge_Kutta_method };

struct polar_parameters {
  double r;
  double phi;
  double zeta;

  friend std::ofstream &operator<<(std::ofstream &os, const polar_parameters &data) {
    os << data.r * cos(data.phi) << "\t" << data.r * sin(data.phi) << "\t" << data.zeta << std::endl;
    return os;
  }

};
//std::ofstream& operator<<(std::ofstream& os, const polar_parameters& data)


struct solution {
  double F;
  double G;
  double H;
};

using vector = std::vector<double>;

class solver {
 public:
  solver(const vector &y_0, double s, const std::function<vector(double, vector &, double)> &f, double ksi0,
         double ksin, int n);

  auto solve_system() -> std::vector<double>;
  void solve_system_RK_method(vector temp_y, bool output_flag, bool stream_flag);
  void solve_system_Shooting_method(vector temp_y, bool output_flag, bool stream_flag);
  double find_infinity();
  void find_streamlines();

 private:
  double ksi_0, ksi_n;
  double s;
  int n_cells;
  double h;
  std::vector<double> y;
  vector y_0_new;
  vector y_0_old;

  double eps = 1e-4;

  std::function<vector(double, vector &, double)> func;

  solver_types type;
  bool plot_flag = true;

  std::vector<solution> res;

  vector Runge_Kutta_solver(vector temp_y, bool output_flag, bool stream_flag);

  vector Shooting_method(const vector &tmp_y, bool output_flag, bool stream_flag);

};
