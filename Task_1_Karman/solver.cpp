#include "solver.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include <math.h>

#include "vector.hpp"
#include "system/system_helpers.hpp"

solver::solver( const vector &_y_0, double s_initial, const std::function<vector( double, vector &, double )> &f,
                double ksi0, double ksin, double _h ) {
  y_0 = _y_0;
  func = f;
  ksi_0 = ksi0;
  ksi_n = ksin;
  h = _h;
  s = s_initial;
  n_cells = (ksi_n - ksi_0) / h;
  type = solver_types::Runge_Kutta_method;
}

auto solver::solve_system(const std::string& dir = "./results/") -> std::vector<double> {
  vector result(y_0.size());
  std::string directory = mk_dir(s, dir);

  if (type == solver_types::Runge_Kutta_method) {
    result = Runge_Kutta_solver(y_0, s, true, false, directory);
  }
  if (type == solver_types::Shooting_method) {
    result = Shooting_method(y_0, true, false);
  }
  std::string python_script = "python3 ./python_scripts/plots.py  " + directory;
  if (system(python_script.c_str()) == -1) {
    throw std::runtime_error("bad script");
  }
  return result;
}

vector solver::find_infinity(const vector &y, double &speed, const std::string& dir, bool output_flag ) {
  vector result(y_0.size());
  vector result_1(y_0.size());
  std::ofstream out_parameters;
  std::ofstream shooting_output;
  std::string directory;
  if (output_flag) {
    directory = mk_dir(s, dir);
    out_parameters.open(directory + "/common_parameters.csv");
    shooting_output.open(directory + "/shooting_parameters.csv");
    out_parameters << std::scientific;
  }

  ksi_n = 5;
  n_cells = (ksi_n - ksi_0) / h;

  result = Shooting_method(y_0, false, false);
  if (output_flag) {
    out_parameters << ksi_n << "\t" << result[0] << "\t" << result[1] << "\t" << result[2] << "\t"
                   << result[3] << "\t" << result[4] << std::endl;
    shooting_output << ksi_n << "\t" << y_0_new[1] << "\t" << y_0_new[3] << std::endl;
  }

  y_0_old = y_0;

  do {
    if (ksi_n >= 10) {
      ksi_n += 1;
    } else {
      ksi_n += 5;
    }
    n_cells = (ksi_n - ksi_0) / h;
    y_0_old = y_0_new;
    result_1 = Shooting_method(y_0_new, false, false);
    if (output_flag) {
      out_parameters << ksi_n << "\t" << result_1[0] << "\t" << result_1[1] << "\t" << result_1[2] << "\t"
                     << result_1[3] << "\t" << result_1[4] << std::endl;
      shooting_output << ksi_n << "\t" << y_0_new[1] << "\t" << y_0_new[3] << std::endl;
    }
    result = result_1;
//    std::cout << ksi_n << " " << fabs(y_0_old[1] - y_0_new[1]) << " " << fabs(y_0_old[3] - y_0_new[3]) << std::endl;
    std::cout << ksi_n << " " << fabs(y_0_old[1] - y_0_new[1]) << " " << fabs(y_0_old[3] - y_0_new[3]) << std::endl;
  } while (fabs(y_0_old[1] - y_0_new[1]) + fabs(y_0_old[3] - y_0_new[3]) > 1e-6);

  std::cout << ksi_n << " " << y_0_new[1] << " " << y_0_new[3] << std::endl;

  if (output_flag) {
    Runge_Kutta_solver(y_0_new, s, true, true, directory);
    std::string python_script = "python3 ./python_scripts/plots.py  " + directory;
    if (system(python_script.c_str()) == -1) {
      throw std::runtime_error("bad script");
    }
  }
  return y_0_new;
}

vector solver::Runge_Kutta_solver( vector temp_y,
                                   double speed,
                                   bool output_flag,
                                   bool stream_flag,
                                   std::string dir = "" ) {
  vector k[5];

  std::ofstream out;
  if (output_flag) {
    out.open(dir + "result.csv");
    out << std::scientific;
    std::ofstream out_par(dir + "result_parameters.csv");
    out_par << ksi_n << " " << h << " " << speed << " " << y_0[1] << " " << y_0[3] << std::endl;
    out << ksi_0 << " " << y_0[0] << " " << y_0[1] << " " << y_0[2] << " " << y_0[3] << " " << y_0[4] << std::endl;
//    out << ksi_0 << "\t" << y_0[0] << "\t" << y_0[1] << "\t" << y_0[2] << "\t" << y_0[3] << "\t" << y_0[4] << std::endl;

  }
  if (stream_flag) {
    res.resize(n_cells + 1);
    res[0].F = y_0[0];
    res[0].G = y_0[2];
    res[0].H = y_0[4];
  }
  vector cur_y = temp_y;

  for (int i = 1; i <= n_cells; ++i) {
    double ksi = ksi_0 + i * h;
    k[1] = h * func(ksi, cur_y, speed);
    cur_y = temp_y + 0.5 * k[1];
    k[2] = h * func(ksi + 0.5 * h, cur_y, speed);
    cur_y = temp_y + 0.5 * k[2];
    k[3] = h * func(ksi + 0.5 * h, cur_y, speed);
    cur_y = temp_y + k[3];
    k[4] = h * func(ksi + h, cur_y, speed);
    k[0] = (k[1] + 2 * k[2] + 2 * k[3] + k[4]) / 6.;
    temp_y = temp_y + k[0];
    if (output_flag) {
/*
      out << ksi << "\t" << temp_y[0] << "\t" << temp_y[1] << "\t" << temp_y[2] << "\t" << temp_y[3] << "\t"
          << temp_y[4] << std::endl;
*/

      out << ksi << " " << temp_y[0] << " " << temp_y[1] << " " << temp_y[2] << " " << temp_y[3] << " "
          << temp_y[4] << std::endl;
    }
    if (stream_flag) {
      res[i].F = temp_y[0];
      res[i].G = temp_y[2];
      res[i].H = temp_y[4];
    }
  }
  if (output_flag) {
    out.close();
  }
  return temp_y;
}

vector solver::Shooting_method( const vector &y_new, bool output_flag, bool stream_flag ) {

  vector result(y_new.size());

  vector derivative_result(y_new.size());
  vector shooting_y(y_new.size());
  vector derivative_y(y_new.size());
  shooting_y = y_new;

  double F, G;

  double h_derivative = 1e-4;
  double f_alpha, f_beta, g_alpha, g_beta;
  double j11, j12, j21, j22;
  double delta_alpha;
  double delta_beta;
  result = Runge_Kutta_solver(shooting_y, s, false, false);
  F = result[0];
  G = result[2];

  while (/*sqrt(pow(F, 2.) + pow(G - s, 2.)) > eps*/ fabs(F) + fabs(G - s) > eps) {

    derivative_y = shooting_y;
    derivative_y[1] += h_derivative;
    derivative_result = Runge_Kutta_solver(derivative_y, s, false, false);
    f_alpha = (derivative_result[0] - F) / h_derivative;
    g_alpha = (derivative_result[2] - G) / h_derivative;

    derivative_y = shooting_y;
    derivative_y[3] += h_derivative;
    derivative_result = Runge_Kutta_solver(derivative_y, s, false, false);
    f_beta = (derivative_result[0] - F) / h_derivative;
    g_beta = (derivative_result[2] - G) / h_derivative;

    double det = f_alpha * g_beta - g_alpha * f_beta;

    j11 = g_beta / det;
    j12 = -f_beta / det;
    j21 = -g_alpha / det;
    j22 = f_alpha / det;

    delta_alpha = j11 * F + j12 * (G - s);
    delta_beta = j21 * F + j22 * (G - s);

    shooting_y[1] -= delta_alpha;
    shooting_y[3] -= delta_beta;

    result = Runge_Kutta_solver(shooting_y, s, false, false);
    F = result[0];
    G = result[2];

    std::cout << ksi_n << " " << shooting_y[1] << " " << shooting_y[3] << " " << result[0] << " " << result[2]
              << std::endl;
  }

  result = Runge_Kutta_solver(shooting_y, s, output_flag, stream_flag);
  y_0_new = shooting_y;
  return result;
}

vector solver::find_infinity_second_branch( const vector &y, double &speed, std::string dir, bool output_flag ) {
  vector result(y.size());
  vector result_1(y.size());
  std::ofstream out_parameters;
  std::ofstream shooting_output;
  std::string directory;

  ksi_n = 10;
  n_cells = (ksi_n - ksi_0) / h;

  result = Shooting_method_second_branch(y, speed, speed, false);
  y_0_old = y;

  do {
      ksi_n += 1;

    n_cells = (ksi_n - ksi_0) / h;
    y_0_old = y_0_new;
    result_1 = Shooting_method_second_branch(y_0_new, speed, false, false);
    result = result_1;
    std::cout << ksi_n << " " << fabs(y_0_old[1] - y_0_new[1]) << " " << fabs(y_0_old[3] - y_0_new[3]) << std::endl;
  } while (fabs(y_0_old[1] - y_0_new[1]) + fabs(y_0_old[3] - y_0_new[3]) > 1e-6 && ksi_n < 25);
  std::cout << result[0] << " " << result[1] << " " << result[2] << " " << result[3] << " " << result[4] << std::endl;

  if (isnan(result[4])) {
    throw std::runtime_error("invalid data");
  }

  if (output_flag) {
    directory = mk_dir(speed, dir);
    out_parameters.open(directory + "/common_parameters.csv");
    shooting_output.open(directory + "/shooting_parameters.csv");
    out_parameters << std::scientific;

    shooting_output << ksi_n << "\t" << y_0_new[1] << "\t" << y_0_new[3] << std::endl;

    Runge_Kutta_solver(y_0_new, speed, true, true, directory);
    std::string python_script = "python3 ./python_scripts/plots.py  " + directory;
    if (system(python_script.c_str()) == -1) {
      throw std::runtime_error("bad script");
    }
  }

  return y_0_new;
}

vector solver::Shooting_method_second_branch( const vector &y_new, double &speed, bool output_flag, bool stream_flag ) {

  vector result(y_new.size());

  vector derivative_result(y_new.size());
  vector shooting_y(y_new.size());
  vector derivative_y(y_new.size());
  shooting_y = y_new;

  double F, G;

  double h_derivative = 1e-5;
  double f_alpha, f_s, g_alpha, g_s;
  double j11, j12, j21, j22;
  double delta_alpha;
  double delta_s;
  double derivative_s = speed;
  result = Runge_Kutta_solver(shooting_y, speed, false, false);
//  std::cout << result[0] << " " << result[1] << " " << result[2] << " " << result[3] << " " << result[4] << std::endl;
  F = result[0];
  G = result[2];

  int k = 0;
  while (/*sqrt(pow(F, 2.) + pow(G - s, 2.)) > eps*/ fabs(F) + fabs(G - s) > eps) {

    derivative_y = shooting_y;
    derivative_y[1] += h_derivative;
    derivative_result = Runge_Kutta_solver(derivative_y, speed, false, false);
    f_alpha = (derivative_result[0] - F) / h_derivative;
    g_alpha = (derivative_result[2] - G) / h_derivative;

    derivative_y = shooting_y;
    derivative_result = Runge_Kutta_solver(derivative_y, speed + h_derivative, false, false);
    f_s = (derivative_result[0] - F) / h_derivative;
    g_s = (derivative_result[2] - G) / h_derivative;

    double det = f_alpha * g_s - g_alpha * f_s;

    j11 = g_s / det;
    j12 = -f_s / det;
    j21 = -g_alpha / det;
    j22 = f_alpha / det;

    delta_alpha = j11 * F + j12 * (G - s);
    delta_s = j21 * F + j22 * (G - s);

    shooting_y[1] -= delta_alpha;
    speed -= delta_s;

    result = Runge_Kutta_solver(shooting_y, speed, false, false);
    F = result[0];
    G = result[2];
    std::cout << ksi_n << " " << shooting_y[1] << " " << speed << std::endl;

    ++k;
    if (k > 15) {
      throw std::runtime_error("invalid data");
    }
//    std::cout << result[0] << " " << result[1] << " " << result[2] << " " << result[3] << " " << result[4] << std::endl;

  }
  result = Runge_Kutta_solver(shooting_y, speed, output_flag, stream_flag);
  y_0_new = shooting_y;
  return result;
}

/*void solver::grid_search(vector &y, double &speed, std::string dir, bool output_flag ) {
  double alpha;
  double beta;
  alpha = 0.47;
  beta  = -0.56;
  while (speed > 0) {
    while (alpha >= 0.47 && alpha <= 0.49 ) {
      while (beta >= -0.56 && beta <= -0.57) {
        y[1] = alpha;
        y[3] = beta;
        find_infinity_second_branch(y, speed, dir, output_flag);
        beta += 0.001;
      }
      alpha += 0.001;
    }
    speed += 0.001;

  }
}*/

/*
void solver::solve_system_RK_method( vector y, bool output_flag, bool stream_flag ) {
  std::string directory = mk_dir(s);
  vector result(y.size());
  Runge_Kutta_solver(y, output_flag, stream_flag, directory);

  std::string python_script = "python3 ./python_scripts/plots.py  " + directory;

  if (system(python_script.c_str()) == -1) {
    throw std::runtime_error("bad script");
  }
}

*/
