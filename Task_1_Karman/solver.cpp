#include "solver.hpp"

#include <utility>
#include <fstream>
#include <iostream>
#include <cmath>

#include "vector.hpp"

#define EPS 1e-5

solver::solver( const vector &y_0, double s_initial, const std::function<vector( double, vector &, double )> &f,
                double ksi0, double ksin, int n ) {
  y = y_0;
  func = f;
  ksi_0 = ksi0;
  ksi_n = ksin;
  n_cells = n;
  s = s_initial;
  h = (ksi_n - ksi_0) / n;
  type = solver_types::Runge_Kutta_method;
}

auto solver::solve_system() -> std::vector<double> {
  vector result(y.size());
  if (type == solver_types::Runge_Kutta_method) {
    result = Runge_Kutta_solver(y);
  }
  if (type == solver_types::Shooting_method) {
    result = Shooting_method(y);
  }
  return result;
}

double solver::find_infinity() {
  vector result(y.size());
  vector result_1(y.size());
  std::ofstream out_parameters("./results/common_parameters.txt");
  std::ofstream shooting_output("../results/shooting_parameters.txt");
  out_parameters << std::scientific;
  ksi_n = 5;
  result = Shooting_method(y);
  out_parameters << ksi_n << "\t" << result[0] << "\t" << result[1] << "\t" << result[2] << "\t"
                                                 << result[3] << "\t" << result[4] << std::endl;
  double norm;
  y_0_old = y;

  do {
    ksi_n += 5;
/*    if (h > 1e-5) {
      h /= 5;
    }*/
    n_cells = (ksi_n - ksi_0) / h;
/*    if (eps >= 1e-) {
      eps /= 10;
    }*/

    y_0_old = y_0_new;
    result_1 = Shooting_method(y_0_new);

    out_parameters << ksi_n << "\t" << result_1[0] << "\t" << result_1[1] << "\t" << result_1[2] << "\t"
                   << result_1[3] << "\t" << result_1[4] << std::endl;
    shooting_output << ksi_n << "\t" << y_0_new[1] << "\t" << y_0_new[3] << std::endl;
    norm = vector_norm(result, result_1);
    result = result_1;

//    std::cout << ksi_n << " " << fabs(y_0_old[1] - y_0_new[1]) << " " << fabs(y_0_old[3] - y_0_new[3]) << std::endl;
  } while (fabs(y_0_old[1] - y_0_new[1]) + fabs(y_0_old[3] - y_0_new[3]) > 1e-6);

  output_flag = true;
  stream_flag = true;
  Runge_Kutta_solver(y_0_old);

  return ksi_n;
}

vector solver::Runge_Kutta_solver( vector temp_y ) {
  vector k[5];

  std::ofstream out;
  if (output_flag) {
    out.open("./results/result.txt");
    out << std::scientific;
    std::ofstream out_par("./results/result_parameters.txt");
    out_par << ksi_n << " " << s << " " << y[1] << " " << y[3] << std::endl;
    out << ksi_0 << "\t" << y[0] << "\t" << y[1] << "\t" << y[2] << "\t" << y[3] << "\t" << y[4] << std::endl;
  }
  if (stream_flag) {
    res.resize(n_cells + 1);
    res[0].F = y[0];
    res[0].G = y[2];
    res[0].H = y[4];
  }
  vector cur_y = temp_y;

  for (int i = 1; i <= n_cells; ++i) {
    double ksi = ksi_0 + i * h;
    k[1] = h * func(ksi, cur_y, s);
    cur_y = temp_y + 0.5 * k[1];
    k[2] = h * func(ksi + 0.5 * h, cur_y, s);
    cur_y = temp_y + 0.5 * k[2];
    k[3] = h * func(ksi + 0.5 * h, cur_y, s);
    cur_y = temp_y + k[3];
    k[4] = h * func(ksi + h, cur_y, s);
    k[0] = (k[1] + 2 * k[2] + 2 * k[3] + k[4]) / 6.;
    temp_y = temp_y + k[0];
    if (output_flag) {
      out << ksi << "\t" << temp_y[0] << "\t" << temp_y[1] << "\t" << temp_y[2] << "\t" << temp_y[3] << "\t"
          << temp_y[4] << std::endl;
    }
    if (stream_flag) {
      res[i].F = temp_y[0];
      res[i].G = temp_y[2];
      res[i].H = temp_y[4];
    }
  }
  if (output_flag && plot_flag) {
    system("python3 ./python_scripts/plots.py");
  }
  return temp_y;
}

vector solver::Shooting_method( const vector &y_new ) {
  output_flag = false;
  vector result(y_new.size());

  vector derivative_result(y_new.size());
  vector shooting_y(y_new.size());
  vector derivative_y(y_new.size());
  shooting_y = y_new;

  double F, G;

  double h_derivative = 1e-3;
  double f_alpha, f_beta, g_alpha, g_beta;
  double j11, j12, j21, j22;
  double delta_alpha;
  double delta_beta;
  result = Runge_Kutta_solver(shooting_y);
  F = result[0];
  G = result[2];

  while (sqrt(pow(F, 2.) + pow(G - s, 2.)) > eps /*fabs(F) + fabs(G - s) > eps*/) {
    derivative_y = shooting_y;
    derivative_y[1] += h_derivative;
    derivative_result = Runge_Kutta_solver(derivative_y);
    f_alpha = (derivative_result[0] - F) / h_derivative;
    g_alpha = (derivative_result[2] - G) / h_derivative;

    derivative_y = shooting_y;
    derivative_y[3] += h_derivative;
    derivative_result = Runge_Kutta_solver(derivative_y);
    f_beta = (derivative_result[0] - F) / h_derivative;
    g_beta = (derivative_result[2] - G) / h_derivative;

    double det = f_alpha * g_beta - g_alpha * f_beta;

    j11 = g_beta / det;
    j12 = -f_beta / det;
    j21 = -g_alpha / det;
    j22 = f_alpha / det;

    delta_alpha = j11 * F + j12 * G;
    delta_beta = j21 * F + j22 * G;

    shooting_y[1] -= delta_alpha;
    shooting_y[3] -= delta_beta;

    result = Runge_Kutta_solver(shooting_y);
    F = result[0];
    G = result[2];

    std::cout << shooting_y[1] << " " << shooting_y[3] << std::endl;
  }

  output_flag = true;
  result = Runge_Kutta_solver(shooting_y);
  y_0_new = shooting_y;
  return result;
}

void solver::find_streamlines() {
  std::string file_name;
  polar_parameters psk{};
  polar_parameters new_psc{};

  find_infinity();

  int N = 1;
  for (int j = 0; j < N; ++j) {
    psk.r = 1;
    psk.phi = 2 * M_PI * j / N;
    psk.zeta = ksi_n;

    file_name = "../results/streamlines";
    file_name += std::to_string(j);
    file_name += ".txt";
    std::ofstream out_stream(file_name);
    out_stream << std::scientific;

    for (int i = n_cells - 1; i > 0; --i) {
//      std::cout << res[i + 1].F << " " << res[i + 1].G << " " << res[i + 1].H << std::endl;
//      std::cout << res[i].F << " " << res[i].G << " " << res[i].H << std::endl;

//      new_psc.r = psk.r * (2 * res[i].H - h * res[i ].F) / (2 * res[i + 1].H - h * res[i + 1].F) * res[i + 1].H / res[i ].H;
//      new_psc.phi = psk.phi + 0.5 *(res[i + 1].G/ res[i + 1].H + res[i].G / res[i].H) * h;
//      new_psc.zeta = psk.zeta - h;

      new_psc.r = psk.r * (2 * res[i + 1].H - h * res[i + 1].F) / (2 * res[i].H - h * res[i].F) * res[i].H / res[i + 1].H;
      new_psc.phi = psk.phi + 0.5 *(res[i + 1].G/ res[i + 1].H + res[i].G / res[i].H) * h;
      new_psc.zeta = psk.zeta - h;

      psk = new_psc;
      if (i % 100 == 0)
        out_stream << psk;
    }
  }



}

