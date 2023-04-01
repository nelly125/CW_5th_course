#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>

void ADT_solver( const uint32_t n,
                 const uint32_t m,
                 const double H,
                 double Mach,
                 std::ofstream &output_psi,
                 std::ofstream &output_cp ) {
  double eps     = 0.1;
  double x_left  = -1;
  double x_right = 1;

  uint32_t N = n;
  double   dx, dy;

  dx = (x_right - x_left) / N;
  uint32_t M = m;
  dy = H / M;

  double M_infty = Mach;

  std::vector<std::vector<double>> psi(N + 1, std::vector<double>(M + 1));

  std::vector<double> A(M + 1);
  std::vector<double> B(M + 1);
  std::vector<double> C(M + 1);
  std::vector<double> F(M + 1);
  std::vector<double> alpha(M + 1);
  std::vector<double> beta(M + 1);

  for (uint32_t i = 0; i < N + 1; ++i) {

    double x = x_left + i * dx;

    psi[i][0] = fabs(x) <= 1 ? eps * (1 - pow(x, 2.)) : 0;

    B[0] = 1.;
    C[0] = 0;
    F[0] = fabs(x) <= 1 ? eps * (1 - pow(x, 2.)) : 0;

    alpha[1] = -C[0] / B[0];
    beta[1]  = F[0] / B[0];

    for (uint32_t j = 1; j < M; ++j) {

      A[j] = 1;
      B[j] = (1 - pow(M_infty, 2.)) * pow(dy, 2.) / pow(dx, 2.) - 2;
      C[j] = 1;
      F[j] = (1 - pow(M_infty, 2.)) * (pow(dy, 2.) / pow(dx, 2.)) *
          (2 * ((i > 0) ? psi[i - 1][j] : 0) - ((i > 1) ? psi[i - 2][j] : 0));

      alpha[j + 1] = -C[j] / (A[j] * alpha[j] + B[j]);
      beta[j + 1]  = (F[j] - A[j] * beta[j]) / (A[j] * alpha[j] + B[j]);
    }

    psi[i][M] = 0;
    for (int j = M - 1; j > 0; --j) {
      psi[i][j] = alpha[j + 1] * psi[i][j + 1] + beta[j + 1];
    }

  }

/*  output_cp << "x\tC_p\texact" << std::endl;
  for (uint32_t i = 0; i < N + 1; ++i) {
    double x          = x_left + i * dx;
    double exact_beta = sqrt(pow(M_infty, 2.) - 1);
    double C_p        = -2 * (1 / (pow(M_infty, 2.) - 1) * (psi[i][1] - psi[i][0]) / (dy));
    double exact_C_p  = (-4 * eps * x) / exact_beta;
    output_cp << x << "\t" << C_p << "\t" << exact_C_p << std::endl;
  }*/

  output_cp << std::scientific << std::setprecision(6);
  output_cp << "x\tC_p\texact" << std::endl;
  for (uint32_t i = 0; i < N + 1; ++i) {
    double x          = x_left + i * dx;
    double exact_beta = sqrt(pow(M_infty, 2.) - 1);
    double C_p        = -2 * (1 / (pow(M_infty, 2.) - 1) * (psi[i][1] - psi[i][0]) / (dy));
    double exact_C_p  = (-4 * eps * x) / exact_beta;
    output_cp << x << "\t" << C_p << "\t" << exact_C_p << std::endl;
  }
  output_psi << "x\ty\tpsi\trho" << std::endl;
  for (uint32_t j = 1; j < M; ++j) {
    for (uint32_t i = 0; i < N + 1; ++i) {
      double x   = x_left + i * dx;
      double y   = j * dy;
      double rho = -pow(M_infty, 2.) * (1 / (pow(M_infty, 2.) - 1) * (psi[i][j + 1] - psi[i][j - 1]) / (2 * dy));
      output_psi << x << "\t" << y << "\t" << psi[i][j] << "\t" << rho << std::endl;

    }
  }
}

void H_dependency() {

  std::vector<double> H = {1. / 3, 0.5, 1, 5, 10};

  for (auto h : H) {
    std::ofstream output_psi("./results/ADT_output_psi_H_" + std::to_string(h) + ".txt");
    std::ofstream output_cp("./results/ADT_output_H_" + std::to_string(h) + ".txt");

    ADT_solver(1000, 1000, h, 20, output_psi, output_cp);
  }
}

void N_M_dependency() {

  std::vector<uint32_t> N = {10000};
  std::vector<uint32_t> M = {100, 500, 1000};

  for (auto n : N) {
    for (auto m : M) {
    std::ofstream output_cp("./results/ADT_output_N_" + std::to_string(n) + "_M_" + std::to_string(m) + ".txt");

    ADT_solver(n, n, 1. / 3, sqrt(2), output_cp, output_cp);
    }
  }
}

int main() {
//  N_M_dependency();
  H_dependency();
  return 0;
}
