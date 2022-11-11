#include <vector>
#include <cmath>
#include <iostream>

#include "solver.hpp"
#include "stream.hpp"
#include <string>
#include <filesystem>

vector func( double ksi, const vector &y, double s ) {
  std::vector<double> res(y.size());
  res[0] = y[1];
  res[1] = y[4] * y[1] + pow(y[0], 2.) - pow(y[2], 2.) + pow(s, 2.);
  res[2] = y[3];
  res[3] = 2 * y[0] * y[2] + y[4] * y[3];
  res[4] = -2 * y[0];
  return res;
}

void solution_eps_0( double &alpha, double &beta, double &s ) {
//  double alpha = 0.510214, beta = -0.61591;
  alpha = 0.5;
  beta = -0.6;
  s = 0;
}

void solution_eps_0_1( double &alpha, double &beta, double &s ) {
  alpha = 0.5;
  beta = -0.6;
  s = 0.1;
}

void solution_eps_0_2( double &alpha, double &beta, double &s ) {
  alpha = 0.5;
  beta = -0.6;
  s = 0.2;
}

void solution_eps_0_3( double &alpha, double &beta, double &s ) {
  alpha = 0.47;
  beta = -0.53;
  s = 0.3;
}

void solution_eps_0_4( double &alpha, double &beta, double &s ) {
  alpha = 0.5;
  beta = -0.5;
  s = 0.4;

}

void solution_eps_0_5( double &alpha, double &beta, double &s ) {
  alpha = 0.5;
  beta = -0.5;
  s = 0.5;
}

void solution_eps_0_6( double &alpha, double &beta, double &s ) {
  alpha = 0.4;
  beta = -0.4;
  s = 0.6;
}

void solution_eps_0_7( double &alpha, double &beta, double &s ) {
  alpha = 0.3;
  beta = -0.3;
  s = 0.7;
}

void solution_eps_0_8( double &alpha, double &beta, double &s ) {
  alpha = 0.18;
  beta = -0.19;
  s = 0.8;
}

void solution_eps_0_9( double &alpha, double &beta, double &s ) {
  alpha = 0.09;
  beta = -0.1;
  s = 0.9;
}

void solution_eps_1( double &alpha, double &beta, double &s ) {
  alpha = 0.001;
  beta = -0.0;
  s = 1;
}

void solution_eps_m_0_1( double &alpha, double &beta, double &s ) {
  alpha = 0.5;
  beta = -0.6;
  s = -0.1;
}

void solution_eps_m_0_075( double &alpha, double &beta, double &s ) {
  alpha = 0.5;
  beta = -0.6;
  s = -0.075;
}

void solution_eps_m_0_05( double &alpha, double &beta, double &s ) {
  alpha = 0.5;
  beta = -0.6;
  s = -0.05;
}

void solution_eps_m_0_15( double &alpha, double &beta, double &s ) {
  alpha = 0.5;
  beta = -0.6;
  s = -0.15;
}
void solution_eps_m_0_1525( double &alpha, double &beta, double &s ) {
  alpha = 0.5;
  beta = -0.6;
  s = -0.1525;
}

void solution_eps_m_0_155( double &alpha, double &beta, double &s ) {
  alpha = 0.5;
  beta = -0.58;
  s = -0.155;
}

void solution_eps_m_0_1575( double &alpha, double &beta, double &s ) {
  alpha = 0.5;
  beta = -0.58;
  s = -0.1575;
}

void solution2_eps_0( double &alpha, double &beta, double &s ) {
//  double alpha = 0.510214, beta = -0.61591;
  alpha = 0.49904295;
  beta = -0.56374644;
  s = 0;
}

void solution2_eps_0_1( double &alpha, double &beta, double &s ) {
//  double alpha = 0.510214, beta = -0.61591;
  alpha = 0.48994137;
  beta = -0.56134736;
  s = -0.1;
}

void solution_eps_0_1575( double &alpha, double &beta, double &s ) {
//  double alpha = 0.510214, beta = -0.61591;
  alpha = 0.47412746;
  beta = -0.58038111;
  s = -0.1575;
}

void solution2_eps_0_1575( double &alpha, double &beta, double &s ) {
//  double alpha = 0.510214, beta = -0.61591;
  alpha = 0.47501806;
  beta = -0.56736644;
  s = -0.1575;
}

void solution2_eps_m_0_159( double &alpha, double &beta, double &s ) {
  alpha = 0.47374;
  beta = -0.578154;
  s = -0.159;
}

void tryuing_to_find_the_second_solution_branch() {
  double alpha = 0, beta = 0;
  double h = 1e-3;
  double s = 1;
  double ksi_0 = 0.;
  double ksi_n = 20;

  std::vector<double> initial_conditions = {0., alpha, 1., beta, 0};

  /*
   * ПОПЫТКА № 1
 * Возьмём значения пристрочных параметров из статьи при том значении s, где две ветки неоднозначного решения
 * больше вссего отличаются (в данном случае s =0.07)
 * */
  ksi_n = 25;
  s = 0.07;
  initial_conditions[1] = 0.49904295;
  initial_conditions[3] = -0.56374644;
  solver Solver1(initial_conditions, s, func, ksi_0, ksi_n, h);
  Solver1.solve_system("./results/second_solution/");


  /*
   * ПОПЫТКА № 2
 * Метод пристрелки по параметрам alpha и s
 * RESULT: выкинул на первую ветку решения (в данном случае подобрал s = 0.2191) с бесконечностью zeta = 21
 * */

  s = 0.07;
  initial_conditions[1] = 0.49904295;
  initial_conditions[3] = -0.56374644;
  solver Solver2(initial_conditions, s, func, ksi_0, ksi_n, h);
  initial_conditions = Solver2.find_infinity_second_branch(initial_conditions, s, "./results/second_solution/", true);

/*
 * ПОПЫТКА № 3
 * Метод пристрелки по параметрам alpha beta
 * RESULT: моментально перекинул на первую ветку решения с zeta =11
 * */
  s = 0.07;
  initial_conditions[1] = 0.49904295;
  initial_conditions[3] = -0.56374644;
  solver Solver3(initial_conditions, s, func, ksi_0, ksi_n, h);
  initial_conditions = Solver3.find_infinity(initial_conditions, s, "./results/second_solution/", true);

/*
 * ПОПЫТКА № 4
 *Попорбуем grid-searchетода пристрелки по параметрам alpha,beta. Будем изменять значения параметров alpha, beta
 *RESULT: выкидывает на первую ветку
 * */
  initial_conditions[1] = 0.49904295;
  initial_conditions[3] = -0.56374644;

  while (initial_conditions[1] > 0.498 && initial_conditions[3] < -0.56) {
    try {
      solver Solver4(initial_conditions, s, func, ksi_0, ksi_n, h);
      Solver4.find_infinity(initial_conditions, s, ".`/results/second_solution/", false);
    } catch (const std::runtime_error &error) {
      std::cout << "not converge" << std::endl;
    }
    initial_conditions[1] -= 0.0001;
    initial_conditions[3] += 0.0001;
  }

  /*
   * ПОПЫТКА № 5
   * IT WORKS
   * Теперь возьмём пристрелочные коэффициенты с первой ветки решения (s=-0.159) и попробуем с маленьким шагом
   * изменять коэффициент beta
   * РЕЗУЛЬТАТ: нашли значение критической точки с точностью до 6-ого знака. -0.160539
   * Нашли часть второй ветви решения до s=-0.1575
   * */

//-0.157224 -0.567054 not converge
  solution2_eps_m_0_159(alpha, beta, s);
  initial_conditions[1] = alpha;
  initial_conditions[3] = beta;
  solver Solver5(initial_conditions, s, func, ksi_0, ksi_n, h);
  double tmp_beta = initial_conditions[3];
  while (tmp_beta < -0.56) {
    initial_conditions[3] = tmp_beta;
    try {
      initial_conditions =
          Solver5.find_infinity_second_branch(initial_conditions, s, "./results/second_solution/", true);
    } catch (const std::runtime_error &error) {
      std::cout << "not converge" << std::endl;
    }
    tmp_beta += 0.0001;
    std::cout << s << " " << tmp_beta << std::endl;
  }
}

void find_first_solution_branch() {
  double alpha = 0, beta = 0;
  double h = 1e-3;
  double s = 1;
  double ksi_0 = 0.;
  double ksi_n = 20;

  solution2_eps_m_0_159(alpha, beta, s);

  std::vector<double> initial_conditions = {0., alpha, 1., beta, 0};

    while (s >= -0.16) {
    if (s <= -0.15) {
      s -= 0.001;
    } else {
      s -= 0.01;
    }
    bool output_f = false;
    if (int(s * 100) % 5 == 0 || ((s < 0.2) && (int(s * 100) % 25 == 0)) || ((s <= -0.15))) {
      output_f = true;
    }
    solver Solver(initial_conditions, s, func, ksi_0, ksi_n, h);
    initial_conditions = Solver.find_infinity(initial_conditions, s, "./results/first_solution/", output_f);
  }
}

// Проходим по всем значениям параметра s и строим линии тока в папке линии тока в директории plots

void streamlines_first_branch() {
  std::string first_path = "./results/first_solution/";

  std::string path = "./results/first_solution/";
  for (auto const &dir_entry : std::filesystem::directory_iterator{first_path}) {
    streamlines::streamlines(1, dir_entry.path(), 15);
//    streamlines::streamlines(sqrt(2), dir_entry.path(), 1);
  }
  system("python3 ./python_scripts/streamlines_plot.py");
}

void streamlines_second_branch() {
  std::string second_path = "./results/second_solution/";

     std::string path  = "./results/second_solution/";
  for (auto const &dir_entry : std::filesystem::directory_iterator{path}) {
    streamlines::streamlines(1, dir_entry.path(), 1);
//    streamlines::streamlines(sqrt(2), dir_entry.path(), 1);
  }
  system("python3 ./python_scripts/streamlines_plot.py");
}

void accuracy (){
  std::string first_path = "./results/first_solution/";
  std::string second_path = "./results/second_solution/";
  std::string epsilon = "eps_-0.1605";
  double alpha = 0, beta = 0;
  double h = 1e-3;
  double s = 1;
  double ksi_0 = 0.;
  double ksi_n = 20;

  std::vector<double> initial_conditions = {0., alpha, 1., beta, 0};

  std::ifstream parameters;

  parameters.open(first_path + epsilon + "/result_parameters.csv");
  if (!parameters.is_open()) {
    throw std::runtime_error("can't open file");
  }
  double zeta_n;
  int n_cells;
  double tmp;
  if (!(parameters >> zeta_n >> h >> s >> alpha >> beta)) {
    throw std::runtime_error("wrong data 1");
  }
  parameters.close();
  h = 1e-4;

  solver Solver1(initial_conditions, s, func, ksi_0, ksi_n, h);
  Solver1.find_infinity(initial_conditions, s, first_path, true);

  /*  parameters.open(second_path + epsilon + "/result_parameters.csv");
  if (!parameters.is_open()) {
    throw std::runtime_error("can't open file");
  }
  if (!(parameters >> zeta_n >> h >> s >> alpha >> beta)) {
    throw std::runtime_error("wrong data 1");
  }
  parameters.close();
  h = 1e-4;

  solver Solver2(initial_conditions, s, func, ksi_0, ksi_n, h);
  Solver2.find_infinity_second_branch(initial_conditions, s, second_path, true);*/
}

int main() {

  double alpha = 0, beta = 0;
  double h = 1e-3;
  double s = 1;
  double ksi_0 = 0.;
  double ksi_n = 20;

  std::vector<double> initial_conditions = {0., alpha, 1., beta, 0};

  std::vector<double> result = func(ksi_0, initial_conditions, s);

  std::string first_path = "./results/first_solution/";
  std::string second_path = "./results/second_solution/";
  std::string epsilon = "eps_-0.1605";

  streamlines_first_branch();

//  streamlines::streamlines(1, first_path + epsilon, 1);
//  streamlines::streamlines(1, second_path + epsilon, 30);

  return 0;

}
