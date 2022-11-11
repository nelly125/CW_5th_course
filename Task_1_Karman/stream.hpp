#ifndef TASK_1_KARMAN_TASK_1_KARMAN_STREAM_HPP_
#define TASK_1_KARMAN_TASK_1_KARMAN_STREAM_HPP_

#include <sstream>
#include <iomanip>
namespace streamlines {
  struct solution {
    double F;
    double G;
    double H;
  };

  struct psc_parameters {
    double r;
    double phi;
    double zeta;

  };

  void streamlines( double r_0, const std::string &path, int count ) {
    std::string line;
    std::string file_name;

    std::ifstream solution;
    std::ifstream parameters;
    solution.open(path + "/result.csv");
    parameters.open(path + "/result_parameters.csv");
    if (!solution.is_open() || !parameters.is_open()) {
      throw std::runtime_error("can't open file");
    }

    double zeta_n;
    double h;
    double s;
    double alpha;
    double beta;
    int n_cells;
    double tmp;
    if (!(parameters >> zeta_n >> h >> s >> alpha >> beta)) {
      throw std::runtime_error("wrong data 1");
    }
    n_cells = (zeta_n) / h;

    std::vector<streamlines::solution> res;
    res.resize(n_cells + 1);
    int i = 0;
    while (getline(solution, line)) {
      std::stringstream ss(line);
      if (!(ss >> tmp >> res[i].F >> tmp >> res[i].G >> tmp >> res[i].H)) {
        throw std::runtime_error("wrong data 2");
      }
      ++i;
    }

    psc_parameters sk0{};
    psc_parameters sk1{};

    int N = count;
    for (int j = 0; j < N; ++j) {
      sk0.r = r_0;
      sk0.phi = 2 * M_PI * j / N;
      sk0.zeta = zeta_n;

      file_name = path + "/streamlines_";
      std::ostringstream temp;
      temp << std::fixed;
      temp << std::setprecision(4);
      temp << r_0;
      file_name += temp.str() + "_";

      file_name += std::to_string(j);
      file_name += ".csv";
      std::ofstream out_stream(file_name);
      out_stream << std::scientific;
      file_name = path + "/rphiz_";
      file_name += temp.str() + "_";
      file_name += std::to_string(j);
      file_name += ".csv";
      std::ofstream out_rphiz_stream(file_name);
      out_rphiz_stream << std::scientific;

      for (i = n_cells - 1; i > 0; --i) {

        sk1.r = sk0.r - sk0.r * 0.5 * (res[i + 1].F / res[i + 1].H + res[i].F / res[i].H) * h;
        sk1.phi = sk0.phi - 0.5 * (res[i + 1].G / res[i + 1].H + res[i].G / res[i].H) * h;
        sk1.zeta = sk0.zeta - h;

        sk0 = sk1;
//        out_stream << sk0.r << " " << sk0.phi << " " << sk0.zeta << std::endl;
        if (i % 1 == 0) {
          out_stream << sk0.r * cos(sk0.phi) << "\t" << sk0.r * sin(sk0.phi) << "\t" << sk0.zeta << std::endl;
          out_rphiz_stream << sk0.r << "\t" << sk0.phi << "\t" << sk0.zeta << "\t" << res[i].F << "\t" << res[i].G
                           << "\t" << res[i].H << std::endl;
        }

      }
    }
  }
}
#endif //TASK_1_KARMAN_TASK_1_KARMAN_STREAM_HPP_

