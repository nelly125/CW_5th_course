#include "system_helpers.hpp"

auto mk_dir( double eps, const std::string& dir) -> std::string {
  std::string directory = dir + "eps_";
  int check;

  std::ostringstream streamObj3;
  streamObj3 << std::fixed;
  streamObj3 << std::setprecision(4);
  streamObj3 << eps;

  directory += streamObj3.str();

  std::cout << directory << std::endl;

  check = mkdir(directory.c_str(), 0777);
  if (!check)
    std::cout << "Directory created\n";
  else {
    std::cout << "Unable to create directory\n";
  }
  directory += "/";

  std::string plot_dir = directory + "plots/";

  check = mkdir(plot_dir.c_str(), 0777);
  if (!check)
    std::cout << "Directory created\n";
  else {
    std::cout << "Unable to create directory\n";
  }

  return directory;
}