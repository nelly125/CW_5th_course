#ifndef TASK_3_SOLVER_HPP
#define TASK_3_SOLVER_HPP

#include "configure.hpp"

#include <cmath>
#include <vector>
#include <string>

using vector = std::vector< double >;
using matrix = std::vector< vector >;


class solver_r_theta
{
public:
  solver_r_theta();

  virtual ~solver_r_theta() = default;

private:
  double T = M_PI;

  matrix rho{ N + 1, vector( M + 1 ) };

  matrix u{ N + 1, vector( M + 1 ) },
         v{ N + 1, vector( M + 1 ) };

  double dr     = (double) R / (N);
  double dtheta = T / (M);

  double cu = 0.5 / dr;
  vector cv;

  double Mach  = MACH;
  double sigma = 0.1;

public:
  void solve();

private:
  virtual void scheme() = 0;

  virtual void coeff() = 0;

  void foutput_Phi( const std::string &filename );

  [[nodiscard]] std::string generate_filename( const std::string &method, size_t step ) const
  {
    return "results/" + method + "_" +
           std::to_string( N ) + "_" +
           std::to_string( M ) + "_" +
           std::to_string( Mach ) + "_" +
           std::to_string( step );
  }

protected:
  void uvp();

  matrix p{ N + 1, vector( M + 1 ) };
  matrix pPhi{ N + 1, vector( M + 1 ) };
  matrix Phi{ N + 1, vector( M + 1 ) };
  vector c;
  vector a;
  vector r;
  vector theta;

  std::string method_name;

  void updateBoundaries();
};

class relax : public solver_r_theta
{
public:
  relax() { method_name = "relax"; }

private:
  void scheme() override;

  void coeff() override;

  matrix B{ N + 1, vector( M + 1 ) },
         D{ N + 1, vector( M + 1 ) },
         oE{ N + 1, vector( M + 1 ) },
         F{ N + 1, vector( M + 1 ) },
         H{ N + 1, vector( M + 1 ) };

};

class chess : public solver_r_theta
{
public:
  chess() { method_name = "chess"; }

  explicit chess( double s ) : sigma( s ) { method_name = "sigma"; }

private:
  void scheme() override;

  void coeff() override;

  matrix B{ N + 1, vector( M + 1 ) },
         D{ N + 1, vector( M + 1 ) },
         F{ N + 1, vector( M + 1 ) },
         H{ N + 1, vector( M + 1 ) },
         sm{ N + 1, vector( M + 1 ) },
         osp{ N + 1, vector( M + 1 ) };

  double sigma = 0.1;
};

class FEM : public solver_r_theta
{
public:
  FEM() { method_name = "FEM"; }

private:
  void scheme() override;

  void coeff() override;

  matrix B{ N + 1, vector( M + 1 ) },
         D{ N + 1, vector( M + 1 ) },
         oE{ N + 1, vector( M + 1 ) },
         F{ N + 1, vector( M + 1 ) },
         H{ N + 1, vector( M + 1 ) };
};

class Beloglazkin : public solver_r_theta
{
public:

  Beloglazkin() { method_name = "beloglazkin"; }

private:
  void scheme() override;

  void coeff() override;

  matrix B{ N + 1, vector( M + 1 ) },
         D{ N + 1, vector( M + 1 ) },
         E{ N + 1, vector( M + 1 ) },
         F{ N + 1, vector( M + 1 ) },
         H{ N + 1, vector( M + 1 ) };

  matrix b{ N + 1, vector( M + 1 ) },
         c{ N + 1, vector( M + 1 ) },
         d{ N + 1, vector( M + 1 ) },
         e{ N + 1, vector( M + 1 ) },
         f{ N + 1, vector( M + 1 ) };

};

class sweep : public solver_r_theta
{
public:
  sweep() { method_name = "sweep"; }

private:
  void scheme() override;

  void coeff() override;

  matrix B{ N + 1, vector( M + 1 ) },
         D{ N + 1, vector( M + 1 ) },
         E{ N + 1, vector( M + 1 ) },
         F{ N + 1, vector( M + 1 ) },
         H{ N + 1, vector( M + 1 ) };
};

#endif //TASK_3_SOLVER_HPP