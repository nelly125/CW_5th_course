#include <vector>
#include <cmath>
#include <fstream>

using vector = std::vector< double >;

#define N 5000
#define M 5000

class solver
{
public:
  solver();

private:
  double L = 4, T = 2. * M_PI;

  std::vector< vector > Phi{ N + 1, vector( M + 1 ) };
//  std::vector< vector > rho{ N + 1, vector( M + 1 ) };

  vector xi, eta;

  double dr = L / (N + 1);
  double dt = T / (M + 1);

  int    steps = 100000;
  double Mach  = 0;

public:
  void relax();

  void fully_explicit();

private:
  void foutput_Phi( const std::string &filename );

  void output_c_p( const std::string &filename );

//    [[nodiscard]] std::string generate_filename() const {return std::to_string(Mach) + "_" + std::to_string(steps);}
};

solver::solver()
{
  xi.resize( N + 1 );
  eta.resize( M + 1 );
  for ( int i = 0; i < N + 1; ++i ) { xi[i] = i * dr; }
  for ( int j = 0; j < M + 1; ++j ) { eta[j] = j * dt; }

  // Initial conditions
  for ( int i = 0; i < N + 1; ++i )
  {
    for ( int j = 0; j < M + 1; ++j )
    {
      Phi[i][j] = std::exp( xi[i] ) * std::sin( eta[j] );
    }
  }

  // Boundary conditions
  for ( int j = 0; j < M; ++j )
  {
    Phi[0][j] = 0;   // \xi = 0
    Phi[N][j] = std::exp( xi[N] ) * std::sin( eta[j] ); // \xi = N
  }

  for ( int i = 0; i < N; ++i )
  {
    Phi[i][0] = 0; // \eta = 0
    Phi[i][M] = 0; // \eta = 2pi
  }
}

void solver::foutput_Phi( const std::string &filename )
{
  std::ofstream out( filename + ".txt" );
  std::ofstream r_out( filename + "_r.txt" );
  std::ofstream theta_out( filename + "_theta.txt" );

  for ( int i = 0; i < N + 1; ++i )
  {
    for ( int j = 0; j < M + 1; ++j )
    {
      out << Phi[i][j] << " ";
    }
    out << "\n";
  }
  out.close();

  for ( int i = 0; i < N + 1; ++i ) { r_out << xi[i] << " "; }
  for ( int j = 0; j < M + 1; ++j ) { theta_out << eta[j] << " "; }
}

void solver::output_c_p( const std::string &filename )
{
  double        dtheta = eta[1] - eta[0];
  std::ofstream out_cp( filename );
  for ( int     j      = 0; j < M + 1; ++j )
  {
    double c_p_j = 1. - std::pow((Phi[1][j] - Phi[0][j]) / dtheta, 2. );
    out_cp << eta[j] << " " << c_p_j << std::endl;
  }
  out_cp.close();
}

void solver::fully_explicit()
{
  vector A{ M + 1 }, B{ M + 1 }, C{ M + 1 }, F{ M + 1 };
  vector alpha{ M + 1 }, beta{ M + 1 };

  for ( int k = 0; k < steps; ++k )
  {
    for ( int i = 1; i < N; ++i )
    {
      for ( int j = 1; j < M; ++j )
      {
        Phi[i][j] = (Phi[i + 1][j] + Phi[i - 1][j] + Phi[i][j + 1] + Phi[i][j - 1]) / 4.;
      }
    }
  }
  foutput_Phi( "fully_explicit_");
}

void solver::relax()
{
  vector A( M + 1 ), B(M + 1 ), C( M + 1 ), F(M + 1);

  vector alpha(M + 1), beta(M + 1);

  for ( int k = 0; k < steps; ++k )
  {
    A[0] = 1.;
    B[0] = -4.;
    C[0] = 0;
    F[0] = 0;

    alpha[0] = - C[0] / B[0];
    beta[0] = F[0] / B[0];

    for ( int i = 1; i < N; ++i )
    {
      for ( int j = 1; j < M; ++j )
      {
        A[j] = 1;
        B[j] = -4;
        C[j] = 1;
        F[j] = -Phi[i - 1][j] - Phi[i + 1][j];

        alpha[j] = - C[j] / (B[j] + A[j] * alpha[j - 1]);
        beta[j] = (F[j] - A[j] * beta[j - 1]) / (B[j] + A[j] * alpha[j - 1]);
      }
      for (int j = 1; j < M; ++j) {
        Phi[i][j] = alpha[j] * Phi[i][j + 1] + beta[j];
      }
    }
  }
  foutput_Phi( "relax");
}

int main()
{
/*  {
    solver s;
    s.fully_explicit();
  }*/

//  solver s;
//  s.fully_explicit();

  solver s;
  s.relax();
}