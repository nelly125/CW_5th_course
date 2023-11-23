#include <stdexcept>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "solver.hpp"

double matrix_diff( const matrix &a, const matrix &b )
{
  if ( a.size() != b.size()) { return 0; }
  double    diff = 0;
  for ( int i    = 0; i < N + 1; ++i )
  {
    for ( int j = 0; j < M + 1; ++j )
    {
      if ( std::isnan( a[i][j] )) { throw std::runtime_error( "nan" ); }
      diff += fabs( a[i][j] - b[i][j] );
    };
  }
  return diff;
}

solver_r_theta::solver_r_theta()
{
  r.resize( N + 1 );
  theta.resize( M + 1 );
  for ( int i = 0; i < N + 1; ++i ) { r[i] = 1 + i * dr; }
  for ( int j = 0; j < M + 1; ++j ) { theta[j] = j * dtheta; }

  a.resize( N + 1 );
  c.resize( N + 1 );
  cv.resize( N + 1 );

  for ( int i = 0; i < N + 1; ++i )
  {
    a[i]  = dr / (2 * r[i]);
    c[i]  = std::pow( dr / (r[i] * dtheta), 2. );
    cv[i] = 0.5 / (r[i] * dtheta);
  }

  // Initial conditions
  for ( int i = 0; i < N + 1; ++i )
  {
    for ( int j = 0; j < M + 1; ++j ) { Phi[i][j] = (r[i]) * std::cos( theta[j] ); }
  }
}

void solver_r_theta::uvp()
{
  double cMM = Mach * Mach * 0.2;

  for ( int j = 1; j < M; ++j )
  {
    for ( int i = 1; i < N; ++i )
    {
      u[i][j] = (Phi[i + 1][j] - Phi[i - 1][j]) * cu;
      v[i][j] = (Phi[i][j + 1] - Phi[i][j - 1]) * cv[i];
      double V = std::pow( u[i][j], 2. ) + std::pow( v[i][j], 2. );
      double q = (1 - V) * cMM + 1;
      p[i][j] = std::sqrt( q ) * std::pow( q, 2. );
    }

    u[0][j] = u[1][j];
    u[N][j] = u[N - 1][j];
    v[0][j] = v[1][j];
    v[N][j] = v[N - 1][j];
    p[0][j] = p[1][j];
    p[N][j] = p[N - 1][j];
  }

  for ( int i = 0; i < N; ++i )
  {
    u[i][0] = u[i][1];
    u[i][M] = u[i][M - 1];
    v[i][0] = v[i][1];
    v[i][M] = v[i][M - 1];
    p[i][0] = p[i][1];
    p[i][M] = p[i][M - 1];
  }
}

void solver_r_theta::updateBoundaries()
{
  for ( int j = 0; j < M + 1; ++j )
  {
    Phi[0][j] = Phi[1][j];
    Phi[N][j] = r[N] * std::cos( theta[j] );
  }
  for ( int i = 0; i < N + 1; ++i )
  {
    Phi[i][0] = Phi[i][1];
    Phi[i][M] = Phi[i][M - 1];
  }
}

void solver_r_theta::solve()
{
  int step;
  for ( step = 0; step < 1000000; ++step )
  {
    if ( step % 1000 == 0 ) { std::cout << step << std::endl; }
    pPhi = Phi;
    coeff();
    scheme();
    if ( std::fabs( matrix_diff( Phi, pPhi )) < 1e-7 ) { break; }
  }
  std::cout << step << std::endl;
  foutput_Phi( generate_filename( method_name, step ));
}


void solver_r_theta::foutput_Phi( const std::string &filename )
{
  std::ofstream out( filename + ".txt" );
  std::ofstream r_out( filename + "_r.txt" );
  std::ofstream theta_out( filename + "_theta.txt" );

  out << std::setprecision( 6 );

  for ( int i = 0; i < N + 1; ++i )
  {
    for ( int j = 0; j < M + 1; ++j )
    {
      out << Phi[i][j] << " ";
    }
    out << "\n";
  }
  out.close();

  for ( int i = 0; i < N + 1; ++i ) { r_out << r[i] << " "; }
  for ( int j = 0; j < M + 1; ++j ) { theta_out << theta[j] << " "; }
}

void relax::coeff()
{
  uvp();
  for ( int j = 1; j < M; ++j )
  {
    for ( int i = 1; i < N; ++i )
    {
      B[i][j] = (p[i][j - 1] + p[i][j]) * 0.5 * c[i];
      H[i][j] = (p[i][j] + p[i][j + 1]) * 0.5 * c[i];
      double Pim = (p[i - 1][j] + p[i][j]) * 0.5;
      double Pip = (p[i][j] + p[i + 1][j]) * 0.5;
      double Eij = Pip + Pim + B[i][j] + H[i][j];
      double aP  = a[i] * p[i][j];
      F[i][j]  = Pip + aP;
      D[i][j]  = Pim - aP;
      oE[i][j] = 1 / Eij;
    }
  }
}

void relax::scheme()
{
  auto calculatePhi = [this]( int i, int j ) {
      return (Phi[i][j - 1] * B[i][j] +
              Phi[i - 1][j] * D[i][j] +
              Phi[i + 1][j] * F[i][j] +
              Phi[i][j + 1] * H[i][j]) * oE[i][j];
  };

  for ( int i = 1; i < N; ++i )
  {
    for ( int j = 1; j < M; ++j ) { Phi[i][j] = calculatePhi( i, j ); }
  }
  for ( int i = N - 1; i > 0; --i )
  {
    for ( int j = M - 1; j > 0; --j ) { Phi[i][j] = calculatePhi( i, j ); }
  }
  updateBoundaries();

  for ( int j = M - 1; j > 0; --j )
  {
    for ( int i = 1; i < N; ++i ) { Phi[i][j] = calculatePhi( i, j ); }
  }
  for ( int j = 1; j < M; ++j )
  {
    for ( int i = N - 1; i > 0; --i ) { Phi[i][j] = calculatePhi( i, j ); }
  }
  updateBoundaries();
}


void chess::coeff()
{
  uvp();
  for ( int j = 1; j < M; ++j )
  {
    for ( int i = 1; i < N; ++i )
    {
      B[i][j] = (p[i][j - 1] + p[i][j]) * 0.5 * c[i];
      H[i][j] = (p[i][j] + p[i][j + 1]) * 0.5 * c[i];
      double Pim = (p[i - 1][j] + p[i][j]) * 0.5;
      double Pip = (p[i][j] + p[i + 1][j]) * 0.5;
      double Eij = Pip + Pim + B[i][j] + H[i][j];
      double aP  = a[i] * p[i][j];
      F[i][j]   = Pip + aP;
      D[i][j]   = Pim - aP;
      sm[i][j]  = sigma - Eij / 2;
      osp[i][j] = 1. / (sigma + Eij / 2);
    }
  }
}

void chess::scheme()
{
  auto calculatePhi = [this]( int i, int j ) {
      return (Phi[i][j] * sm[i][j] +
              Phi[i][j - 1] * B[i][j] +
              Phi[i - 1][j] * D[i][j] +
              Phi[i + 1][j] * F[i][j] +
              Phi[i][j + 1] * H[i][j]) * osp[i][j];
  };

  for ( int i = 1; i < N; ++i )
  {
    for ( int j = 1; j < M; ++j )
    {
      if ((i + j) % 2 == 0 ) { Phi[i][j] = calculatePhi( i, j ); }
    }
  }
  for ( int i = 1; i < N; ++i )
  {
    for ( int j = 1; j < M; ++j )
    {
      if ((i + j) % 2 == 1 ) { Phi[i][j] = calculatePhi( i, j ); }

    }
  }
  updateBoundaries();
}

void FEM::coeff()
{
  uvp();
  for ( int j = 1; j < M; ++j )
  {
    for ( int i = 1; i < N; ++i )
    {
      B[i][j] = (p[i][j - 1] + p[i][j]) * 0.5 * c[i];
      H[i][j] = (p[i][j] + p[i][j + 1]) * 0.5 * c[i];
      double Pim = (p[i - 1][j] + p[i][j]) * 0.5;
      double Pip = (p[i][j] + p[i + 1][j]) * 0.5;
      double Eij = Pip + Pim + B[i][j] + H[i][j];
      double aP  = a[i] * p[i][j];
      F[i][j]  = Pip + aP;
      D[i][j]  = Pim - aP;
      oE[i][j] = 1. / Eij;
    }
  }
}


void FEM::scheme()
{
  for ( int i = 1; i < N; ++i )
  {
    for ( int j = 1; j < M; ++j )
    {
      Phi[i][j] = (Phi[i][j - 1] * B[i][j] +
                   Phi[i - 1][j] * D[i][j] +
                   Phi[i + 1][j] * F[i][j] +
                   Phi[i][j + 1] * H[i][j]) * oE[i][j];
    }

  }
  updateBoundaries();
}

void sweep::coeff()
{
  uvp();
  for ( int j = 1; j < M; ++j )
  {
    for ( int i = 1; i < N; ++i )
    {
      B[i][j] = (p[i][j - 1] + p[i][j]) * 0.5 * c[i];
      H[i][j] = (p[i][j] + p[i][j + 1]) * 0.5 * c[i];
      double Pim = (p[i - 1][j] + p[i][j]) * 0.5;
      double Pip = (p[i][j] + p[i + 1][j]) * 0.5;
      double Eij = Pip + Pim + B[i][j] + H[i][j];
      double aP  = a[i] * p[i][j];
      F[i][j] = Pip + aP;
      D[i][j] = Pim - aP;
      E[i][j] = Eij;
    }
  }
}

void sweep::scheme()
{
  vector c_A( M + 1 ), c_B( M + 1 ), c_C( M + 1 ), c_F( M + 1 );
  vector alpha( M + 1 ), beta( M + 1 );

  for ( int j = 1; j < M; ++j )
  {
    c_A[0] = 0;
    c_C[0] = (p[1][j] + p[0][j]) / 2 + a[0] * p[0][j];
    c_B[0] = -(p[1][j] + p[0][j]) / 2 - (p[0][j]) -
             c[0] * ((p[0][j + 1] + p[0][j]) / 2 + (p[0][j] + p[0][j - 1]) / 2);
    c_F[0] = 0;

    alpha[1] = -c_C[0] / c_B[0];
    beta[1]  = c_F[0] / c_B[0];
    for ( int i = 1; i < N; ++i )
    {
      c_A[i] = D[i][j];
      c_B[i] = -E[i][j];
      c_C[i] = F[i][j];
      c_F[i] = -B[i][j] * Phi[i][j - 1] - H[i][j] * Phi[i][j + 1];

      alpha[i + 1] = -c_C[i] / (c_B[i] + c_A[i] * alpha[i]);
      beta[i + 1]  = (c_F[i] - c_A[i] * beta[i]) / (c_B[i] + c_A[i] * alpha[i]);
    }
    for ( int i = N - 1; i > 0; --i )
    {
      Phi[i][j] = alpha[i + 1] * Phi[i + 1][j] + beta[i + 1];
    }
  }
  updateBoundaries();
}

void Beloglazkin::coeff() {
  uvp();
  for ( int j = 1; j < M; ++j )
  {
    for ( int i = 1; i < N; ++i )
    {
      B[i][j] = (p[i][j - 1] + p[i][j]) * 0.5 * c[i];
      H[i][j] = (p[i][j] + p[i][j + 1]) * 0.5 * c[i];
      double Pim = (p[i - 1][j] + p[i][j]) * 0.5;
      double Pip = (p[i][j] + p[i + 1][j]) * 0.5;
      double Eij = Pip + Pim + B[i][j] + H[i][j];
      double aP  = a[i] * p[i][j];
      F[i][j] = Pip + aP;
      D[i][j] = Pim - aP;
      E[i][j] = Eij;
    }
  }


}

void Beloglazkin::scheme() {

}

