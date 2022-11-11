#pragma once

#include <vector>
#include <stdexcept>

using vector = std::vector<double>;

vector operator+( const vector &, double );
vector operator*( double y, const vector &x );
vector operator/( const vector &x, double y );
vector operator+( vector x, vector y );

vector operator+( const vector &x, double y ) {
  vector res(x.size());
  for (uint32_t i = 0; i < x.size(); ++i) {
    res[i] = x[i] + y;
  }
  return res;
}

vector operator*( double y, const vector &x ) {
  vector res(x.size());
  for (uint32_t i = 0; i < x.size(); ++i) {
    res[i] = x[i] * y;
  }
  return res;
}

vector operator/( const vector &x, double y ) {
  vector res(x.size());
  for (uint32_t i = 0; i < x.size(); ++i) {
    res[i] = x[i] / y;
  }
  return res;
}

vector operator+( vector x, vector y ) {
  if (x.size() != y.size()) {
    throw std::runtime_error("(((");
  }
  vector res(x.size());
  for (uint32_t i = 0; i < x.size(); ++i) {
    res[i] = x[i] + y[i];
  }
  return res;
}

double vector_norm( vector x, vector y ) {
  double norm = 0;
  if (x.size() != y.size()) {
    throw std::runtime_error("(((");
  }
  for (uint32_t i = 0; i < x.size(); ++i) {
    norm += fabs(x[i] - y[i]);
  }
  return sqrt(norm);
}