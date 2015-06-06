#pragma once

#include "Kolmogorov2D.hpp"
#include "Converter.hpp"
#include "diff.hpp"

namespace Kolmogorov2D {

template <typename T> class Advection {
  Coefficient<T> tmp;
  Field<T> f, g;
  ConverterC2R<T> c2r;
  ConverterR2C<T> r2c;

public:
  Advection(int Nx, int Ny)
      : tmp(Nx, Ny), f(Nx, Ny), g(Nx, Ny), c2r(Nx, Ny), r2c(Nx, Ny) {}

  void operator()(Coefficient<T> &u);
  void operator()(const Coefficient<T> &, Coefficient<T> &);
};

} // namespace Kolmogorov2D
