#pragma once

#include "cujak/cufft.hpp"
#include <thrust/device_ptr.h>

namespace Kolmogorov2D {

/*!
 * @class Coefficient
 * @headerfile Kolmogorov2D.hpp "Kolmogorov2D.hpp"
 *
 * @brief ポインタをFourier係数の列として解釈する
 */
template <typename Float = float> class Coefficient {
public:
  using typename cujak::fft2d::traits<Float>::Complex;
  typedef thrust::device_ptr<Complex> pComplex;

  Coefficient(int Nx_, int Ny_, pComplex u_)
      : Nx(Nx_), Ny(Ny_), stride(Ny / 2 + 1), N(Nx * stride), u(u_) {}

  /* accesors */
  Complex &operator()(int i, int j) { return u[stride * i + j]; }
  Complex get(int i, int j) const { return u[stride * i + j]; }
  void set(int i, int j, Complex v) { u[stride * i + j] = v; }

  pComplex get() const { return u; }

private:
  const int Nx, Ny, stride, N /** Complexとしてのuの個数 */;
  pComplex u;
};

/*!
 * @class Field
 *
 * @brief 実空間の場を保持する
 */
template <typename Float = float> class Field {
public:
  using typename cujak::fft2d::traits<Float>::Real;
  typedef thrust::device_ptr<Real> pReal;

  Field(int Nx_, int Ny_, pReal u_) : Nx(Nx_), Ny(Ny_), u(u_) {}

  /* accesors */
  Real &operator()(int i, int j) { return u[Ny * i + j]; }
  Real get(int i, int j) const { return u[Ny * i + j]; }
  void set(int i, int j, Real v) { u[Ny * i + j] = v; }

  pReal get() const { return u; }

private:
  const int Nx, Ny;
  pReal u;
};

} // namespace Kolmogorov2D
