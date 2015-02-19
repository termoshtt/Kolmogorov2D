#pragma once

#include "cujak/cufft.hpp"
#include <thrust/device_vector.h>

namespace Kolmogorov2D {

/*!
 * @class Coefficient
 *
 * @brief Fourier係数として解釈する
 */
template <typename Float = float, typename Int = unsigned int>
class Coefficient {
public:
  using typename cujak::fft2d::traits<Float>::Real;
  using typename cujak::fft2d::traits<Float>::Complex;

private:
  const Int Nx, Ny, stride, N /** Complexとしてのuの個数 */;
  Complex *u;

public:
  Coefficient(Int Nx_, Int Ny_, Complex *u_)
      : Nx(Nx_), Ny(Ny_), stride(Ny / 2 + 1), N(Nx * stride), u(u_) {}

  /* accesors */
  Complex &operator()(Int i, Int j) { return u[stride * i + j]; }
  Complex get(Int i, Int j) const { return u[stride * i + j]; }
  void set(Int i, Int j, Complex v) { u[stride * i + j] = v; }
};

/*!
 * @class Field
 *
 * @brief 実空間の場を保持する
 */
template <typename Float = float, typename Int = unsigned int> class Field {
public:
  Field(Int Nx, Int Ny);
};

/*!
 * @class ConverterC2R
 *
 * @brief Fourier係数から場を計算する
 */
class ConverterC2R {};

/*!
 * @class ConverterR2C
 *
 * @brief 場から係数を計算する
 */
class ConverterR2C {};

} // namespace Kolmogorov2D
