#pragma once

#include "Kolmogorov2D.hpp"

#include "cujak/exception.hpp"

namespace Kolmogorov2D {

/*!
 * @class ConverterC2R
 * @headerfile Converter.hpp "Converter.hpp"
 * @brief Fourier係数から場を計算する
 *
 * cufftHandleの薄いラッパー。
 *
 * cuFFTのマニュアルによりout-placeで変換する方
 * が推奨されている事に注意し、基本的にout-placeで変換する。
 */
template <typename Float> class ConverterC2R {
public:
  typedef typename cujak::fft2d::traits<Float>::Real Real;
  typedef typename cujak::fft2d::traits<Float>::Complex Complex;

  ConverterC2R(int Nx, int Ny);
  ~ConverterC2R() { cufftDestroy(plan); }
  void operator()(const Complex *uf, Real *u) const;
  void operator()(const Coefficient<Float> &uf, Field<Float> &u) const {
    operator()(uf.get(), u.get());
  }

private:
  cufftHandle plan;
};

/*!
 * @class ConverterR2C
 *
 * @brief 場から係数を計算する
 */
class ConverterR2C {};

} // namespace Kolmogorov2D
