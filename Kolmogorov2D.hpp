#pragma once

#include <thrust/device_vector.h>

namespace Kolmogorov2D {

/*!
 * @class Coefficient
 *
 * @brief Fourier係数を保持する
 */
class Coefficient {
public:
  Coefficient(unsigned int Nx, unsigned int Ny);
};

/*!
 * @class Field
 *
 * @brief 実空間の場を保持する
 */
class Field {
public:
  Field(unsigned int Nx, unsigned int Ny);
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
