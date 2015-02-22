
#include "Converter.hpp"

namespace Kolmogorov2D {

template <> ConverterC2R<float>::ConverterC2R(int Nx, int Ny) {
  cujak::exec(cufftPlan2d(&plan, Nx, Ny, CUFFT_C2R));
}
template <> ConverterC2R<double>::ConverterC2R(int Nx, int Ny) {
  cujak::exec(cufftPlan2d(&plan, Nx, Ny, CUFFT_Z2D));
}

template <>
void ConverterC2R<float>::operator()(const Complex *uf, Real *u) const {
  cujak::exec(cufftExecC2R(plan, const_cast<Complex *>(uf), u));
}

template <>
void ConverterC2R<double>::operator()(const Complex *uf, Real *u) const {
  cujak::exec(cufftExecZ2D(plan, const_cast<Complex *>(uf), u));
}

template <> ConverterR2C<float>::ConverterR2C(int Nx, int Ny) {
  cujak::exec(cufftPlan2d(&plan, Nx, Ny, CUFFT_R2C));
}
template <> ConverterR2C<double>::ConverterR2C(int Nx, int Ny) {
  cujak::exec(cufftPlan2d(&plan, Nx, Ny, CUFFT_D2Z));
}

template <>
void ConverterR2C<float>::operator()(const Real*u, Complex *uf) const {
  cujak::exec(cufftExecR2C(plan, const_cast<Real *>(u), uf));
}

template <>
void ConverterR2C<double>::operator()(const Real *u, Complex *uf) const {
  cujak::exec(cufftExecD2Z(plan, const_cast<Real *>(u), uf));
}

} // namespace Kolmogorov2D
