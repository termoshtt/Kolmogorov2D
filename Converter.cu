
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

} // namespace Kolmogorov2D
