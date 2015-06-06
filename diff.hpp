#pragma once

#include "Kolmogorov2D.hpp"

#include <thrust/transform.h>
#include <thrust/complex.h>

namespace Kolmogorov2D {

template <typename T> __device__ void dx(Coefficient<T> &u, T Lx) {
  int stride = u.get_stride();
  auto b = u.get();
  auto e = u.get() + u.size();
  thrust::transform(b, e, thrust::counting_iterator<int>(0), b,
                    [stride, Lx](thrust::complex<T> a, int idx) {
                      int i = idx / stride;
                      thrust::complex<T> kx(0, 2 * M_PI * i / Lx);
                      return kx * a;
                    });
}

} // namespace Kolmogorov2D
