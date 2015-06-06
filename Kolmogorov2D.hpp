#pragma once

#include "cujak/cufft.hpp"

#include <thrust/device_vector.h>

namespace Kolmogorov2D {

/*!
 * @class Coefficient
 * @headerfile Kolmogorov2D.hpp "Kolmogorov2D.hpp"
 */
template <typename Float> class Coefficient {
public:
  typedef typename cujak::traits<Float>::Complex Complex;
  typedef thrust::device_vector<Complex> cVector;
  typedef thrust::device_ptr<Complex> iterator;
  typedef Complex value_type;

  Coefficient(int Nx_, int Ny_)
      : Nx(Nx_), Ny(Ny_), stride(Ny / 2 + 1), N(Nx * stride), u(N) {}

  /* accesors */
  Complex get(int i, int j) const { return u[stride * i + j]; }
  void set(int i, int j, Complex v) { u[stride * i + j] = v; }

  thrust::device_ptr<Complex> get() { return u.data(); }
  const thrust::device_ptr<Complex> get() const { return u.data().get(); }

  int size_x() const { return Nx; }
  int size_y() const { return Ny; }
  int size() const { return N; }
  int get_stride() const { return stride; }

private:
  const int Nx, Ny, stride, N /** Complexとしてのuの個数 */;
  cVector u;
};

/*!
 * @class Field
 * @headerfile Kolmogorov2D.hpp "Kolmogorov2D.hpp"
 */
template <typename Float> class Field {
public:
  typedef typename cujak::traits<Float>::Real Real;
  typedef thrust::device_vector<Real> Vector;
  typedef thrust::device_ptr<Real> iterator;
  typedef Real value_type;

  Field(int Nx_, int Ny_) : Nx(Nx_), Ny(Ny_), u(Nx * Ny) {}

  /* accesors */
  Real get(int i, int j) const { return u[Ny * i + j]; }
  void set(int i, int j, Real v) { u[Ny * i + j] = v; }

  Real *get() { return u.data().get(); }
  const Real *get() const { return u.data().get(); }

  int size_x() const { return Nx; }
  int size_y() const { return Ny; }

  iterator begin() { return u.begin(); }
  iterator end() { return u.end(); }
  const iterator begin() const { return u.begin(); }
  const iterator end() const { return u.end(); }

private:
  const int Nx, Ny;
  Vector u;
};

} // namespace Kolmogorov2D
