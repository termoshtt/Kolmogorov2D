
#include "../Converter.hpp"
#include <iostream>

typedef typename cujak::fft2d::traits<float>::Real Real;
typedef typename cujak::fft2d::traits<float>::Complex Complex;
typedef thrust::device_vector<Real> Vector;
typedef thrust::device_vector<Complex> cVector;
typedef thrust::host_vector<Real> hVector;

using namespace Kolmogorov2D;

int main(int argc, char const *argv[]) {
  const int Nx = 128, Ny = 128;

  cVector uf(Nx*Ny);
  Coefficient<float> C(Nx, Ny, uf);
  Complex c = {0.0, 1.0};
  C.set(0, 1, c);

  Vector u(Nx * Ny);
  Field<float> F(Nx, Ny, u);

  ConverterC2R<float> c2r(Nx, Ny);
  c2r(C, F);

  hVector h(u);
  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
      std::cout << i << " " << j << " " << h[Ny*i+j] << "\n";
    }
    std::cout << std::endl;
  }
  return 0;
}
