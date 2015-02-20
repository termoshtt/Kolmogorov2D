
#include "../Converter.hpp"
#include <iostream>

typedef thrust::device_vector<float> Vector;
using Kolmogorov2D::Coefficient;
using Kolmogorov2D::Field;

int main(int argc, char const *argv[]) {
  const int Nx = 128, Ny = 128;
  Vector u(Nx * Ny);
  Field<float> F(Nx, Ny, u);

  F.set(0, 1, 1.0);
  std::cout << F.get(0, 1) << std::endl;

  return 0;
}
