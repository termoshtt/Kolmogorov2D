
#include "../Converter.hpp"
#include "../output.hpp"

#include <fstream>
#include <iomanip>

typedef typename cujak::traits<float>::Real Real;
typedef typename cujak::traits<float>::Complex Complex;

using namespace Kolmogorov2D;

const int Nx = 128, Ny = 128;

int main(int argc, char const *argv[]) {
  Coefficient<float> C(Nx, Ny);
  Field<float> F(Nx, Ny);

  ConverterC2R<float> c2r(Nx, Ny);
  ConverterR2C<float> r2c(Nx, Ny);

  Complex c = { 0.0, 1.0 };
  C.set(0, 1, c);
  C.set(1, 0, c);

  output_ascii(C, "coef1.dat");
  c2r(C, F);
  output_ascii(F, "field.dat");
  r2c(F, C);
  output_ascii(C, "coef2.dat");

  return 0;
}
