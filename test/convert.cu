
#include "../Converter.hpp"
#include <fstream>
#include <iomanip>

typedef typename cujak::fft2d::traits<float>::Real Real;
typedef typename cujak::fft2d::traits<float>::Complex Complex;
typedef thrust::device_vector<Real> Vector;
typedef thrust::device_vector<Complex> cVector;
typedef thrust::host_vector<Real> hVector;
typedef thrust::host_vector<Complex> hcVector;

using namespace Kolmogorov2D;

const int Nx = 128, Ny = 128;

void output_field(Vector& u, std::string filename){
  hVector h(u);
  std::ofstream ofs1(filename.c_str());
  ofs1 << std::scientific << std::setprecision(7);
  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
      ofs1 << i << " " << j << " " << h[Ny * i + j] << "\n";
    }
    ofs1<< std::endl;
  }
}

void output_coef(cVector& uf, std::string filename){
  hcVector hf(uf);
  std::ofstream ofs2(filename.c_str());
  ofs2 << std::scientific << std::setprecision(7);
  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny/2 + 1; j++) {
      Complex c = hf[Ny * i + j];
      ofs2 << i << " " << j << " " << c.x << " " << c.y << "\n";
    }
    ofs2 << std::endl;
  }
}

int main(int argc, char const *argv[]) {
  cVector uf(Nx * (Ny/2+1));
  Vector u(Nx * Ny);
  Coefficient<float> C(Nx, Ny, uf);
  Field<float> F(Nx, Ny, u);

  ConverterC2R<float> c2r(Nx, Ny);
  ConverterR2C<float> r2c(Nx, Ny);

  Complex c = { 0.0, 1.0 };
  C.set(0, 1, c);
  C.set(1, 0, c);

  output_coef(uf, "coef1.dat");
  c2r(C, F);
  output_field(u, "field.dat");
  r2c(F, C);
  output_coef(uf, "coef2.dat");

  return 0;
}
