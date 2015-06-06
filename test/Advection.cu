
#include "../Advection.hpp"

using namespace Kolmogorov2D;

const int Nx = 128, Ny = 128;

int main(int argc, char const* argv[]) {
  
  Coefficient<float> C(Nx, Ny);
  C.set(0, 1, {0, 1});
  dx<float>(C, 1.0);

  return 0;
}
