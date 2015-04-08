#pragma once

#include "Kolmogorov2D.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>

namespace Kolmogorov2D {

template <class Data>
void output_ascii(const Data &data, std::string filename) {
  std::ofstream ofs(filename.c_str());
  ofs << std::scientific << std::setprecision(7);
  output_ascii(data, ofs);
}

template <typename T>
void output_ascii(const Coefficient<T> &C, std::ostream &ost) {
  int Nx = C.size_x();
  int stride = C.get_stride();
  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < stride; j++) {
      auto c = C.get(i, j);
      ost << i << " " << j << " " << c.x << " " << c.y << "\n";
    }
    ost << '\n';
  }
  ost << std::flush;
}

template <typename T> void output_ascii(const Field<T> &F, std::ostream &ost) {
  int Nx = F.size_x();
  int Ny = F.size_y();
  for (int i = 0; i < Nx; i++) {
    for (int j = 0; j < Ny; j++) {
      ost << i << " " << j << " " << F.get(i, j) << "\n";
    }
    ost << '\n';
  }
  ost << std::flush;
}

} // namespace Kolmogorov2D
