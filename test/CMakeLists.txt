cmake_minimum_required(VERSION 2.8)

cujak_add_test(convert S convert.cu cuKolmogorov2D)
cujak_add_test(Advection S Advection.cu cuKolmogorov2D)
