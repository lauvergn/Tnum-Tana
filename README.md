# Tnum-Tana


Tnum-Tana is a Fortran library which deals with kinetic energy operator (KEO), T, or the metric tensor (covariant, g, or contravariant, G, components) using curvilinear coordinates for a molecular system of any size, without build in limitation.

1. The **Tnum** code enables to calculate numerically and exactly the KEO, g and G. With respect to similar approaches or codes, it is highly flexible in terms coordinates. One can add several coordinates transformations (no limitation) between the coordinates used for the dynamics (Qdyn or Qact) and the Cartesian coordinates (body-fixed).
The exactness is guaranteed by the automatic differentiation procedures used for all coordinate transformations (except one or two).
2. The **Tana** code enables to compute the analitical expression of the KEO and G with polysherical coordinates. It has no limitation in terms of the molecular size of the number of dregrees of freedom.


**Tnum** and **Tana** share the same imput file.



## 1) Installation

The installation is simple with a makefile. However, we do not have an fully automatic procedure (like configure ...). The program uses some Fortran 2003 features. Therefore, the compilers gfortran or ifort need to be recent.

