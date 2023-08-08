# Tnum-Tana


Tnum-Tana is a Fortran library which deals with kinetic energy operators (KEO) and coordinates transformations.
Maore precisely, this library enables one to compute :

- the KEO or the metric tensor (covariant, g, or contravariant, G, components) using curvilinear coordinates for a molecular system of any size, without build in limitation.
- Coordinates transformations from curvilinear coordinates to Cartesian coordinates or the reverse (when possible). 

1. The **Tnum** code enables to calculate numerically and exactly the KEO, g and G. With respect to similar approaches or codes, it is highly flexible in terms coordinates. One can add several coordinates transformations (no limitation) between the coordinates used for the dynamics (Qdyn or Qact) and the Cartesian coordinates (body-fixed).
The exactness is guaranteed by the automatic differentiation procedures used for all coordinate transformations (except one or two).
2. The **Tana** code enables to compute the analitical expression of the KEO and G with polysherical coordinates. It has no limitation in terms of the molecular size of the number of dregrees of freedom.

**Tnum** and **Tana** share the same imput file.

## 1) Installation

The installation is simple with a makefile. However, we do not have an fully automatic procedure (like configure ...). The program uses some Fortran 2003 features. Therefore, the compilers gfortran or ifort need to be recent.

For instance, main executables can be built:

```bash
# Main Tnum/Tana executable:
make Tnum
# Example of Fortran driver:
make Tnum_FDriver
# Example of c driver:
make Tnum_cDriver
# Special main Tnum/Tana executable for MCTDH:
make Tnum_MCTDH
# Special main Tnum/Tana executable for MidasCpp:
make Tnum_MidasCpp
```

If one to build the TnumTana library:

```bash
make lib
```

Furthermore, one can add options:

```bash
make FC=ifort OMP=0 OPT=0 LAPACK=0 INT=4
  # FC=ifort to change the compiller to ifort
  # OMP=0/1 to turn off/on the OpenMP fortran flag.
  # OPT=0/1 to turn off/on the fortran optimization.
  # LAPACK=0/1 to turn off/on the lapack use
  # INT=4/8 to change the default integer
```

The library, **libTnum-Tana_dnSVM_XXX_oppY_ompZ_lapackW_intV.a** is created in the main directory.
Remarks : 

- XXX is the compiller (gfortran, ifort ...)
- Y is 0 or 1 (opt0 / opt1: compiler optimization)
- Z is 0 or 1 (omp0 / omp1: whitout/with OpenMP)
- W is 0 or 1 (lapack0 / lapack1: whitout/with lapack)
- V is 4 or 8 (int4 / int8)

If needed, the .mod files are in the **obj/obj_XXX_oppY_ompZ_lapackW_intV** directory.
