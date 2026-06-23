# Tnum-Tana


Tnum-Tana is a Fortran library which deals with kinetic energy operators (KEO) and coordinates transformations.
More precisely, this library enables one to compute :

- the KEO or the metric tensor (covariant, g, or contravariant, G, components) using curvilinear coordinates for a molecular system of any size, without build in limitation.
- Coordinates transformations from curvilinear coordinates to Cartesian coordinates or the reverse (when possible). 

1. The **Tnum** code enables to calculate numerically and exactly the KEO, g and G. With respect to similar approaches or codes, **Tnum** is highly flexible in terms coordinates. One can add several coordinates transformations (no limitation) between the coordinates used for the dynamics (Qdyn or Qact) and the Cartesian coordinates (body-fixed).
The exactness is guaranteed by the automatic differentiation procedures used for all coordinate transformations (except one or two).
2. The **Tana** code enables to compute the analitical expression of the KEO and G with polysherical coordinates. It has no limitation in terms of the molecular size of the number of dregrees of freedom.

**Tnum** and **Tana** share the same imput file.

## 1) Installation

### 1a) Installation with a makefile

The installation is simple with a makefile. However, we do not have an fully automatic procedure (like configure ...). The program uses some Fortran 2003 features. Therefore, the compilers gfortran need to be recent.
The code works with the following compilers:

- gfortran-11 and above
- ifx 
- nagfor (nAG compiler version 7.1)

For instance, main executables can be built:

```bash
make app
# Main Tnum/Tana:                         Tnum90.exe
# Main Tnum/Tana for MCTDH and Quantics:  Tnum90_MCTDH.exe
# Main Tnum/Tana executable for MidasCpp: Tnum90_MidasCpp.exe
#
# Example of Fortran driver:              Main_TnumTana_FDriver.exe 
# Example of c driver:                    Main_TnumTana_cDriver.exe
#
# Example of Fortran driver (cartesian to curvilinear transformation): Main_X2Q.exe
```

The libraries can be built:

```bash
make lib
```

The libraries with the name starting with **libTnum-Tana** contain only the objects of Tnum-Tana.
The libraries with the name starting with **libTnum-TanaFull** contain all objects from Tnum-Tana and from the external libraries.


Furthermore, one can add options:

```bash
make FC=ifort OMP=0 OPT=0 LAPACK=0 INT=4 RKIND=real64
  # FC=ifort to change the compiller to ifort
  # OMP=0/1 to turn off/on the OpenMP fortran flag.
  # OPT=0/1 to turn off/on the fortran optimization.
  # LAPACK=0/1 to turn off/on the lapack use
  # INT=4/8 to change the compiler default integer
  # RKIND=real64/real32/real128 to change the real kind. 
```

The library, **libTnum-Tana_XXX_optY_ompZ_lapackW_intV_realR.a** is created in the main directory.
Remarks: 

- XXX is the compiler (gfortran, ifort ...)
- Y is 0 or 1 (opt0 / opt1: compiler optimization)
- Z is 0 or 1 (omp0 / omp1: whitout/with OpenMP)
- W is 0 or 1 (lapack0 / lapack1: whitout/with lapack)
- V is 4 or 8 (int4 / int8)
- R is 64 or 32 or 128 (eal64 / real32 / real128 

If needed, the .mod files are in the **obj/obj_XXX_optY_ompZ_lapackW_intV_realR** directory.

### 1b) Installation with fpm

**fpm** is a Fortran Package Manager (see https://fpm.fortran-lang.org).
If it is not done yet, the external library directory links (in Ext_lib directory) need to be set by:

```bash
make fpm
```

Then to build 

```bash
fpm build
```

To run some tests:

```bash
fpm run cDriver --< TESTS/exa_TnumDriver/dat_driver0
fpm run cDriver --< TESTS/exa_TnumDriver/dat_driver1

fpm run X2Q --< TESTS/exa_TnumDriver/dat_driver0
fpm run X2Q --< TESTS/exa_TnumDriver/dat_driver1

fpm run FDriver -- -i TESTS/exa_TnumDriver/dat_driver0
fpm run FDriver -- -i TESTS/exa_TnumDriver/dat_driver1

fpm run Tnum90          --< TESTS/exa_TnumDriver/dat_driver0
fpm run Tnum90_MidasCpp --< TESTS/exa_TnumDriver/dat_driver0
fpm run Tnum90_MCTDH    --< TESTS/exa_TnumDriver/dat_driver0
```

## 2) Documentation

See [Tnum-Tana](https://github.com/lauvergn/Tnum-Tana/wiki) Wiki (uncomplete)

## 3) Some driver examples

Main programs are in the APP directory.

### 3a) Eckart transformation

The Fortran main is "APP/Main_Eckart.f90".
To compile the code :

```
make OPT=0 app
```
or
```
make OPT=0 app APPSRC=Main_Eckart.f90
```

The code reads two geometries: the reference one and then the current one. 
Remark: If masses are read (--read_masses t), the list of masses must be given after the reference geometry (the first one).

Then, the code computes the Eckart rotation matrix.
Depending of the argumant options:

- masses can be read (-m or --read_masses). They are read in g.mol$^{-1}$ or unified atomic mass unit (u or AMU) or Dalton (Da).
- the current transformed geometry is printed (-t or --coord_transfo)

Example, without reading the masses. Therefore, Tnum masses are used (CODATA2006 physical constants and NIST2012 masses).

```
./Main_Eckart.exe -i APP_InputFiles/dat_Eckart -o APP_OutputFiles/res --coord_transfo .TRUE.
```
or
```
./Main_Eckart.exe -i APP_InputFiles/dat_Eckart -o APP_OutputFiles/res --coord_transfo .TRUE. --read_masses f
```
The Eckart rotation matrix and the transformed Cartessian geometry are given in "APP_OutputFiles/res"

Example: reading the masses.

```
./Main_Eckart.exe -i APP_InputFiles/dat_Eckart2 -o APP_OutputFiles/res --coord_transfo .TRUE. --read_masses t
```

The Eckart rotation matrix and the transformed Cartessian geometry are given in "APP_OutputFiles/res"


### 3b) Inter-fragment coordinates

The Fortran main is "APP/Main_NMdimer.f90".
To compile the code :

```
make OPT=0 app
```
or
```
make OPT=0 app APPSRC=Main_NMdimer.f90
```

The code reads two fragment geometries (A and B) to form a dimer (AB). 
Remark: If masses are read (--read_masses t), the list of masses must be given after each geometry.

Then, the code computes the coordinates associated to the inter-fragment coordinates.
Depending of the argument options:

- masses can be read (-m or --read_masses). They are read in g.mol$^{-1}$ or unified atomic mass unit (u or AMU) or Dalton (Da).

Example, without reading the masses. Therefore, Tnum masses are used (CODATA2006 physical constants and NIST2012 masses).

```
./Main_NMdimer.exe -i APP_InputFiles/dat_dimer -o APP_OutputFiles/res
```
or
```
./Main_NMdimer.exe -i APP_InputFiles/dat_dimer -o APP_OutputFiles/res --read_masses f
```

The code also generates several coordinate motions (can be visual with jmol as vibrations):

- The translational and rotational motions of fragments, A and B.
- The inter-fragment motions (6 coordinates) of dimer AB.

Example: reading the masses.

```
./Main_NMdimer.exe -i APP_InputFiles/dat_dimer2 -o APP_OutputFiles/res --read_masses f
```

The inter-fragment coordinates of dimer AB are given at the end of "APP_OutputFiles/res" (search *Inter_dimAB*) as 6 linear combinations (in column) of Cartesian coordinates.