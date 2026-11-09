# Tnum-Tana


Tnum-Tana is a Fortran library which deals with kinetic energy operators (KEO) and coordinates transformations.
More precisely, this library enables one to compute:

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

To run unitary tests:
```bash
make ut
# The 32 tests should pass.
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

The libraries, **libTnum-Tana_XXX_optY_ompZ_lapackW_intV_realR.a** and **libTnum-TanaFull_XXX_optY_ompZ_lapackW_intV_realR.a** are created in the main directory.
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


The available Fortran drivers are (in SRC/TnumTana_Driver.f90):

```fortran
CALL Init_TnumTana_Driver(nb_act,nb_cart,init_sub)
  ! This subroutine reads the input file and initializes the coordinates and the reference geometry.
  !   nb_act:   is an integer. It is initialized to the number of active coordinates.
  !   nb_cart:  is an integer. It is initialized to the number of cartesian coordinates.
  !   init_sub: is an integer. init_sub=1, when the initialization is done.
  ! remarks: by default, the input and output files are the standard ones (unit 5 or 6 in Fortran). 
  !          This can be changed with Init_InputUnit_Driver and Init_OutputUnit_Driver subroutines (see below).
```

```fortran
CALL Tnum_get_Qact0(Qact,nb_act)
  ! This subroutine gets the reference geometry from the initialization.
  !   Qact(:):  is a table of real of the active coordinates.
  !   nb_act:   is the number of active coordinates.
  ! Remark: useful when the reference geometry is given with Cartesian coordinates.
```

```fortran
CALL Tnum_Qact_TO_cart(Qact,nb_act,Qcart,nb_cart)
  ! This subroutine performs the initialization (if required) and generates the cartesian coordinates, Qcart, from the active ones, Qact.
  !   Qact(:):  is a table of real of the active coordinates.
  !   nb_act:   is the number of active coordinates. Its value needs to be the one obtained in the Init_TnumTana_Driver subroutine.
  !   Qcart(:): is a table of real of the active coordinates. It will contain the Cartesian coordinates
  !   nb_cart:  is the number of Cartesian coordinates. Its value needs to be the one obtained in the Init_TnumTana_Driver subroutine
  ! remarks: The Cartesian coordinates are not recentered with respect to the center of mass.
```

```fortran
CALL Tnum_Qact_TO_cartCOM(Qact,nb_act,Qcart,nb_cart)
  ! This subroutine performs the initialization (if required) and generates the active coordinates, Qact, from the Cartesian ones, Qcart.
  ! The variables are the same as the previous subroutine.
```

```fortran
CALL Tnum_cart_TO_Qact(Qact,nb_act,Qcart,nb_cart)
  ! It generates the active coordinates, Qact, from the Cartesian ones, Qcart.
  ! The variables are the same as the previous subroutines.
```

```fortran
CALL Tnum_get_d0GG(Qact,nb_act,d0GG,ndimG,def)
CALL Tnum_get_d0d1GG(Qact,nb_act,d0GG,d1GG,ndimG,def)
CALL Tnum_get_d0d1d2GG(Qact,nb_act,d0GG,d1GG,d2GG,ndimG,def)
  ! These subroutines compute the metric tensor (d0GG) and its derivatives (d1GG, d2GG) with respect to the actives coordinates.
  !   def: logical variable. When def=.TRUE. only the deformation part of the metric tensor is retrived (its size is (nb_act,nb_act)). When def=.FALSE., the deformation and the rotational parts are retrived (its size is (nb_act+3,nb_act+3))
  !   nb_act: as in the Init_TnumTana_Driver subroutine.
  !   ndimG:  is an integer. When def=.TRUE., its value must be nb_act. When def=.FALSE., its value must be nb_act+3, but it better to get it from the Tnum_get_ndimG(ndimG) subroutine.
  !   d0GG,d1GG,d2GG: they are, respectively, rank 2,3,4 tables. 
  !   The dimensions are d0GG(ndimG,ndimG), d1GG(ndimG,ndimG,nb_act), d2GG(ndimG,ndimG,nb_act,nb_act)
  ! Remark: see the main Fortran is TEST_Quantics.
```

### 3a) Link TnumTana library

In the TEST_Quantics directory, check the script **comp**.

1) Install Tnum-Tana package in directory (different the one of the main)
2) Compilation of the Library (see [section 1a](#installation-with-a-makefile)):

```bash
make OPT=0 lib
```
This will create several libraries (lib...a). The one required is: 
  **libTnum-TanaFull_gfortran_opt0_omp1_lapack1_int4_real64.a**

3) Compile of your main (main.f90) and linking with the Tnum-Tana library:

```bash
gfortran main.f90 -fopenmp libTnum-TanaFull_gfortran_opt0_omp1_lapack1_int4_real64.a -llapack -lblas
```

### 3b) Eckart transformation

The Fortran main is "APP/Main_Eckart.f90".
To compile the code :

```bash
make OPT=0 app
```
or
```bash
make OPT=0 app APPSRC=Main_Eckart.f90
```

The code reads two geometries: the reference one and then the current one. 
Remark: If masses are read (--read_masses t), the list of masses must be given after the reference geometry (the first one).

Then, the code computes the Eckart rotation matrix.
Depending of the argumant options:

- masses can be read (-m or --read_masses). They are read in g.mol$^{-1}$ or unified atomic mass unit (u or AMU) or Dalton (Da).
- the current transformed geometry is printed (-t or --coord_transfo)

Example, without reading the masses. Therefore, Tnum masses are used (CODATA2006 physical constants and NIST2012 masses).

```bash
./Main_Eckart.exe -i APP_InputFiles/dat_Eckart -o APP_OutputFiles/res --coord_transfo .TRUE.
```
or
```bash
./Main_Eckart.exe -i APP_InputFiles/dat_Eckart -o APP_OutputFiles/res --coord_transfo .TRUE. --read_masses f
```
The Eckart rotation matrix and the transformed Cartessian geometry are given in "APP_OutputFiles/res"

Example: reading the masses.

```bash
./Main_Eckart.exe -i APP_InputFiles/dat_Eckart2 -o APP_OutputFiles/res --coord_transfo .TRUE. --read_masses t
```

The Eckart rotation matrix and the transformed Cartessian geometry are given in "APP_OutputFiles/res"


### 3c) Inter-fragment coordinates

The Fortran main is "APP/Main_NMdimer.f90".
To compile the code :

```bash
make OPT=0 app
```
or
```bash
make OPT=0 app APPSRC=Main_NMdimer.f90
```

The code reads two fragment geometries (A and B) to form a dimer (AB). 
Remark: If masses are read (--read_masses t), the list of masses must be given after each geometry.

Then, the code computes the coordinates associated to the inter-fragment coordinates.
Depending of the argument options:

- masses can be read (-m or --read_masses). They are read in g.mol$^{-1}$ or unified atomic mass unit (u or AMU) or Dalton (Da).

Example, without reading the masses. Therefore, Tnum masses are used (CODATA2006 physical constants and NIST2012 masses).

```bash
./Main_NMdimer.exe -i APP_InputFiles/dat_dimer -o APP_OutputFiles/res
```
or
```bash
./Main_NMdimer.exe -i APP_InputFiles/dat_dimer -o APP_OutputFiles/res --read_masses f
```

The code also generates several coordinate motions (can be visual with jmol as vibrations):

- The translational and rotational motions of fragments, A and B.
- The inter-fragment motions (6 coordinates) of dimer AB.

Example: reading the masses.

```bash
./Main_NMdimer.exe -i APP_InputFiles/dat_dimer2 -o APP_OutputFiles/res --read_masses f
```

The inter-fragment coordinates of dimer AB are given at the end of "APP_OutputFiles/res" (search *Inter_dimAB*) as 6 linear combinations (in column) of Cartesian coordinates.