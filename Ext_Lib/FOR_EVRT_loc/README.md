# Some library for ElVibRot and Tnum/Tana

## 1) Installation with make

Make with optional options:

```bash
make FC=ifort OMP=0 OPT=0 LAPACK=0 INT=4
  # FC=ifort to change the compiller to ifort
  # OMP=0/1 to turn off/on the OpenMP fortran flag.
  # OPT=0/1 to turn off/on the fortran optimization.
  # LAPACK=0/1 to turn off/on the lapack use
  # INT=4/8 to change the default integer
```

It will run the test as well

To clean the library:

```bash
make cleanall
```


## 1) Installation with fpm

If it is not done yet, the external library directory links (in Ext_lib directory) need to be set by:

```bash
make fpm
```

Then 

```bash
fpm test
```

