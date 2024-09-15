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

In the code, we use the c-preprocessing to turn on/off the MPI.
But the actual fpm (0.9.0, alpha) does the c-preprocessing after checking the module uses. 
Therefore, when MPI is not installed on your computer, the fpm intallation fails when MPI is turn off by the c-preprocessing.

To overcome this problem run first:

```bash
make fpm
````

Then 

```bash
fpm test
```

