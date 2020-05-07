Matrices with tunable ∞-norm condition number and no need for pivoting in LU factorization
============================================================

This repository contains the code for producing the data in the tables in [sect. 5, 1].


Dependencies
------------

The code requires the GNU Scientific Library (GSL), the Intel Math Kernel Library (MKL), and OpenMPI.


Building the sources
--------------------

Issuing `make all` should compile the source code and produce the two executables `test_lowprec` and `test_norminf` used in the experiments. The Makefile assumes that the code is being compiled on a Linux machine and that Intel MKL 2019 is available in the default location `/opt/intel/compilers_and_libraries_2019`. The variables `CMAKE_PATH` and `CMAKE_INCLUDE` in the Makefile should be modified to run the experiment on a system with a different configuration.


Running the experiments
-----------------------

The binary file `test_lowprec` is a standalone program and can be run with a simple `./test_lowprec`. The binary file `test_norminf` should be used with the `mpirun` command, and the number of row and column processes to be used should be explicitly stated. The command

```
mpirun -n $(($nprows * $npcols)) ./test_norminf -m $nprows -n $npcols
```

should run to completion as long as the variables `nprows` and `npcols` are set and there are at least `$(($nprows * $npcols))` available slots in the system (if not, the command line option `--oversubscribe` can be used).


References
----------
If you use the code in this repository, please cite the preprint:

[1] M. Fasi and N. J. Higham. [*Matrices with tunable infinity-norm condition number and no need for pivoting in LU factorization*](http://eprints.maths.manchester.ac.uk/2775/).  MIMS EPrint 2020.17, Manchester Institute for Mathematical Sciences, The University of Manchester, UK, July 2020.


License
-------

This software is distributed under the terms of the GNU GPL v. 2 software license (see [LICENSE.md](./LICENSE.md)).
