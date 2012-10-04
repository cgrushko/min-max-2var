A 2-variable linear min-max solver, © (2012) Carmi Grushko.

License: GPL v3.

# Usage

  1. Edit Makefile to choose which external solvers you want. (e.g., EXTERNAL_SOLVER_LPSOLVE55). 
  2. Make sure the LIBRARY_PATH environment variable points to the different solvers (CGAL, GMP, MOSEK, etc.)
  3. Compile built-in external solvers
       
        cd lp_solve_5.5/lpsolve55
        sh ccc # or sh ccc.osx, and os on
        cd ../bfp/bfp_etaPFI/
        sh ccc # or sh ccc.osx, and os on
        cd ../../../seidel-lp-solver/
        gmake       
       
  2. Call `gmake' in the source root.
  3. Call `./a.out`. 
     a. If you want a specific solver, say `./a.out <number of solver>`.
     b. If you want to skip some problems, say `./a.out <number of solver> skip-list-file.txt`. See `seidel-skip-list.txt` for the format.
  4. Have a look at `settings.h`.

# Built-in solvers

The following external solvers are provided with the source code, and do not require any external dependencies.

  1. `lpsolve55`
  2. R. Seidel's LP solver
  
 
======

October 4, 2012