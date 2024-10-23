Ash3d-OptionalModules
==========

volcano-ash3d-optionalmodules is a repository the illustrates how to build user-created
optional modules to link into the core Ash3d source code. In this case, the module
built is `TestCases`.

To build the module, change to the `src` directory and edit `makefile_optmod`.
You will need to ensure that the path to the Ash3d repository is set
correctly for your system.

  `ASH3dCCSRC=~/work/USGS/Software/GIT/volcano-ash3d/src`

Once this path is set, you can test the build process by typing:

  `make -f makefile_optmod`

This will build the module in:

`volcano-ash3d-optionalmodules/src/Optional_Modules/TestCases/Testcases.f90`

Note that there two other files that replace the counter parts in `volcano-ash3d`:

1. `Ash3d_TC.F90` 
This file is a local copy of the the top-level file of the Ash3d code, but with the
calls to the optional module subroutines as needed. The design goal is that the
optional modules should be minimally invasive and hopefully only require editing
this top-level file.

2. `Set_BC_TC.f90`
Unfortunately, sometimes a module requires editing additional file of the core-code.
In this case, the method of manufactured solutions test case requires a more general
treatment of boundary conditions.

As with the makefile from the core-code repository, you can edit the settings for
`SYSTEM` to change the compiler, `RUN` to switch between debugging, profiling, or
optimized code, among various other settings.

Usage
-----

Other than being an example of how to build custom modules against the core-code,
this repository tests the performance of the Ash3d routines using as series of
convergence tests:

1. Linear horizontal advection (both Cartesian and Spherical)
    1. x+y0
    2. x-y0
    3. x0y+
    4. x0y-
    5. x+y+
    6. x-y-
    7. x-y+
    8. x+y-

2. Linear vertical advection
    1. z+ vf0
    2. z- vf0
    3. z0 vf+
    4. z0 vf-

3. Horizontal rotation of block and cone (both Cartesian and Spherical)

4. Diffusion (Cartesian only)
    1. Explicit Diffusion in x
    2. Explicit Diffusion in y
    3. Explicit Diffusion in z
    4. Crank-Nicolson in x
    5. Crank-Nicolson in y
    6. Crank-Nicolson in z

5. Horizontal rotational shear (both Cartesian and Spherical)

6. Method of manufactured solutions (Cartesian only)

To run all the test cases, change to the directory:  
`volcano-ash3d-optionalmodules/examples/Testcases`
and edit the `run_all.sh` script.  

To control how many limiters to use, set `il` to a number from 0-6, where
0=No limiter, 1=Superbee, 2=Lax-Wendroff, 3=Beam-Warming, 4=Fromm, 5=Minmod,
and 6=MC. The default is `il=1` which means both the No-limiter and Superbee
cases are run.  

To control the number of refinement steps, set `ix` to be 1-5.  

To control which test cases are run, toggle the corresponding elements of `cases=()`
from 0 (do not run) to 1 (run case).  

Note that post-processing requires `octave`, both for plotting and to calculate convergence rate.

Authors
-------

Hans F. Schwaiger <hschwaiger@usgs.gov>  

