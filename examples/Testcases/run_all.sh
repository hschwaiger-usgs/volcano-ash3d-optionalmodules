#!/bin/bash
TestcasesDir=`pwd`
il=1 # 0-6: 0 is just no-limiter, 1 is NO and Superbee, 6 is all (limiters only apply to advection)
ix=3 # 1-5: 5 is the max prepared but that takes overnight

rc=0
which octave > /dev/null
rc=$((rc + $?))
if [[ "$rc" -gt 0 ]] ; then
  echo "ERROR: octave not found in path.  Post-processing requires octave."
  exit
fi

# Specify which cases to turn on/off, by setting the corresponding value to 1/0
    #  Case 1 : Horizontal advection
    #  Case 2 : Vertical advection
    #  Case 3 : Horizontal Rigid Rotation (Cone/block)
    #  Case 4 : Diffusion
    #  Case 5 : Horizontal Shear Rotation
    #  Case 6 : Method of Manufactured Solutions
    #  1 2 3 4 5 6
cases=(1 1 1 1 1 1)

if [ "$[cases[0]]" -ne 0 ]; then
  cd ${TestcasesDir}/TC1
  limmax=$il idx=$ix time ./TC1_XY_run
  ln -s TC1_XY_process.m TC1_XY_process
  ./TC1_XY_process
  limmax=$il idx=$ix time ./TC1_LL_run
  ln -s TC1_LL_process.m TC1_LL_process
  ./TC1_LL_process
  cd ${TestcasesDir}
fi

if [ "$[cases[1]]" -ne 0 ]; then
  cd ${TestcasesDir}/TC2
  limmax=$il idx=$ix time ./TC2_XY_run
  ln -s TC2_XY_process.m TC2_XY_process
  ./TC2_XY_process
  cd ${TestcasesDir}
fi

if [ "$[cases[2]]" -ne 0 ]; then
  cd ${TestcasesDir}/TC3
  limmax=$il idx=$ix time ./TC3_XY_run
  ln -s TC3_XY_process.m TC3_XY_process
  ./TC3_XY_process
  limmax=$il idx=$ix time ./TC3_LL_run
  ln -s TC3_LL_process.m TC3_LL_process
  ./TC3_LL_process
  cd ${TestcasesDir}
fi

if [ "$[cases[3]]" -ne 0 ]; then 
  cd ${TestcasesDir}/TC4
  limmax=$il idx=$ix time ./TC4_XY_run
  ln -s TC4_XY_process.m TC4_XY_process
  ./TC4_XY_process
  cd ${TestcasesDir}
fi

if [ "$[cases[4]]" -ne 0 ]; then
  cd ${TestcasesDir}/TC5
  limmax=$il idx=$ix time ./TC5_XY_run
  ln -s TC5_XY_process.m TC5_XY_process
  ./TC5_XY_process
  limmax=$il idx=$ix time ./TC5_LL_run
  ln -s TC5_LL_process.m TC5_LL_process
  ./TC5_LL_process
  cd ${TestcasesDir}
fi

if [ "$[cases[5]]" -ne 0 ]; then
  cd ${TestcasesDir}/TC6
  limmax=$il idx=$ix time ./TC6_XY_run
  ln -s TC6_XY_process.m TC6_XY_process
  ./TC6_XY_process
  cd ${TestcasesDir}
fi

# Reporting on how we did
if [ "$[cases[0]]" -ne 0 ]; then
  echo "------------------------------------------------------------------------------"
  echo "     Test Case 1: Horizontal advection"
  echo "------------------------------------------------------------------------------"
  echo "        XY"
  cat TC1/DATA/TC1_ConvRate_XY.dat
  echo "        LL"
  cat TC1/DATA/TC1_ConvRate_LL.dat
fi

if [ "$[cases[1]]" -ne 0 ]; then
  echo "------------------------------------------------------------------------------"
  echo "     Test Case 2: Vertical advection"
  echo "------------------------------------------------------------------------------"
  echo "        XY"
  cat TC2/DATA/TC2_ConvRate_XY.dat
fi

if [ "$[cases[2]]" -ne 0 ]; then
  echo "------------------------------------------------------------------------------"
  echo "     Test Case 3: Horizontal box/cone rigid rotation"
  echo "------------------------------------------------------------------------------"
  echo "        XY"
  cat TC3/DATA/TC3_ConvRate_XY.dat
  echo "        LL"
  cat TC3/DATA/TC3_ConvRate_LL.dat
fi

if [ "$[cases[3]]" -ne 0 ]; then
  echo "------------------------------------------------------------------------------"
  echo "     Test Case 4: Diffusion (Expl. and Crank-Nic.)"
  echo "------------------------------------------------------------------------------"
  echo "        XY"
  cat TC4/DATA/TC4_ConvRate_XY.dat
fi

if [ "$[cases[4]]" -ne 0 ]; then
  echo "------------------------------------------------------------------------------"
  echo "     Test Case 5: Horizontal shear rotation"
  echo "------------------------------------------------------------------------------"
  echo "        XY"
  cat TC5/DATA/TC5_ConvRate_XY.dat
  echo "        LL"
  cat TC5/DATA/TC5_ConvRate_LL.dat
fi

if [ "$[cases[5]]" -ne 0 ]; then
  echo "------------------------------------------------------------------------------"
  echo "     Test Case 6: Method of Manufactured Solutions"
  echo "------------------------------------------------------------------------------"
  echo "        XY"
  cat TC6/DATA/TC6_ConvRate_XY.dat
fi

