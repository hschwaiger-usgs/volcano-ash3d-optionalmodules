#!/bin/bash
TestcasesDir=`pwd`
il=1 # 0 is just no-limiter, 1 in NO and Superbee, 6 is all
ix=3 # 5 is the max coded but that takes overnight

cd ${TestcasesDir}/TC1
limmax=$il idx=$ix time ./TC1_XY_run
ln -s TC1_XY_process.m TC1_XY_process
./TC1_XY_process
limmax=$il idx=$ix time ./TC1_LL_run
ln -s TC1_LL_process.m TC1_LL_process
./TC1_LL_process

cd ${TestcasesDir}/TC2
limmax=$il idx=$ix time ./TC2_XY_run
ln -s TC2_XY_process.m TC2_XY_process
./TC2_XY_process

cd ${TestcasesDir}/TC3
limmax=$il idx=$ix time ./TC3_XY_run
ln -s TC3_XY_process.m TC3_XY_process
./TC3_XY_process
limmax=$il idx=$ix time ./TC3_LL_run
ln -s TC3_LL_process.m TC3_LL_process
./TC3_LL_process

cd ${TestcasesDir}/TC4
limmax=$il idx=$ix time ./TC4_XY_run
ln -s TC4_XY_process.m TC4_XY_process
./TC4_XY_process

cd ${TestcasesDir}/TC5
limmax=$il idx=$ix time ./TC5_XY_run
ln -s TC5_XY_process.m TC5_XY_process
./TC5_XY_process
limmax=$il idx=$ix time ./TC5_LL_run
ln -s TC5_LL_process.m TC5_LL_process
./TC5_LL_process

echo "------------------------------------------------------------------------------"
echo "     Test Case 1: Horizontal advection"
echo "------------------------------------------------------------------------------"
echo "        XY"
cat TC1/DATA/TC1_ConvRate_XY.dat
echo "        LL"
cat TC1/DATA/TC1_ConvRate_LL.dat
echo "------------------------------------------------------------------------------"
echo "     Test Case 2: Vertical advection"
echo "------------------------------------------------------------------------------"
echo "        XY"
cat TC2/DATA/TC2_ConvRate_XY.dat
echo "------------------------------------------------------------------------------"
echo "     Test Case 3: Horizontal box/cone rigid rotation"
echo "------------------------------------------------------------------------------"
echo "        XY"
cat TC3/DATA/TC3_ConvRate_XY.dat
echo "        LL"
cat TC3/DATA/TC3_ConvRate_LL.dat
echo "------------------------------------------------------------------------------"
echo "     Test Case 4: Diffusion (Expl. and Crank-Nic."
echo "------------------------------------------------------------------------------"
echo "        XY"
cat TC4/DATA/TC4_ConvRate_XY.dat
echo "------------------------------------------------------------------------------"
echo "     Test Case 5: Horizontal shear rotation"
echo "------------------------------------------------------------------------------"
echo "        XY"
cat TC5/DATA/TC5_ConvRate_XY.dat
echo "        LL"
cat TC5/DATA/TC5_ConvRate_LL.dat







