#!/bin/bash
#cd ../data/example1/
#gfortran -O3 ../../code/cg.f -o cg
#./cg inputfile1 > test1.log
#mkdir -p ../../results/example1/
#mv test1* ../../results/example1/

# example2 requires a few days to run
cd ../data/example2/
gfortran -O3 ../../code/cg.f -o cg
./cg inputfile2 > test2.log
mkdir -p ../../results/example2/
mv test2* ../../results/example2/

# example3 requires more than a month to run
# if openmp is available, using 4 cores is advised
#cd ../../data/example3/
#ulimit -s unlimited
#gfortran -O3 -fopenmp cg.f -o pid
#env OMP_NUM_THREADS=4 ./cg inputfile3 > test3.log
#mkdir -p ../../results/example3/
#mv test3* ../../results/example3/
