#!/bin/bash

# Select Fortran compiler based on hostname.
llapack='-llapack'
ilapack=''
case "$(hostname 2>/dev/null)" in
pc*.tcm.phy.cam.ac.uk)
  F=mpgfortran ;;
allogin*)
  . /usr/share/modules/init/bash
  module load ifort mpi.intel
  F=mpiifort ;;
vortex*)
  . /etc/profile.d/modules.sh
  module load intel impi mkl
  F=mpiifort
  llapack="${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl"
  ilapack="-I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include" ;;
*)
  F=mpif90 ;;
esac

# Choose compiler options.
case "$F" in
mpif90*|*gfortran*)
  if [ "$1" = -d ] ; then
    opts="-Wall -Wextra -fimplicit-none -O0 -fbounds-check -g -pg -pedantic\
       -fbacktrace -ffpe-trap=invalid,zero,overflow"
  else
    opts="-Ofast -fprotect-parens -march=native"
  fi
  case "$(hostname 2>/dev/null)" in
  pcal*) : ;;
  *) opts="$opts -fdiagnostics-color=always" ;;
  esac ;;
mpifort*|mpiifort*)
  if [ "$1" = -d ] ; then
    opts="-check all -check noarg_temp_created -fp-stack-check -g -O0\
       -implicitnone -std95 -traceback -warn all,nounused -debug all -ftrapuv\
       -assume buffered_io"
  else
    opts="-O3 -no-prec-div -no-prec-sqrt -funroll-loops -no-fp-port -ip\
       -complex-limited-range -assume protect_parens -assume buffered_io"
  fi ;;
mpnagfor*)
  if [ "$1" = -d ] ; then
    opts="-g -C=all -colour -f95 -strict95 -u -v -gline -nan\
       -no_underflow_warning -w=uda -O0"
  else
    opts="-O4 -Ounsafe"
  fi ;;
*)
  echo "Compiler '$F' not configured."
  exit 1 ;;
esac

$F $opts $ilapack -o treat treat.f90 $llapack 2>&1 | grep -vE "Unused parameter.*mpi_|mpif|^$"
