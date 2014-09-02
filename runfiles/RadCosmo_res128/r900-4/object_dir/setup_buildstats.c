/* !!! DO NOT EDIT, FILES WRITTEN BY SETUP SCRIPT !!!
   
!!****f* object/setup_buildstats
!!
!! NAME
!!
!!  setup_buildstats
!!
!!
!! SYNOPSIS
!!
!!  call setup_buildstats(build_date, build_dir, build_machine, setup_call)
!!
!!  call setup_buildstats(character, character, character, character)
!!
!! 
!! DESCRIPTION
!!
!!  Simple subroutine generated at build time that returns the build date
!!  build directory, build machine, c and f flags and the full setup 
!!  command used to assemble the FLASH executable.
!!
!!
!!
!!***
*/

#include "Flash.h"
#include "constants.h"
#include "mangle_names.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

void FTOC(setup_buildstats)(char* build_date, 
		    char* build_dir, 
		    char* build_machine, 
		    char* setup_call, 
		    char* c_flags, 
		    char* f_flags){



     strncpy(build_date, "Tue Jun  3 13:31:10 CDT 2014",80);
     strncpy(build_dir, "/data1/r900-4/mepa/repos/flash4.2mepa/object", 80);
     strncpy(build_machine, "Linux r900-4.as.utexas.edu 2.6.32-279.2.1.el6.x86_64 #1 SMP Thu Jul 5 21:08:58 E", 80);
     strncpy(setup_call, "/data1/r900-4/mepa/repos/flash4.2mepa/bin/setup.py RadCosmo -auto -3d -maxblocks=600 ",400);
     strncpy(c_flags, "mpicc -I /home/r900-3/milos/hdf5-1.6.7/include -O2 -c -DMAXBLOCKS=600 -DNXB=8 -DNYB=8 -DNZB=8 -DN_DIM=3", 400);
     strncpy(f_flags, "mpif90 -f90=ifort -I/opt/r900-4/intel/composerxe-2011.1.107/composerxe-2011.1.107/mkl/include -c -i4 -r8 -O3 -static -traceback -DMAXBLOCKS=600 -DNXB=8 -DNYB=8 -DNZB=8 -DN_DIM=3", 400);


}

