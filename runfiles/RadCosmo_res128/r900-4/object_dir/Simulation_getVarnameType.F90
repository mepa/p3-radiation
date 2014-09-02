!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! File created at setup time.  DO NOT EDIT !!!
!!

subroutine Simulation_getVarnameType(varname,vartype)
implicit none 

#include "constants.h"
#include "Flash.h"

   integer, intent(out) :: vartype
   integer, intent(in) :: varname
   integer :: ctr

   !! varnames is an array of possible input strings (with possibly extra spaces)
   !! vartypes is an array of corresponding answers

   integer, dimension(1:55), parameter :: varnames = (/&
   & H2IO_VAR, HEAT_VAR, HEIO_VAR, HEPI_VAR, HION_VAR, &
     AUX_VAR, C011_VAR, C012_VAR, C013_VAR, C014_VAR, &
     C015_VAR, C016_VAR, C041_VAR, C042_VAR, C043_VAR, &
     C044_VAR, C045_VAR, C046_VAR, C071_VAR, C072_VAR, &
     C073_VAR, C074_VAR, C075_VAR, C076_VAR, C081_VAR, &
     C082_VAR, C083_VAR, C084_VAR, C085_VAR, C086_VAR, &
     CHEM_VAR, DENS_VAR, EINT_VAR, ENER_VAR, FLH2_VAR, &
     FLX1_VAR, FLX2_VAR, FLX3_VAR, GAMC_VAR, GAME_VAR, &
     GPOL_VAR, GPOT_VAR, GRAC_VAR, ICOR_VAR, IMGM_VAR, &
     IMGP_VAR, ISLS_VAR, PDE_VAR, PDEN_VAR, PRES_VAR, &
     SIM1_VAR, TEMP_VAR, VELX_VAR, VELY_VAR, VELZ_VAR &
   &/)

   integer, dimension(1:55), parameter :: vartypes = (/&
   & VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC,&
     VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC,&
     VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC,&
     VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC,&
     VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC,&
     VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC,&
     VARTYPE_GENERIC, VARTYPE_PER_VOLUME, VARTYPE_PER_MASS, VARTYPE_PER_MASS, VARTYPE_GENERIC,&
     VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC,&
     VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC,&
     VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_PER_VOLUME, VARTYPE_GENERIC,&
     VARTYPE_GENERIC, VARTYPE_GENERIC, VARTYPE_PER_MASS, VARTYPE_PER_MASS, VARTYPE_PER_MASS &
   &/)

   ! Species and mass scalars are kept in primitive form. - KW
   if (varname .ge. SPECIES_BEGIN .and. varname .le. MASS_SCALARS_END) then
      vartype = VARTYPE_PER_MASS
      return
   end if

   vartype = VARTYPE_ERROR  ! By default vartype is an ERROR
   do ctr = LBOUND(varnames,1), UBOUND(varnames,1)
      if ( varnames(ctr) .eq. varname ) vartype = vartypes(ctr)
   end do

end subroutine Simulation_getVarnameType

