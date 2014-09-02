!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! File created at setup time.  DO NOT EDIT !!!

subroutine Simulation_mapIntToStr(key, str, map)
    use Grid_interface, only: Grid_formatNonRep
    use Driver_interface, only: Driver_getMype, Driver_getNumProcs
    use RuntimeParameters_interface, only: RuntimeParameters_get
    implicit none 
#include "constants.h"
#include "Flash.h"
    integer, intent(in) :: key, map
    character(len=*), intent(inout) :: str

    integer, parameter :: locunk1(0:NONREP_COUNT) = NONREP_LOCUNK1
    integer, parameter :: maxlocs(0:NONREP_COUNT) = NONREP_MAXLOCS
    character(len=*), parameter :: rpcount_flat = NONREP_RPCOUNT_FLAT
    integer, parameter :: rpcount_start(NONREP_COUNT+1) = NONREP_RPCOUNT_START
    
    integer :: mesh, meshes
    integer :: k, nonrep, nglob, iglob, iloc
    
    call Driver_getMype(MESH_ACROSS_COMM, mesh)
    call Driver_getNumProcs(MESH_ACROSS_COMM, meshes)
 
    k = key + map*MAPBLOCKSIZE
    select case(k)
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  1); str="h2io"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  2); str="heat"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  3); str="heio"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  4); str="hepi"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  5); str="hion"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  6); str="aux"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  7); str="c011"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  8); str="c012"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  9); str="c013"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 10); str="c014"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 11); str="c015"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 12); str="c016"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 13); str="c041"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 14); str="c042"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 15); str="c043"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 16); str="c044"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 17); str="c045"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 18); str="c046"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 19); str="c071"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 20); str="c072"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 21); str="c073"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 22); str="c074"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 23); str="c075"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 24); str="c076"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 25); str="c081"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 26); str="c082"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 27); str="c083"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 28); str="c084"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 29); str="c085"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 30); str="c086"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 31); str="chem"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 32); str="dens"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 33); str="eint"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 34); str="ener"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 35); str="flh2"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 36); str="flx1"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 37); str="flx2"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 38); str="flx3"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 39); str="gamc"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 40); str="game"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 41); str="gpol"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 42); str="gpot"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 43); str="grac"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 44); str="icor"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 45); str="imgm"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 46); str="imgp"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 47); str="isls"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 48); str="pde"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 49); str="pden"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 50); str="pres"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 51); str="sim1"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 52); str="temp"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 53); str="velx"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 54); str="vely"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 55); str="velz"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 56); str="deut"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 57); str="dplu"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 58); str="elec"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 59); str="h"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 60); str="hd"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 61); str="hel"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 62); str="hep"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 63); str="hepp"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 64); str="hmin"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 65); str="hplu"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 66); str="htwo"
    case((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 67); str="htwp"
    case((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  1); str="e"
    case((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  2); str="eint"
    case((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  3); str="p"
    case((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  4); str="rho"
    case((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  5); str="u"
    case((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  6); str="ut"
    case((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  7); str="utt"
    case((MAPBLOCK_PART * MAPBLOCKSIZE)+  1); str="accx"
    case((MAPBLOCK_PART * MAPBLOCKSIZE)+  2); str="accy"
    case((MAPBLOCK_PART * MAPBLOCKSIZE)+  3); str="accz"
    case((MAPBLOCK_PART * MAPBLOCKSIZE)+  4); str="blk"
    case((MAPBLOCK_PART * MAPBLOCKSIZE)+  5); str="gpot"
    case((MAPBLOCK_PART * MAPBLOCKSIZE)+  6); str="grid_dens"
    case((MAPBLOCK_PART * MAPBLOCKSIZE)+  7); str="mass"
    case((MAPBLOCK_PART * MAPBLOCKSIZE)+  8); str="mass_save"
    case((MAPBLOCK_PART * MAPBLOCKSIZE)+  9); str="oacx"
    case((MAPBLOCK_PART * MAPBLOCKSIZE)+ 10); str="oacy"
    case((MAPBLOCK_PART * MAPBLOCKSIZE)+ 11); str="oacz"
    case((MAPBLOCK_PART * MAPBLOCKSIZE)+ 12); str="posx"
    case((MAPBLOCK_PART * MAPBLOCKSIZE)+ 13); str="posy"
    case((MAPBLOCK_PART * MAPBLOCKSIZE)+ 14); str="posz"
    case((MAPBLOCK_PART * MAPBLOCKSIZE)+ 15); str="proc"
    case((MAPBLOCK_PART * MAPBLOCKSIZE)+ 16); str="smoothtag"
    case((MAPBLOCK_PART * MAPBLOCKSIZE)+ 17); str="tag"
    case((MAPBLOCK_PART * MAPBLOCKSIZE)+ 18); str="velx"
    case((MAPBLOCK_PART * MAPBLOCKSIZE)+ 19); str="vely"
    case((MAPBLOCK_PART * MAPBLOCKSIZE)+ 20); str="velz"
    case default; str = "err"
    end select

    do nonrep=1, NONREP_COUNT
        if(locunk1(nonrep) <= k .and. k - locunk1(nonrep) < maxlocs(nonrep)) then
            iloc = k - locunk1(nonrep) + 1
            ! get the size of this nonrep array
            call RuntimeParameters_get(rpcount_flat(rpcount_start(nonrep):rpcount_start(nonrep+1)-1), nglob)
            iglob = NONREP_LOC2GLOB(iloc, mesh, meshes)
            if(iglob .gt. nglob) then
                str = "err"
                return
            end if
            call Grid_formatNonRep(nonrep, iglob, str)
            exit
        end if
    end do
end subroutine Simulation_mapIntToStr
