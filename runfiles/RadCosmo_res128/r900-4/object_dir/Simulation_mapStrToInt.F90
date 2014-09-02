!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! File created at setup time.  DO NOT EDIT !!!
!!

subroutine Simulation_mapStrToInt(str,key,map)
    use Grid_interface, only: Grid_parseNonRep
    use Driver_interface, only: Driver_getMype, Driver_getNumProcs
    use RuntimeParameters_interface, only: RuntimeParameters_get
    implicit none
#include "constants.h"
#include "Flash.h"
    character(len=*), intent(in) :: str
    integer, intent(out) :: key 
    integer, intent(in) :: map

    integer, parameter :: locunk1(0:NONREP_COUNT) = NONREP_LOCUNK1
    character(*), parameter :: rpcount_flat = NONREP_RPCOUNT_FLAT
    integer, parameter :: rpcount_start(NONREP_COUNT+1) = NONREP_RPCOUNT_START
    
    integer :: mesh, meshes
    character(len=MAX_STRING_LENGTH) :: strlwr
    integer :: nonrep, glob, nglob
    
    call Driver_getMype(MESH_ACROSS_COMM, mesh)
    call Driver_getNumProcs(MESH_ACROSS_COMM, meshes)
    key = NONEXISTENT
    strlwr = str
    call makeLowercase(strlwr)
    
    call Grid_parseNonRep(strlwr(1:len(str)), nonrep, glob)
    if(nonrep .gt. 0) then
        call RuntimeParameters_get(rpcount_flat(rpcount_start(nonrep):rpcount_start(nonrep+1)-1), nglob)
        if(glob .gt. nglob .or. mesh .ne. NONREP_MESHOFGLOB(glob, meshes)) return ! NONEXISTENT
        key = locunk1(nonrep)-1 + NONREP_GLOB2LOC(glob, mesh, meshes)
        return
    end if
    
    select case(strlwr)
    case("ener")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 34)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 34
        end select
    case("htwp")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 67)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 67
        end select
    case("htwo")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 66)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 66
        end select
    case("c085")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 29)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 29
        end select
    case("c084")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 28)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 28
        end select
    case("c086")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 30)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 30
        end select
    case("c081")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 25)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 25
        end select
    case("c083")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 27)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 27
        end select
    case("c082")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 26)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 26
        end select
    case("gpol")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 41)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 41
        end select
    case("aux")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  6)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  6
        end select
    case("c071")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 19)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 19
        end select
    case("c072")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 20)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 20
        end select
    case("c073")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 21)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 21
        end select
    case("c074")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 22)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 22
        end select
    case("c075")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 23)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 23
        end select
    case("c076")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 24)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 24
        end select
    case("game")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 40)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 40
        end select
    case("gamc")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 39)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 39
        end select
    case("h")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 59)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 59
        end select
    case("hplu")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 65)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 65
        end select
    case("p")
        select case(map)
        case(((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  3)/MAPBLOCKSIZE); key = (MAPBLOCK_FLUX * MAPBLOCKSIZE)+  3
        end select
    case("grac")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 43)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 43
        end select
    case("hion")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  5)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  5
        end select
    case("flx2")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 37)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 37
        end select
    case("flx3")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 38)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 38
        end select
    case("flx1")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 36)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 36
        end select
    case("velz")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 55)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 55
        case(((MAPBLOCK_PART * MAPBLOCKSIZE)+ 20)/MAPBLOCKSIZE); key = (MAPBLOCK_PART * MAPBLOCKSIZE)+ 20
        end select
    case("vely")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 54)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 54
        case(((MAPBLOCK_PART * MAPBLOCKSIZE)+ 19)/MAPBLOCKSIZE); key = (MAPBLOCK_PART * MAPBLOCKSIZE)+ 19
        end select
    case("velx")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 53)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 53
        case(((MAPBLOCK_PART * MAPBLOCKSIZE)+ 18)/MAPBLOCKSIZE); key = (MAPBLOCK_PART * MAPBLOCKSIZE)+ 18
        end select
    case("chem")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 31)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 31
        end select
    case("grid_dens")
        select case(map)
        case(((MAPBLOCK_PART * MAPBLOCKSIZE)+  6)/MAPBLOCKSIZE); key = (MAPBLOCK_PART * MAPBLOCKSIZE)+  6
        end select
    case("blk")
        select case(map)
        case(((MAPBLOCK_PART * MAPBLOCKSIZE)+  4)/MAPBLOCKSIZE); key = (MAPBLOCK_PART * MAPBLOCKSIZE)+  4
        end select
    case("proc")
        select case(map)
        case(((MAPBLOCK_PART * MAPBLOCKSIZE)+ 15)/MAPBLOCKSIZE); key = (MAPBLOCK_PART * MAPBLOCKSIZE)+ 15
        end select
    case("oacy")
        select case(map)
        case(((MAPBLOCK_PART * MAPBLOCKSIZE)+ 10)/MAPBLOCKSIZE); key = (MAPBLOCK_PART * MAPBLOCKSIZE)+ 10
        end select
    case("oacx")
        select case(map)
        case(((MAPBLOCK_PART * MAPBLOCKSIZE)+  9)/MAPBLOCKSIZE); key = (MAPBLOCK_PART * MAPBLOCKSIZE)+  9
        end select
    case("oacz")
        select case(map)
        case(((MAPBLOCK_PART * MAPBLOCKSIZE)+ 11)/MAPBLOCKSIZE); key = (MAPBLOCK_PART * MAPBLOCKSIZE)+ 11
        end select
    case("posz")
        select case(map)
        case(((MAPBLOCK_PART * MAPBLOCKSIZE)+ 14)/MAPBLOCKSIZE); key = (MAPBLOCK_PART * MAPBLOCKSIZE)+ 14
        end select
    case("posx")
        select case(map)
        case(((MAPBLOCK_PART * MAPBLOCKSIZE)+ 12)/MAPBLOCKSIZE); key = (MAPBLOCK_PART * MAPBLOCKSIZE)+ 12
        end select
    case("posy")
        select case(map)
        case(((MAPBLOCK_PART * MAPBLOCKSIZE)+ 13)/MAPBLOCKSIZE); key = (MAPBLOCK_PART * MAPBLOCKSIZE)+ 13
        end select
    case("flh2")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 35)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 35
        end select
    case("pde")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 48)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 48
        end select
    case("eint")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 33)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 33
        case(((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  2)/MAPBLOCKSIZE); key = (MAPBLOCK_FLUX * MAPBLOCKSIZE)+  2
        end select
    case("c041")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 13)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 13
        end select
    case("c043")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 15)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 15
        end select
    case("c042")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 14)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 14
        end select
    case("c045")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 17)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 17
        end select
    case("c044")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 16)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 16
        end select
    case("c046")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 18)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 18
        end select
    case("gpot")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 42)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 42
        case(((MAPBLOCK_PART * MAPBLOCKSIZE)+  5)/MAPBLOCKSIZE); key = (MAPBLOCK_PART * MAPBLOCKSIZE)+  5
        end select
    case("accy")
        select case(map)
        case(((MAPBLOCK_PART * MAPBLOCKSIZE)+  2)/MAPBLOCKSIZE); key = (MAPBLOCK_PART * MAPBLOCKSIZE)+  2
        end select
    case("accx")
        select case(map)
        case(((MAPBLOCK_PART * MAPBLOCKSIZE)+  1)/MAPBLOCKSIZE); key = (MAPBLOCK_PART * MAPBLOCKSIZE)+  1
        end select
    case("accz")
        select case(map)
        case(((MAPBLOCK_PART * MAPBLOCKSIZE)+  3)/MAPBLOCKSIZE); key = (MAPBLOCK_PART * MAPBLOCKSIZE)+  3
        end select
    case("dplu")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 57)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 57
        end select
    case("c011")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  7)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  7
        end select
    case("hep")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 62)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 62
        end select
    case("hel")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 61)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 61
        end select
    case("pres")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 50)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 50
        end select
    case("pden")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 49)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 49
        end select
    case("heio")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  3)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  3
        end select
    case("heat")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  2)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  2
        end select
    case("mass_save")
        select case(map)
        case(((MAPBLOCK_PART * MAPBLOCKSIZE)+  8)/MAPBLOCKSIZE); key = (MAPBLOCK_PART * MAPBLOCKSIZE)+  8
        end select
    case("hepi")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  4)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  4
        end select
    case("hd")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 60)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 60
        end select
    case("temp")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 52)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 52
        end select
    case("ut")
        select case(map)
        case(((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  6)/MAPBLOCKSIZE); key = (MAPBLOCK_FLUX * MAPBLOCKSIZE)+  6
        end select
    case("hmin")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 64)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 64
        end select
    case("tag")
        select case(map)
        case(((MAPBLOCK_PART * MAPBLOCKSIZE)+ 17)/MAPBLOCKSIZE); key = (MAPBLOCK_PART * MAPBLOCKSIZE)+ 17
        end select
    case("utt")
        select case(map)
        case(((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  7)/MAPBLOCKSIZE); key = (MAPBLOCK_FLUX * MAPBLOCKSIZE)+  7
        end select
    case("hepp")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 63)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 63
        end select
    case("sim1")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 51)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 51
        end select
    case("smoothtag")
        select case(map)
        case(((MAPBLOCK_PART * MAPBLOCKSIZE)+ 16)/MAPBLOCKSIZE); key = (MAPBLOCK_PART * MAPBLOCKSIZE)+ 16
        end select
    case("isls")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 47)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 47
        end select
    case("mass")
        select case(map)
        case(((MAPBLOCK_PART * MAPBLOCKSIZE)+  7)/MAPBLOCKSIZE); key = (MAPBLOCK_PART * MAPBLOCKSIZE)+  7
        end select
    case("deut")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 56)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 56
        end select
    case("icor")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 44)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 44
        end select
    case("imgm")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 45)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 45
        end select
    case("h2io")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  1)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  1
        end select
    case("rho")
        select case(map)
        case(((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  4)/MAPBLOCKSIZE); key = (MAPBLOCK_FLUX * MAPBLOCKSIZE)+  4
        end select
    case("c012")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  8)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  8
        end select
    case("c013")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) +  9)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) +  9
        end select
    case("imgp")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 46)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 46
        end select
    case("c016")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 12)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 12
        end select
    case("c014")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 10)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 10
        end select
    case("c015")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 11)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 11
        end select
    case("elec")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 58)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 58
        end select
    case("e")
        select case(map)
        case(((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  1)/MAPBLOCKSIZE); key = (MAPBLOCK_FLUX * MAPBLOCKSIZE)+  1
        end select
    case("dens")
        select case(map)
        case(((MAPBLOCK_UNK  * MAPBLOCKSIZE) + 32)/MAPBLOCKSIZE); key = (MAPBLOCK_UNK  * MAPBLOCKSIZE) + 32
        end select
    case("u")
        select case(map)
        case(((MAPBLOCK_FLUX * MAPBLOCKSIZE)+  5)/MAPBLOCKSIZE); key = (MAPBLOCK_FLUX * MAPBLOCKSIZE)+  5
        end select
    end select

    if(key .ne. NONEXISTENT) key = mod(key,MAPBLOCKSIZE)
end subroutine Simulation_mapStrToInt
