!!****f* object/setup_flashUnits
!!
!! NAME
!!
!!  setup_getFlashUnits
!!
!!
!! SYNOPSIS
!!
!!
!!  setup_getFlashUnits(unit_names)
!!
!!  setup_getFlashUnits(character())
!!
!!
!! DESCRIPTION
!!
!!  Return a character array of size NUM_UNITS containing
!!  the names of all of the FLASH units used to assemble
!!  the current executable
!!
!!  The unit_names variable should be declared as
!!
!!    use flashUnits
!!
!!  
!!    character (len=MAX_STRING_LENGTH) :: flash_units(NUM_UNITS) 
!!
!!
!!  The length of each character string is set to MAX_STRING_LENGTH,
!!  which is defined in the automatically generated flash_defines.fh
!!
!!***

  subroutine setup_getFlashUnits(unit_names)

#include "constants.h"
    implicit none

    integer, PARAMETER :: NUM_UNITS = 137
    character (len=MAX_STRING_LENGTH) :: unit_names(NUM_UNITS)
    integer :: i

    i = 0

    i = i + 1; unit_names(i) = &
"Driver"
    i = i + 1; unit_names(i) = &
"Driver/DriverMain"
    i = i + 1; unit_names(i) = &
"Driver/DriverMain/Split"
    i = i + 1; unit_names(i) = &
"Grid"
    i = i + 1; unit_names(i) = &
"Grid/GridBoundaryConditions"
    i = i + 1; unit_names(i) = &
"Grid/GridMain"
    i = i + 1; unit_names(i) = &
"Grid/GridMain/paramesh"
    i = i + 1; unit_names(i) = &
"Grid/GridMain/paramesh/interpolation"
    i = i + 1; unit_names(i) = &
"Grid/GridMain/paramesh/interpolation/Paramesh4"
    i = i + 1; unit_names(i) = &
"Grid/GridMain/paramesh/interpolation/Paramesh4/prolong"
    i = i + 1; unit_names(i) = &
"Grid/GridMain/paramesh/interpolation/prolong"
    i = i + 1; unit_names(i) = &
"Grid/GridMain/paramesh/paramesh4"
    i = i + 1; unit_names(i) = &
"Grid/GridMain/paramesh/paramesh4/Paramesh4dev"
    i = i + 1; unit_names(i) = &
"Grid/GridMain/paramesh/paramesh4/Paramesh4dev/PM4_package"
    i = i + 1; unit_names(i) = &
"Grid/GridMain/paramesh/paramesh4/Paramesh4dev/PM4_package/headers"
    i = i + 1; unit_names(i) = &
"Grid/GridMain/paramesh/paramesh4/Paramesh4dev/PM4_package/mpi_source"
    i = i + 1; unit_names(i) = &
"Grid/GridMain/paramesh/paramesh4/Paramesh4dev/PM4_package/source"
    i = i + 1; unit_names(i) = &
"Grid/GridMain/paramesh/paramesh4/Paramesh4dev/PM4_package/utilities"
    i = i + 1; unit_names(i) = &
"Grid/GridMain/paramesh/paramesh4/Paramesh4dev/PM4_package/utilities/multigrid"
    i = i + 1; unit_names(i) = &
"Grid/GridParticles"
    i = i + 1; unit_names(i) = &
"Grid/GridParticles/GridParticlesMapFromMesh"
    i = i + 1; unit_names(i) = &
"Grid/GridParticles/GridParticlesMapToMesh"
    i = i + 1; unit_names(i) = &
"Grid/GridParticles/GridParticlesMapToMesh/Paramesh"
    i = i + 1; unit_names(i) = &
"Grid/GridParticles/GridParticlesMapToMesh/Paramesh/MoveSieve"
    i = i + 1; unit_names(i) = &
"Grid/GridParticles/GridParticlesMove"
    i = i + 1; unit_names(i) = &
"Grid/GridParticles/GridParticlesMove/Sieve"
    i = i + 1; unit_names(i) = &
"Grid/GridParticles/GridParticlesMove/Sieve/BlockMatch"
    i = i + 1; unit_names(i) = &
"Grid/GridParticles/GridParticlesMove/paramesh"
    i = i + 1; unit_names(i) = &
"Grid/GridSolvers"
    i = i + 1; unit_names(i) = &
"Grid/GridSolvers/IsoBndMultipole"
    i = i + 1; unit_names(i) = &
"Grid/GridSolvers/Multigrid"
    i = i + 1; unit_names(i) = &
"Grid/GridSolvers/Multigrid/fft"
    i = i + 1; unit_names(i) = &
"Grid/localAPI"
    i = i + 1; unit_names(i) = &
"IO"
    i = i + 1; unit_names(i) = &
"IO/IOMain"
    i = i + 1; unit_names(i) = &
"IO/IOMain/hdf5"
    i = i + 1; unit_names(i) = &
"IO/IOMain/hdf5/serial"
    i = i + 1; unit_names(i) = &
"IO/IOMain/hdf5/serial/PM"
    i = i + 1; unit_names(i) = &
"IO/IOParticles"
    i = i + 1; unit_names(i) = &
"IO/IOParticles/hdf5"
    i = i + 1; unit_names(i) = &
"IO/IOParticles/hdf5/serial"
    i = i + 1; unit_names(i) = &
"IO/localAPI"
    i = i + 1; unit_names(i) = &
"Multispecies"
    i = i + 1; unit_names(i) = &
"Multispecies/MultispeciesMain"
    i = i + 1; unit_names(i) = &
"Particles"
    i = i + 1; unit_names(i) = &
"Particles/ParticlesForces"
    i = i + 1; unit_names(i) = &
"Particles/ParticlesForces/longRange"
    i = i + 1; unit_names(i) = &
"Particles/ParticlesForces/longRange/gravity"
    i = i + 1; unit_names(i) = &
"Particles/ParticlesForces/longRange/gravity/ParticleMesh"
    i = i + 1; unit_names(i) = &
"Particles/ParticlesInitialization"
    i = i + 1; unit_names(i) = &
"Particles/ParticlesMain"
    i = i + 1; unit_names(i) = &
"Particles/ParticlesMain/active"
    i = i + 1; unit_names(i) = &
"Particles/ParticlesMain/active/massive"
    i = i + 1; unit_names(i) = &
"Particles/ParticlesMain/active/massive/LeapfrogCosmo"
    i = i + 1; unit_names(i) = &
"Particles/ParticlesMapping"
    i = i + 1; unit_names(i) = &
"Particles/ParticlesMapping/meshWeighting"
    i = i + 1; unit_names(i) = &
"Particles/ParticlesMapping/meshWeighting/CIC"
    i = i + 1; unit_names(i) = &
"Particles/ParticlesMapping/meshWeighting/MapToMesh"
    i = i + 1; unit_names(i) = &
"Particles/localAPI"
    i = i + 1; unit_names(i) = &
"PhysicalConstants"
    i = i + 1; unit_names(i) = &
"PhysicalConstants/PhysicalConstantsMain"
    i = i + 1; unit_names(i) = &
"RuntimeParameters"
    i = i + 1; unit_names(i) = &
"RuntimeParameters/RuntimeParametersMain"
    i = i + 1; unit_names(i) = &
"Simulation"
    i = i + 1; unit_names(i) = &
"Simulation/SimulationMain"
    i = i + 1; unit_names(i) = &
"Simulation/SimulationMain/RadCosmo"
    i = i + 1; unit_names(i) = &
"flashUtilities"
    i = i + 1; unit_names(i) = &
"flashUtilities/contiguousConversion"
    i = i + 1; unit_names(i) = &
"flashUtilities/general"
    i = i + 1; unit_names(i) = &
"flashUtilities/interpolation"
    i = i + 1; unit_names(i) = &
"flashUtilities/interpolation/oneDim"
    i = i + 1; unit_names(i) = &
"flashUtilities/jeans_overdens_sink_refinement"
    i = i + 1; unit_names(i) = &
"flashUtilities/nameValueLL"
    i = i + 1; unit_names(i) = &
"flashUtilities/prolong_particles"
    i = i + 1; unit_names(i) = &
"flashUtilities/sorting"
    i = i + 1; unit_names(i) = &
"flashUtilities/sorting/quicksort"
    i = i + 1; unit_names(i) = &
"flashUtilities/system"
    i = i + 1; unit_names(i) = &
"flashUtilities/system/memoryUsage"
    i = i + 1; unit_names(i) = &
"flashUtilities/system/memoryUsage/legacy"
    i = i + 1; unit_names(i) = &
"monitors"
    i = i + 1; unit_names(i) = &
"monitors/Debugger"
    i = i + 1; unit_names(i) = &
"monitors/Logfile"
    i = i + 1; unit_names(i) = &
"monitors/Logfile/LogfileMain"
    i = i + 1; unit_names(i) = &
"monitors/Profiler"
    i = i + 1; unit_names(i) = &
"monitors/Timers"
    i = i + 1; unit_names(i) = &
"monitors/Timers/TimersMain"
    i = i + 1; unit_names(i) = &
"monitors/Timers/TimersMain/MPINative"
    i = i + 1; unit_names(i) = &
"physics"
    i = i + 1; unit_names(i) = &
"physics/Cosmology"
    i = i + 1; unit_names(i) = &
"physics/Cosmology/CosmologyMain"
    i = i + 1; unit_names(i) = &
"physics/Cosmology/CosmologyMain/MatterLambdaKernel"
    i = i + 1; unit_names(i) = &
"physics/Diffuse"
    i = i + 1; unit_names(i) = &
"physics/Eos"
    i = i + 1; unit_names(i) = &
"physics/Eos/EosMain"
    i = i + 1; unit_names(i) = &
"physics/Eos/EosMain/Multigamma"
    i = i + 1; unit_names(i) = &
"physics/Eos/localAPI"
    i = i + 1; unit_names(i) = &
"physics/Gravity"
    i = i + 1; unit_names(i) = &
"physics/Gravity/GravityMain"
    i = i + 1; unit_names(i) = &
"physics/Gravity/GravityMain/Poisson"
    i = i + 1; unit_names(i) = &
"physics/Gravity/GravityMain/Poisson/Multigrid"
    i = i + 1; unit_names(i) = &
"physics/Hydro"
    i = i + 1; unit_names(i) = &
"physics/Hydro/HydroMain"
    i = i + 1; unit_names(i) = &
"physics/Hydro/HydroMain/split"
    i = i + 1; unit_names(i) = &
"physics/Hydro/HydroMain/split/PPM"
    i = i + 1; unit_names(i) = &
"physics/Hydro/HydroMain/split/PPM/PPMKernel"
    i = i + 1; unit_names(i) = &
"physics/ImBound"
    i = i + 1; unit_names(i) = &
"physics/IncompNS"
    i = i + 1; unit_names(i) = &
"physics/RadTrans"
    i = i + 1; unit_names(i) = &
"physics/SolidMechanics"
    i = i + 1; unit_names(i) = &
"physics/TreeCol"
    i = i + 1; unit_names(i) = &
"physics/materialProperties/Conductivity"
    i = i + 1; unit_names(i) = &
"physics/materialProperties/MagneticResistivity"
    i = i + 1; unit_names(i) = &
"physics/materialProperties/MassDiffusivity"
    i = i + 1; unit_names(i) = &
"physics/materialProperties/NSE"
    i = i + 1; unit_names(i) = &
"physics/materialProperties/Opacity"
    i = i + 1; unit_names(i) = &
"physics/materialProperties/Viscosity"
    i = i + 1; unit_names(i) = &
"physics/sourceTerms"
    i = i + 1; unit_names(i) = &
"physics/sourceTerms/Burn"
    i = i + 1; unit_names(i) = &
"physics/sourceTerms/Chem"
    i = i + 1; unit_names(i) = &
"physics/sourceTerms/Chem/ChemMain"
    i = i + 1; unit_names(i) = &
"physics/sourceTerms/ColumnDensity"
    i = i + 1; unit_names(i) = &
"physics/sourceTerms/ColumnDensity/ColumnDensityMain"
    i = i + 1; unit_names(i) = &
"physics/sourceTerms/Cool"
    i = i + 1; unit_names(i) = &
"physics/sourceTerms/Deleptonize"
    i = i + 1; unit_names(i) = &
"physics/sourceTerms/EnergyDeposition"
    i = i + 1; unit_names(i) = &
"physics/sourceTerms/Flame"
    i = i + 1; unit_names(i) = &
"physics/sourceTerms/Heat"
    i = i + 1; unit_names(i) = &
"physics/sourceTerms/Heatexchange"
    i = i + 1; unit_names(i) = &
"physics/sourceTerms/Ionize"
    i = i + 1; unit_names(i) = &
"physics/sourceTerms/Photoionization"
    i = i + 1; unit_names(i) = &
"physics/sourceTerms/Photoionization/PhotoionizationMain"
    i = i + 1; unit_names(i) = &
"physics/sourceTerms/Polytrope"
    i = i + 1; unit_names(i) = &
"physics/sourceTerms/PrimordialChemistry"
    i = i + 1; unit_names(i) = &
"physics/sourceTerms/RadiationField"
    i = i + 1; unit_names(i) = &
"physics/sourceTerms/RadiationField/RadiationFieldMain"
    i = i + 1; unit_names(i) = &
"physics/sourceTerms/Stir"
    i = i + 1; unit_names(i) = &
"physics/sourceTerms/Turb"


    return

  end subroutine setup_getFlashUnits

  subroutine setup_getNumFlashUnits(numUnits)

#include "constants.h"
    implicit none

    integer, intent(out) :: numUnits
    integer, PARAMETER :: NUM_UNITS = 137

    numUnits = NUM_UNITS

    return

  end subroutine setup_getNumFlashUnits

