!###########################################################################################################################################################
module ModuleSireType

    implicit none


    !###########################################################################################################################################################

    type Sire

    	integer :: ID, PosID
    	integer :: nOffspring
        double precision :: YieldP1, YieldP2
        double precision :: PhaseAccP1, PhaseAccP2
        integer, allocatable, dimension(:) :: OffspringIds
        real, allocatable, dimension(:,:) :: OffspringGenotypes, MyPhase
        real, allocatable, dimension(:,:,:) :: OffspringPhase
        

    end type Sire

    !###########################################################################################################################################################

end module ModuleSireType

!###########################################################################################################################################################