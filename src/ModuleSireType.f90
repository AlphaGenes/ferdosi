!###########################################################################################################################################################
module ModuleSireType
    use PedigreeModule
    implicit none


    !###########################################################################################################################################################

    type Sire

        double precision :: YieldP1, YieldP2
        double precision :: PhaseAccP1, PhaseAccP2
        type(Individual), pointer :: ind
        

    end type Sire

    !###########################################################################################################################################################

end module ModuleSireType

!###########################################################################################################################################################