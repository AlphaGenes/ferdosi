!###########################################################################################################################################################

module ModuleRunFerdosi
    
    implicit none

    !###########################################################################################################################################################

    contains
        
        !###########################################################################################################################################################


        subroutine RunFerdosiPerSire(AllParameters, ped)

            use ModuleParameters
            use Ferdosi
            use PedigreeModule
            implicit none

            type(Parameters), intent(inout) :: AllParameters
            type(PedigreeHolder) :: ped
            type(IndividualLinkedListNode),pointer :: ind

            integer :: i
            
            print*, "Running Ferdosi for each sire..."
            ind => ped%sireList%first
            do i=1, ped%sireList%length
                call FerdosiRunner(ind%item, AllParameters%nSnp, AllParameters%overwriteHalfSibPhase)
                ind => ind%next
            enddo

            
        end subroutine RunFerdosiPerSire

        !###########################################################################################################################################################

 
 
end module ModuleRunFerdosi

!###########################################################################################################################################################
