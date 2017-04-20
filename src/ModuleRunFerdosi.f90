!###########################################################################################################################################################

module ModuleRunFerdosi
    
    implicit none

    !###########################################################################################################################################################

    contains
        
        !###########################################################################################################################################################


        subroutine RunFerdosiPerSire(AllParameters)

            use ModuleParameters
            use Ferdosi

            implicit none

            type(Parameters), intent(inout) :: AllParameters

            integer :: i, nHalfSib, ParentGender
            
            print*, "Running Ferdosi for each sire..."

            do i=1, AllParameters%nSire
                nHalfSib = AllParameters%SireArray(i)%nOffSpring
                ParentGender = 1
                print*, AllParameters%SireArray(i)%ID
                call FerdosiRunner(nHalfSib, AllParameters%nSnp, ParentGender, AllParameters%SireArray(i)%MyPhase, AllParameters%SireArray(i)%OffspringGenotypes, AllParameters%SireArray(i)%OffspringPhase)

            enddo

            
        end subroutine RunFerdosiPerSire

        !###########################################################################################################################################################

 
 
end module ModuleRunFerdosi

!###########################################################################################################################################################
